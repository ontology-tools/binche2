from load_chebi import load_ontology, load_chebi, download_chebi
from collections import defaultdict
import json
import xml.etree.ElementTree as ET

def _normalize_iri(value):
    return str(value).strip("<>")


def find_has_role_connections_from_owl(owl_file, has_role_property, deprecated_property, roles_map_json):
    """Parse OWL XML directly to find has_role restrictions (streaming to save memory)."""
    roles_map = defaultdict(set)

    ns = {
        "owl": "http://www.w3.org/2002/07/owl#",
        "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
        "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
    }

    has_role_iri = _normalize_iri(has_role_property)
    deprecated_iri = _normalize_iri(deprecated_property)

    class_tag = f"{{{ns['owl']}}}Class"
    deprecated_tag = f"{{{ns['owl']}}}deprecated"
    subclass_tag = f"{{{ns['rdfs']}}}subClassOf"
    restriction_tag = f"{{{ns['owl']}}}Restriction"
    on_property_tag = f"{{{ns['owl']}}}onProperty"
    some_values_tag = f"{{{ns['owl']}}}someValuesFrom"
    rdf_about = f"{{{ns['rdf']}}}about"
    rdf_resource = f"{{{ns['rdf']}}}resource"

    context = ET.iterparse(owl_file, events=("end",))
    processed = 0

    for event, elem in context:
        if elem.tag != class_tag:
            continue

        cls_iri = elem.get(rdf_about)
        if not cls_iri:
            elem.clear()
            continue

        # Skip deprecated classes if flagged
        deprecated_elem = elem.find(deprecated_tag)
        if deprecated_elem is not None and (deprecated_elem.text or "").strip().lower() == "true":
            elem.clear()
            continue

        # Look for restrictions in rdfs:subClassOf
        for sub in elem.findall(subclass_tag):
            restriction = None
            if len(sub):
                for child in sub:
                    if child.tag == restriction_tag:
                        restriction = child
                        break
            if restriction is None:
                continue

            on_prop = restriction.find(on_property_tag)
            if on_prop is None:
                continue
            prop_iri = on_prop.get(rdf_resource)
            if not prop_iri or _normalize_iri(prop_iri) != has_role_iri:
                continue

            filler = restriction.find(some_values_tag)
            if filler is None:
                continue
            filler_iri = filler.get(rdf_resource)
            if filler_iri:
                roles_map[cls_iri].add(filler_iri)

        processed += 1
        if processed % 50000 == 0:
            print(f"Parsed {processed} classes (OWL XML fallback)...")

        elem.clear()

    roles_map_serializable = {k: list(v) for k, v in roles_map.items()}
    with open(roles_map_json, "w") as f:
        json.dump(roles_map_serializable, f, indent=2)

    print(f"Saved roles map with {len(roles_map_serializable)} classes to {roles_map_json} (OWL XML fallback).")
    return roles_map

def create_leaves_to_all_roles_map(roles_map_json, leaves_to_all_parents_json, leaves_to_all_roles_json, parent_map_json):
    """Create a map from each leaf class to all roles associated with it via its ancestors and directly.

    Includes role ancestors (via parent map) so if a leaf has role X, then all ancestors of X are also included.
    """

    with open(roles_map_json, "r") as f:
        roles_map = json.load(f)

    with open(leaves_to_all_parents_json, "r") as f:
        leaves_to_all_parents = json.load(f)

    with open(parent_map_json, "r") as f:
        parent_map = json.load(f) # All classes to direct parents
    
    role_ancestor_cache = {}

    def get_role_ancestors(role_iri):
        if role_iri in role_ancestor_cache:
            return role_ancestor_cache[role_iri]

        ancestors = set()
        stack = [role_iri]
        visited = set()

        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)

            for parent in parent_map.get(current, []):
                if parent not in ancestors:
                    ancestors.add(parent)
                    stack.append(parent)

        role_ancestor_cache[role_iri] = ancestors
        return ancestors

    leaves_to_all_roles = {}

    for leaf_iri, all_parents in leaves_to_all_parents.items():
        all_roles = set()
        # Include direct roles on the leaf itself, if any
        if leaf_iri in roles_map:
            all_roles.update(roles_map[leaf_iri])
        for parent in all_parents:
            parent = parent.strip()
            if parent in roles_map:
                all_roles.update(roles_map[parent])

        # Include ancestors of each role
        expanded_roles = set(all_roles)
        for role in list(all_roles):
            expanded_roles.update(get_role_ancestors(role))

        leaves_to_all_roles[leaf_iri] = list(expanded_roles)

    print(f"Built leaf to all roles map with {len(leaves_to_all_roles)} leaf classes.")

    # Calc the number of empty role lists (i.e. leaf classes with no roles via ancestors)
    num_empty_roles = sum(1 for roles in leaves_to_all_roles.values() if not roles)
    print(f"Number of leaf classes with no roles via ancestors: {num_empty_roles}")

    with open(leaves_to_all_roles_json, "w") as f:
        json.dump(leaves_to_all_roles, f, indent=2)

def create_class_to_all_roles_map(roles_map_json, parent_map_json, class_to_all_roles_json):
    """Create a map from each class to all roles associated with it (direct + inherited from ancestors).

    Includes direct roles on the class, direct roles on ancestor classes,
    and all descendants of those roles (role hierarchy expansion).
    """

    with open(roles_map_json, "r") as f:
        roles_map = json.load(f)

    with open(parent_map_json, "r") as f:
        parent_map = json.load(f)
    
    # Build ancestor cache for classes
    class_ancestor_cache = {}
    
    def get_class_ancestors(class_iri):
        if class_iri in class_ancestor_cache:
            return class_ancestor_cache[class_iri]

        ancestors = set()
        stack = [class_iri]
        visited = set()

        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)

            for parent in parent_map.get(current, []):
                if parent not in ancestors:
                    ancestors.add(parent)
                    stack.append(parent)

        class_ancestor_cache[class_iri] = ancestors
        return ancestors
    
    # Build role descendant cache
    role_descendant_cache = {}

    # Build child map for quick role descendant lookup
    child_map = defaultdict(list)
    for child, parents in parent_map.items():
        for parent in parents:
            child_map[parent].append(child)

    def get_role_descendants(role_iri):
        if role_iri in role_descendant_cache:
            return role_descendant_cache[role_iri]

        descendants = set()
        stack = [role_iri]
        visited = set()

        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)

            for child in child_map.get(current, []):
                if child not in descendants:
                    descendants.add(child)
                    stack.append(child)

        role_descendant_cache[role_iri] = descendants
        return descendants

    # Collect all unique classes (from roles_map and parent_map)
    all_classes = set(roles_map.keys()) | set(parent_map.keys())
    
    class_to_all_roles = {}

    for class_iri in all_classes:
        all_roles = set()
        
        # Include direct roles on the class itself
        if class_iri in roles_map:
            all_roles.update(roles_map[class_iri])
        
        # Include roles from all ancestors
        for ancestor in get_class_ancestors(class_iri):
            if ancestor in roles_map:
                all_roles.update(roles_map[ancestor])

        # Include descendants of each role
        expanded_roles = set(all_roles)
        for role in list(all_roles):
            expanded_roles.update(get_role_descendants(role))

        if expanded_roles:  # Only include classes that have at least one role
            class_to_all_roles[class_iri] = list(expanded_roles)

    print(f"Built class to all roles map with {len(class_to_all_roles)} classes.")

    with open(class_to_all_roles_json, "w") as f:
        json.dump(class_to_all_roles, f, indent=2)

def create_roles_to_all_leaves_map(leaves_to_all_roles_json, roles_to_all_leaves_json):
    """Create a map from each role to all leaf classes associated with it via descendants."""

    with open(leaves_to_all_roles_json, "r") as f:
        leaves_to_all_roles = json.load(f)

    roles_to_all_leaves = defaultdict(set)

    for leaf_iri, all_roles in leaves_to_all_roles.items():
        for role in all_roles:
            roles_to_all_leaves[role].add(leaf_iri)

    # Convert sets to lists for JSON serialization
    roles_to_all_leaves_json_serializable = {
        role: list(leaves)
        for role, leaves in roles_to_all_leaves.items()
    }

    print(f"Built roles to all leaves map with {len(roles_to_all_leaves_json_serializable)} roles.")

    with open(roles_to_all_leaves_json, "w") as f:
        json.dump(roles_to_all_leaves_json_serializable, f, indent=2)

if __name__ == "__main__":

    task = "build class to all roles map" # Options: "find has_role connections", "build leaf to all roles map", "build class to all roles map"
    # Order of tasks:
    # task "find has_role connections" will parse the OWL XML directly to find has_role restrictions and build the roles map. The ouput file is then used in the next step.
    # task "build leaf to all roles map" will use the roles map and the leaf->all parents map to build a leaf->all roles map (via ancestors). The reverse map from roles->all leaves is also be built here.
    # task "build class to all roles map" will use the roles map and the parent map to build a class->all roles map (via ancestors) for ALL classes.

    has_role_property = "http://purl.obolibrary.org/obo/RO_0000087"
    deprecated_property = "http://www.w3.org/2002/07/owl#deprecated"

    roles_map_json = "data/class_to_direct_roles_map.json" # output file for roles map
    leaves_to_all_parents_json = "data/removed_leaf_classes_to_ALL_parents_map.json"
    leaves_to_all_roles_json = "data/removed_leaf_classes_to_ALL_roles_map.json" # output file for leaves to all roles map
    parent_map_json = "data/chebi_parent_map.json"
    roles_to_all_leaves_json = "data/roles_to_leaves_map.json" # output file for roles to all leaves map
    class_to_all_roles_json = "data/class_to_all_roles_map.json" # output file for class to all roles map

    if task == "find has_role connections":

        print("="*80 + "\n")
        
        # Run OWL XML parsing directly (more reliable for has_role restrictions)
        owl_file = download_chebi()
        roles_map = find_has_role_connections_from_owl(owl_file, has_role_property, deprecated_property, roles_map_json)
        
        # Summary statistics
        print(f"\n{'='*80}")
        print("SUMMARY STATISTICS")
        print(f"{'='*80}")
        print(f"Total classes with roles: {len(roles_map)}")
        
        if len(roles_map) > 0:
            total_role_connections = sum(len(roles) for roles in roles_map.values())
            print(f"Total role connections: {total_role_connections}")
            avg_roles = total_role_connections / len(roles_map)
            print(f"Average roles per class: {avg_roles:.2f}")
            
            max_roles_class = max(roles_map.items(), key=lambda x: len(x[1]))
            print(f"Class with most roles: {max_roles_class[0]} ({len(max_roles_class[1])} roles)")
        
        # Build leaf -> all roles map using leaf->parents JSON
        create_leaves_to_all_roles_map(roles_map_json, leaves_to_all_parents_json, leaves_to_all_roles_json, parent_map_json)

        # Verification test for function find_has_role_connections_from_owl
        test_class_iri = "http://purl.obolibrary.org/obo/CHEBI_10002"
        
        if test_class_iri in roles_map:
            print(f"\n{'='*80}")
            print(f"VERIFICATION TEST: CHEBI_10002")
            print(f"{'='*80}")
            print(f"Found {len(roles_map[test_class_iri])} roles:")
            for role in roles_map[test_class_iri]:
                print(f"  - {role}")
            
            expected = {
                "http://purl.obolibrary.org/obo/CHEBI_140378",
                "http://purl.obolibrary.org/obo/CHEBI_35620",
                "http://purl.obolibrary.org/obo/CHEBI_35674",
                "http://purl.obolibrary.org/obo/CHEBI_38231"
            }
            
            found = set(roles_map[test_class_iri])
            
            if found == expected:
                print(f"\n✅ SUCCESS! All expected roles found.")
            elif found.issuperset(expected):
                print(f"\n✅ SUCCESS! All expected roles found (plus {len(found - expected)} more).")
            else:
                print(f"\n⚠️  WARNING: Missing roles!")
                print(f"Expected but not found: {expected - found}")
                print(f"Found but not expected: {found - expected}")
        else:
            print(f"\n❌ ERROR: CHEBI_10002 not found in roles_map!")
        
    elif task == "build leaf to all roles map":
        
        create_leaves_to_all_roles_map(roles_map_json, leaves_to_all_parents_json, leaves_to_all_roles_json, parent_map_json)
        create_roles_to_all_leaves_map(leaves_to_all_roles_json, roles_to_all_leaves_json)

    elif task == "build class to all roles map":
        
        create_class_to_all_roles_map(roles_map_json, parent_map_json, class_to_all_roles_json)

    else:
        print(f"Unknown task: {task}. Please select a valid task.")