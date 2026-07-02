import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from flask import Flask, render_template, request, redirect, url_for, session, Response
from fishers_calculations import run_enrichment_analysis, run_enrichment_analysis_plain_enrich_pruning_strategy
from wikidata.narrow_background_fishers import run_narrow_background_enrichment_analysis, run_narrow_background_enrichment_analysis_plain_enrich_pruning_strategy
from weighted_calculations import run_weighted_enrichment_analysis, run_weighted_enrichment_analysis_plain_enrich_pruning_strategy, run_weighted_narrow_background_enrichment_analysis, run_weighted_narrow_background_enrichment_analysis_plain_enrich_pruning_strategy, auto_scale_weights
from visualitations_and_pruning import graph_to_cytospace_json
import re
import csv
from io import StringIO
import uuid
import time
import glob
import uuid
import requests
from rdkit import Chem
from rdkit.Chem import inchi


app = Flask(__name__)
app.secret_key = 'your_secret_key'  # Replace with a secure secret key

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
LOCAL_LOOKUP_FILE = os.path.join(BASE_DIR, 'data', 'removed_leaf_classes_with_inchikeys.csv')

# Maps a "background" form value to its narrow-background leaves JSON file.
NARROW_BACKGROUND_LEAVES_JSON = {
    'human': 'data/human_entities_leaves.json',
    'arabidopsis_thaliana': 'data/arabidopsis_thaliana_leaves.json',
    'endogenous_human': 'data/recon3d_leaves.json',
}


def _clean_lookup_value(value):
    if value is None:
        return ''
    return str(value).strip().strip('"')


def _canonical_smiles(smiles):
    """RDKit-canonical form of a SMILES string, or None if it can't be parsed.

    ChEBI's asserted SMILES aren't guaranteed to be in any particular canonical
    form, so comparing raw strings misses matches between differently-written
    SMILES for the same molecule. Canonicalizing both the lookup table and any
    incoming SMILES through RDKit makes exact-match comparison toolkit-consistent.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return None
    return Chem.MolToSmiles(mol) if mol is not None else None


def _normalize_chebi_id(raw_value):
    value = _clean_lookup_value(raw_value)
    if not value:
        return ''

    if value.startswith('http://purl.obolibrary.org/obo/CHEBI_'):
        value = value.rsplit('/', 1)[-1]

    if value.startswith('CHEBI:'):
        return value.replace(':', '_', 1)

    if value.startswith('CHEBI_'):
        return value

    if value.isdigit():
        return f'CHEBI_{value}'

    return value


def _chebi_sort_key(chebi_id):
    """Numeric sort key for 'CHEBI_12345' ids, so collision tie-breaks are by
    ChEBI ID number rather than CSV row order (which is incidental).
    """
    try:
        return int(chebi_id.rsplit('_', 1)[-1])
    except (ValueError, AttributeError):
        return float('inf')


def _resolve_lookup_collisions(candidates, label):
    """Pick a deterministic winner (lowest ChEBI ID) for each key asserted by more
    than one ChEBI term, and return both the winner map and a {key: [all_ids]} map
    of only the keys that actually collided, so callers can surface the ones that
    were dropped. Logs a single summary line rather than one line per collision,
    since this table has thousands of them; per-match ambiguity is surfaced to
    end users in the Processing Summary instead (see convert_smiles_to_chebi).
    """
    resolved = {}
    collisions = {}
    for key, chebi_ids in candidates.items():
        ordered = sorted(chebi_ids, key=_chebi_sort_key)
        resolved[key] = ordered[0]
        if len(ordered) > 1:
            collisions[key] = ordered
    if collisions:
        print(
            f"Warning: {len(collisions)} distinct {label} values in {LOCAL_LOOKUP_FILE} "
            f"are each asserted by more than one ChEBI term; using the lowest ChEBI ID "
            f"as a deterministic tie-break for each."
        )
    return resolved, collisions


def _load_local_smiles_and_inchikey_maps():
    smiles_candidates = {}
    inchikey_candidates = {}

    with open(LOCAL_LOOKUP_FILE, 'r', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []

        smiles_column = None
        inchikey_column = None
        iri_column = None

        for candidate in ('SMILES', 'smiles'):
            if candidate in fieldnames:
                smiles_column = candidate
                break

        for candidate in ('InChIKey', 'InChIkey', 'inchikey', 'InChIKEY'):
            if candidate in fieldnames:
                inchikey_column = candidate
                break

        for candidate in ('IRI', 'iri'):
            if candidate in fieldnames:
                iri_column = candidate
                break

        if smiles_column is None or iri_column is None:
            raise KeyError(
                f"{LOCAL_LOOKUP_FILE} must contain SMILES and IRI columns to build local lookup maps"
            )

        for row in reader:
            chebi_id = _normalize_chebi_id(row.get(iri_column))
            smiles = _clean_lookup_value(row.get(smiles_column))
            inchikey_value = _clean_lookup_value(row.get(inchikey_column)) if inchikey_column else ''

            if smiles:
                ids = smiles_candidates.setdefault(smiles, [])
                if chebi_id not in ids:
                    ids.append(chebi_id)
            if inchikey_value:
                inchikey_key = inchikey_value.upper()
                ids = inchikey_candidates.setdefault(inchikey_key, [])
                if chebi_id not in ids:
                    ids.append(chebi_id)

    smiles_to_chebi, smiles_collisions = _resolve_lookup_collisions(smiles_candidates, 'SMILES')
    inchikey_to_chebi, inchikey_collisions = _resolve_lookup_collisions(inchikey_candidates, 'InChIKey')

    return smiles_to_chebi, inchikey_to_chebi, smiles_collisions, inchikey_collisions


(LOCAL_SMILES_TO_CHEBI, LOCAL_INCHIKEY_TO_CHEBI,
 LOCAL_SMILES_COLLISIONS, LOCAL_INCHIKEY_COLLISIONS) = _load_local_smiles_and_inchikey_maps()

def cleanup_old_graph_files(max_age_hours=24):
    """Remove graph files older than max_age_hours"""
    cutoff = time.time() - (max_age_hours * 3600)
    for filepath in glob.glob('website/static/data/graph_*.json'):
        try:
            if os.path.getmtime(filepath) < cutoff:
                os.remove(filepath)
        except OSError:
            pass


# Load last update time
def get_data_version():
    try:
        with open(os.path.join(BASE_DIR, 'data_version.txt'), 'r') as f:
            return f.read().strip()
    except FileNotFoundError:
        return 'unknown'

@app.context_processor
def inject_data_version():
    return dict(data_version=get_data_version())

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/submission', methods=['GET', 'POST'])
def submission():
    if request.method == 'POST':
        study_set = request.form.get('study_set')
        session['study_set'] = study_set # Store the study set in session
        classification = request.form.get('classification')
        session['classification'] = classification # Store classification in session
        smiles_option = request.form.get('smiles_option')
        session['smiles_option'] = smiles_option # Store smiles option in session
        background = request.form.get('background') 
        session['background'] = background # Store background in session
        # Store user's preference for expanding the human background (checkbox)
        expand_background = bool(request.form.get('expand_background'))
        session['expand_background'] = expand_background
        
        return render_template('submission.html', user_study_set=study_set)
    
    return render_template('submission.html', user_study_set=None)

# Used in run_analysis route to parse user input
def parse_studyset(studyset: str):
    # Remove surrounding quotes if present
    studyset = studyset.strip()

    studyset_list = []
    weights_dict = {}
    unresolved_smiles = []
    ambiguous_smiles_matches = []

    if not studyset:
        return studyset_list, weights_dict, unresolved_smiles, ambiguous_smiles_matches

    def normalize_id(raw_id: str) -> str:
        value = raw_id.strip().replace('"', "")
        if value.startswith("http://") or value.startswith("https://"):
            return value
        if value.startswith("CHEBI:"):
            return value.replace(":", "_")
        return value.replace(":", "_")
    
    def is_smiles(value: str) -> bool:
        """Detect if a string is likely a SMILES string."""
        value = value.strip()
        # Not a SMILES if it starts with CHEBI: or http
        if value.startswith("CHEBI:") or value.startswith("http://") or value.startswith("https://"):
            return False
        # SMILES typically contain lowercase letters, parentheses, or specific chars
        smiles_chars = set('cCnNoOpPsSFfIiBbr[]()=#@+-\\/')
        return any(c in smiles_chars for c in value)

    def record_ambiguous(smiles, ambiguous_match):
        if ambiguous_match is None:
            return
        chosen, all_ids = ambiguous_match
        ambiguous_smiles_matches.append({
            'smiles': smiles,
            'chosen': chosen,
            'alternatives': [cid for cid in all_ids if cid != chosen],
        })

    # Split by lines first to support optional weights per line
    for line in studyset.splitlines():
        line = line.strip()
        if not line:
            continue

        parts = [p for p in re.split(r'[\s,]+', line) if p]

        if len(parts) >= 2:
            try:
                weight = float(parts[1])
                # Check if first part is SMILES
                if is_smiles(parts[0]):
                    chebi_ids, was_resolved, ambiguous_match = convert_smiles_to_chebi(parts[0])
                    if not was_resolved:
                        unresolved_smiles.append(parts[0])
                    record_ambiguous(parts[0], ambiguous_match)
                    # Apply the same weight to all resulting ChEBI IDs
                    for chebi_id in chebi_ids:
                        class_id = normalize_id(chebi_id)
                        studyset_list.append(class_id)
                        weights_dict[class_id] = weight
                else:
                    class_id = normalize_id(parts[0])
                    studyset_list.append(class_id)
                    weights_dict[class_id] = weight
                continue
            except ValueError:
                pass

        # Fallback: treat all parts as IDs without weights
        for part in parts:
            if is_smiles(part):
                chebi_ids, was_resolved, ambiguous_match = convert_smiles_to_chebi(part)
                if not was_resolved:
                    unresolved_smiles.append(part)
                record_ambiguous(part, ambiguous_match)
                for chebi_id in chebi_ids:
                    class_id = normalize_id(chebi_id)
                    studyset_list.append(class_id)
            else:
                class_id = normalize_id(part)
                studyset_list.append(class_id)

    return studyset_list, weights_dict, unresolved_smiles, ambiguous_smiles_matches

def convert_smiles_to_chebi(smiles_string):
    """Convert a single SMILES string to ChEBI IDs.

    Returns (chebi_ids_list, was_resolved, ambiguous_match). ambiguous_match is
    None unless the matched SMILES/InChIKey is asserted by more than one ChEBI
    term in the local lookup table, in which case it's (chosen_chebi_id,
    all_chebi_ids) so the caller can surface the ambiguity to the user.
    """
    chebi_ids = []
    was_resolved = False
    ambiguous_match = None
    cleaned_smiles = _clean_lookup_value(smiles_string)
    try:
        mol = Chem.MolFromSmiles(cleaned_smiles)
    except Exception as error:
        print(f"Warning: failed to parse SMILES {cleaned_smiles}: {error}")
        mol = None
    canonical_smiles = Chem.MolToSmiles(mol) if mol is not None else None

    # First try a direct SMILES lookup against the local leaf-class table, comparing
    # canonical forms so the match doesn't depend on how either SMILES was written.
    lookup_key = canonical_smiles or cleaned_smiles
    local_chebi_id = LOCAL_SMILES_TO_CHEBI.get(lookup_key)
    if local_chebi_id:
        chebi_ids.append(local_chebi_id.replace('_', ':', 1))
        was_resolved = True
        if lookup_key in LOCAL_SMILES_COLLISIONS:
            ambiguous_match = (local_chebi_id, LOCAL_SMILES_COLLISIONS[lookup_key])
        print(f"Found local ChEBI ID from exact SMILES match: {local_chebi_id} for SMILES {cleaned_smiles}")
        return chebi_ids, was_resolved, ambiguous_match

    # If no direct SMILES match exists, try InChIKey -> ChEBI using RDKit to
    # compute the InChIKey for the submitted SMILES.
    try:
        if mol is not None:
            user_inchikey = inchi.MolToInchiKey(mol).upper()
            local_chebi_id = LOCAL_INCHIKEY_TO_CHEBI.get(user_inchikey)
            if local_chebi_id:
                chebi_ids.append(local_chebi_id.replace('_', ':', 1))
                was_resolved = True
                if user_inchikey in LOCAL_INCHIKEY_COLLISIONS:
                    ambiguous_match = (local_chebi_id, LOCAL_INCHIKEY_COLLISIONS[user_inchikey])
                print(
                    f"Found local ChEBI ID from InChIKey match: {local_chebi_id} "
                    f"for SMILES {cleaned_smiles} (InChIKey {user_inchikey})"
                )
                return chebi_ids, was_resolved, ambiguous_match
    except Exception as error:
        print(f"Warning: failed to compute InChIKey for SMILES {cleaned_smiles}: {error}")

    # Get details from ChEBI lookup to check for a direct match to a ChEBI ID
    response = requests.post("https://chebifier.hastingslab.org/api/details", json={
        "type": "type",
        "smiles": cleaned_smiles,
        "selectedModels": {
            "ChEBI Lookup": True
        }
    })
    
    lookup_infotext = response.json().get("models", {}).get("ChEBI Lookup", {}).get("highlights", [])
    
    # If the lookup highlights contain a ChEBI ID, use that for classification instead of the SMILES string. 
    if lookup_infotext and "CHEBI:" in lookup_infotext[0][1]:
        chebi_id = lookup_infotext[0][1].split("CHEBI:")[1].split()[0].rstrip('.')
        chebi_ids.append(f"CHEBI:{chebi_id}")
        was_resolved = True
        print(f"Found ChEBI ID from lookup: CHEBI:{chebi_id} for SMILES {cleaned_smiles}")
    else:
        # Get direct parents from classification
        smiles_option = session.get('smiles_option')

        if smiles_option == 'use_parents':
            print(f"No direct ChEBI ID found from lookup for SMILES {cleaned_smiles}, attempting classification...")
            response = requests.post("https://chebifier.hastingslab.org/api/classify", json={
                "smiles": cleaned_smiles,
                "ontology": False,
                "selectedModels": {
                    "ELECTRA (ChEBI50-3STAR)": True,
                }
            })
            
            direct_parents = response.json().get("direct_parents")
            if direct_parents:
                # Extract ChEBI IDs from all parent lists
                for parent_list in direct_parents:
                    if parent_list is not None:
                        parent_ids = [f"CHEBI:{parent[0]}" for parent in parent_list]
                        chebi_ids.extend(parent_ids)
                        # Print the parent IDs found
                        print(f"Found direct parent ChEBI IDs from classification for SMILES {cleaned_smiles}: {parent_ids}")
                    else:
                        print(f"No parents found in one of the classification results for SMILES {cleaned_smiles}")
                        # print response content for debugging
                        print(f"Classification response content: {response.content}")
                if chebi_ids:
                    was_resolved = True
                
                print(f"Found {len(chebi_ids)} ChEBI IDs from classification for SMILES {cleaned_smiles}")

        else:
            print(f"No direct ChEBI ID found from lookup for SMILES {cleaned_smiles}, excluding from analysis.")


    return chebi_ids, was_resolved, ambiguous_match
 
def map_p_value_correction_method(method_name):
    if method_name == 'bonferroni':
        return True, False
    elif method_name == 'benjamini_hochberg':
        return False, True
    else:
        return False, False # No correction method selected

@app.route('/run_analysis', methods=['GET', 'POST'])
# option to choose pruning methods
def run_analysis():
    # Cleanup old graph files
    cleanup_old_graph_files(max_age_hours=24)
    
    # Generate or retrieve session ID
    if 'session_id' not in session:
        session['session_id'] = str(uuid.uuid4())

    raw_studyset = session.get('study_set')
    if not raw_studyset:
        return redirect(url_for('submission'))
    # Convert multi-line or comma-separated input into a list
    studyset_list, weights_dict, unresolved_smiles, ambiguous_smiles_matches = parse_studyset(raw_studyset)
    
    # Auto-scale weights if present (only scales up if max < 1000)
    # if weights_dict:
    #     weights_dict = auto_scale_weights(weights_dict, target_max=1000)
    
    session['weights_dict'] = weights_dict

    # Get classification from form (allow changing it during re-run)
    classification = request.form.get('classification')
    if classification:
        session['classification'] = classification
    
    # Get background from form (allow changing it during re-run)
    background_from_form = request.form.get('background')
    if background_from_form:
        session['background'] = background_from_form

    # Get expand_background from form (allow changing it during re-run)
    if request.method == 'POST':
        session['expand_background'] = bool(request.form.get('expand_background'))

    # Keep the correction selection stable across re-runs.
    previous_correction = session.get('correction_method', {
        'bonferroni_correct': False,
        'benjamini_hochberg_correct': False
    })
    method = request.form.get("p_value_correction_method")
    if method:
        bonferroni_correct, benjamini_hochberg_correct = map_p_value_correction_method(method)
    else:
        bonferroni_correct = previous_correction.get('bonferroni_correct', False)
        benjamini_hochberg_correct = previous_correction.get('benjamini_hochberg_correct', False)

    session['correction_method'] = {
        'bonferroni_correct': bonferroni_correct,
        'benjamini_hochberg_correct': benjamini_hochberg_correct
    }

    background = session.get('background')
    looping_prune_method = request.form.get('looping_prune_method')


    if looping_prune_method == 'no_loop_prune':
        # User chose no looping pruning, proceed to other pruning options
        pass
    elif looping_prune_method == 'plain_enrich':
        if background == 'full':
            # Use weighted analysis if weights are present, otherwise use standard analysis
            if weights_dict:
                results, pruned_G = run_weighted_enrichment_analysis(weights_dict,
                                                    classification=session.get('classification'))
            else:
                results, pruned_G = run_enrichment_analysis_plain_enrich_pruning_strategy(studyset_list,
                                                    classification=session.get('classification'))
            # Save JSON representation of pruned_G in session for graph visualization
            graph_json_file = f'website/static/data/graph_{session["session_id"]}.json'
            graph_to_cytospace_json(pruned_G, graph_json_file, results)
            session['graph_file'] = f'graph_{session["session_id"]}.json'

            session['pruning'] = {
                'method': 'plain_enrich'
            }

            session['unresolved_smiles'] = unresolved_smiles
            session['ambiguous_smiles_matches'] = ambiguous_smiles_matches

            return render_template('results.html', results=results, graph_json_file=graph_json_file, unresolved_smiles=unresolved_smiles, ambiguous_smiles_matches=ambiguous_smiles_matches, smiles_option=session.get('smiles_option'), expand_background=session.get('expand_background', True), background=session.get('background', 'full'))

        elif background in NARROW_BACKGROUND_LEAVES_JSON:
            if weights_dict:
                results, pruned_G, leaves_to_expand_background, parents_to_expand_background = run_weighted_narrow_background_enrichment_analysis_plain_enrich_pruning_strategy(
                    weights_dict,
                    classification=session.get('classification'),
                    narrow_background_leaves_json=NARROW_BACKGROUND_LEAVES_JSON[background],
                    expand_background=session.get('expand_background', True)
                )
            else:
                results, pruned_G, leaves_to_expand_background, parents_to_expand_background = run_narrow_background_enrichment_analysis_plain_enrich_pruning_strategy(
                    studyset_list,
                    classification=session.get('classification'),
                    narrow_background_leaves_json=NARROW_BACKGROUND_LEAVES_JSON[background],
                    expand_background=session.get('expand_background', True)
                )

            graph_json_file = f'website/static/data/graph_{session["session_id"]}.json'
            graph_to_cytospace_json(pruned_G, graph_json_file, results)
            session['graph_file'] = f'graph_{session["session_id"]}.json'

            session['pruning'] = {
                'method': 'plain_enrich'
            }

            session['unresolved_smiles'] = unresolved_smiles
            session['ambiguous_smiles_matches'] = ambiguous_smiles_matches

            return render_template('results.html', results=results, graph_json_file=graph_json_file, unresolved_smiles=unresolved_smiles, ambiguous_smiles_matches=ambiguous_smiles_matches, smiles_option=session.get('smiles_option'), leaves_to_expand_background=leaves_to_expand_background, parents_to_expand_background=parents_to_expand_background, expand_background=session.get('expand_background', True), background=session.get('background', 'full'))

    # PRUNING OPTIONS
    root_children_prune = request.form.get('root_children_prune') == 'true'
    levels = int(request.form.get('levels', 2))

    linear_branch_prune = request.form.get('linear_branch_prune') == 'true'
    linear_branch_n = int(request.form.get('linear_branch_n', 2))

    high_p_value_prune = request.form.get('high_p_value_prune') == 'true'
    p_value_threshold = float(request.form.get('p_value_threshold', 0.05))

    zero_degree_prune = request.form.get('zero_degree_prune') == 'true'
    
    background = session.get('background')
    
    # Initialize expanded leaves/parents tracking (for narrow backgrounds)
    leaves_to_expand_background = set()
    parents_to_expand_background = set()

    # Use weighted analysis if weights are present, otherwise use standard analysis
    if weights_dict:
        if background in NARROW_BACKGROUND_LEAVES_JSON:
            results, pruned_G, leaves_to_expand_background, parents_to_expand_background = run_weighted_narrow_background_enrichment_analysis(
                                                weights_dict,
                                                levels=levels,
                                                n=linear_branch_n,
                                                p_value_threshold=p_value_threshold,
                                                classification=session.get('classification'),
                                                root_children_prune=root_children_prune,
                                                linear_branch_prune=linear_branch_prune,
                                                high_p_value_prune=high_p_value_prune,
                                                zero_degree_prune=zero_degree_prune,
                                                bonferroni_correct=bonferroni_correct,
                                                benjamini_hochberg_correct=benjamini_hochberg_correct,
                                                narrow_background_leaves_json=NARROW_BACKGROUND_LEAVES_JSON[background],
                                                expand_background=session.get('expand_background', True))
        else:
            results, pruned_G = run_weighted_enrichment_analysis(weights_dict,
                                            levels=levels,
                                            n=linear_branch_n,
                                            p_value_threshold=p_value_threshold,
                                            classification=session.get('classification'),
                                            root_children_prune=root_children_prune,
                                            linear_branch_prune=linear_branch_prune,
                                            high_p_value_prune=high_p_value_prune,
                                            zero_degree_prune=zero_degree_prune,
                                            bonferroni_correct=bonferroni_correct,
                                            benjamini_hochberg_correct=benjamini_hochberg_correct)
    else:
        if background in NARROW_BACKGROUND_LEAVES_JSON:
            results, pruned_G, leaves_to_expand_background, parents_to_expand_background = run_narrow_background_enrichment_analysis(studyset_list,
                                              bonferroni_correct=bonferroni_correct,
                                              benjamini_hochberg_correct=benjamini_hochberg_correct,
                                               root_children_prune=root_children_prune,levels=levels,
                                               linear_branch_prune=linear_branch_prune, n=linear_branch_n,
                                               high_p_value_prune=high_p_value_prune, p_value_threshold=p_value_threshold,
                                               zero_degree_prune=zero_degree_prune,
                                               classification=session.get('classification'),
                                               narrow_background_leaves_json=NARROW_BACKGROUND_LEAVES_JSON[background],
                                               expand_background=session.get('expand_background', True))
        else:
            results, pruned_G = run_enrichment_analysis(studyset_list,
                                              bonferroni_correct=bonferroni_correct,
                                              benjamini_hochberg_correct=benjamini_hochberg_correct,
                                               root_children_prune=root_children_prune,levels=levels,
                                               linear_branch_prune=linear_branch_prune, n=linear_branch_n,
                                               high_p_value_prune=high_p_value_prune, p_value_threshold=p_value_threshold,
                                               zero_degree_prune=zero_degree_prune,
                                               classification=session.get('classification'))
                                       
    # Save JSON representation of pruned_G in session for graph visualization
    graph_json_file = f'website/static/data/graph_{session["session_id"]}.json'
    graph_to_cytospace_json(pruned_G, graph_json_file, results)
    session['graph_file'] = f'graph_{session["session_id"]}.json'

    session['pruning'] = {
        'method': 'custom',
        'root_children_prune': root_children_prune,
        'levels': levels,
        'linear_branch_prune': linear_branch_prune,
        'linear_branch_n': linear_branch_n,
        'high_p_value_prune': high_p_value_prune,
        'p_value_threshold': p_value_threshold,
        'zero_degree_prune': zero_degree_prune
    }

    session['unresolved_smiles'] = unresolved_smiles
    session['ambiguous_smiles_matches'] = ambiguous_smiles_matches

    return render_template('results.html', results=results, graph_json_file=graph_json_file, unresolved_smiles=unresolved_smiles, ambiguous_smiles_matches=ambiguous_smiles_matches, smiles_option=session.get('smiles_option'), leaves_to_expand_background=leaves_to_expand_background, parents_to_expand_background=parents_to_expand_background, expand_background=session.get('expand_background', True), background=session.get('background', 'full'))

@app.route('/graph')
def graph():
    session_id = session.get('session_id')
    graph_file = f'graph_{session_id}.json' if session_id else 'graph.json'
    pruning = session.get('pruning', {})
    correction_method = session.get('correction_method', {})
    classification = session.get('classification', 'structural')
    background = session.get('background', 'full')
    return render_template('graph.html',
                         graph_file=graph_file,
                         pruning=pruning,
                         correction_method=correction_method,
                         classification=classification,
                         background=background,
                         smiles_option=session.get('smiles_option'),
                         expand_background=session.get('expand_background', True))



if __name__ == '__main__':
    app.run(debug=True)
