import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from flask import Flask, render_template, request, redirect, url_for, session, Response
from fishers_calculations import run_enrichment_analysis, run_enrichment_analysis_plain_enrich_pruning_strategy
from weighted_calculations import run_weighted_enrichment_analysis, run_weighted_enrichment_analysis_plain_enrich_pruning_strategy, auto_scale_weights
from visualitations_and_pruning import graph_to_cytospace_json
import re
import csv
from io import StringIO
import uuid
import time
import glob
import uuid
import requests


app = Flask(__name__)
app.secret_key = 'your_secret_key'  # Replace with a secure secret key

def cleanup_old_graph_files(max_age_hours=24):
    """Remove graph files older than max_age_hours"""
    cutoff = time.time() - (max_age_hours * 3600)
    for filepath in glob.glob('website/static/data/graph_*.json'):
        try:
            if os.path.getmtime(filepath) < cutoff:
                os.remove(filepath)
        except OSError:
            pass


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
        return render_template('submission.html', user_study_set=study_set)
    
    return render_template('submission.html', user_study_set=None)

# Used in run_analysis route to parse user input
def parse_studyset(studyset: str):
    # Remove surrounding quotes if present
    studyset = studyset.strip()

    studyset_list = []
    weights_dict = {}

    if not studyset:
        return studyset_list, weights_dict

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
                    chebi_ids = convert_smiles_to_chebi(parts[0])
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
                chebi_ids = convert_smiles_to_chebi(part)
                for chebi_id in chebi_ids:
                    class_id = normalize_id(chebi_id)
                    studyset_list.append(class_id)
            else:
                class_id = normalize_id(part)
                studyset_list.append(class_id)

    return studyset_list, weights_dict

def convert_smiles_to_chebi(smiles_string):
    """Convert a single SMILES string to ChEBI IDs. Returns a list of ChEBI IDs."""
    chebi_ids = []

    # Get details from ChEBI lookup to check for a direct match to a ChEBI ID
    response = requests.post("https://chebifier.hastingslab.org/api/details", json={
        "type": "type",
        "smiles": smiles_string,
        "selectedModels": {
            "ChEBI Lookup": True
        }
    })
    
    lookup_infotext = response.json().get("models", {}).get("ChEBI Lookup", {}).get("highlights", [])
    
    # If the lookup highlights contain a ChEBI ID, use that for classification instead of the SMILES string. 
    if lookup_infotext and "CHEBI:" in lookup_infotext[0][1]:
        chebi_id = lookup_infotext[0][1].split("CHEBI:")[1].split()[0].rstrip('.')
        chebi_ids.append(f"CHEBI:{chebi_id}")
        print(f"Found ChEBI ID from lookup: CHEBI:{chebi_id} for SMILES {smiles_string}")
    else:
        # Get direct parents from classification
        response = requests.post("https://chebifier.hastingslab.org/api/classify", json={
            "smiles": smiles_string,
            "ontology": False,
            "selectedModels": {
                "ELECTRA (ChEBI50-3STAR)": True,
            }
        })
        
        direct_parents = response.json().get("direct_parents")
        if direct_parents:
            # Extract ChEBI IDs from all parent lists
            for parent_list in direct_parents:
                parent_ids = [f"CHEBI:{parent[0]}" for parent in parent_list]
                chebi_ids.extend(parent_ids)
            print(f"Found {len(chebi_ids)} ChEBI IDs from classification for SMILES {smiles_string}")
    
    return chebi_ids

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
    studyset_list, weights_dict = parse_studyset(raw_studyset)
    
    # Auto-scale weights if present (only scales up if max < 1000)
    # if weights_dict:
    #     weights_dict = auto_scale_weights(weights_dict, target_max=1000)
    
    session['weights_dict'] = weights_dict

    # Get classification from form (allow changing it during re-run)
    classification = request.form.get('classification')
    if classification:
        session['classification'] = classification
    
    looping_prune_method = request.form.get('looping_prune_method')
    if looping_prune_method == 'no_loop_prune':
        # User chose no looping pruning, proceed to other pruning options
        pass
    elif looping_prune_method == 'plain_enrich':
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

        return render_template('results.html', results=results, graph_json_file=graph_json_file)

    # P-VALUE CORRECTION METHOD
    method = request.form.get("p_value_correction_method")
    bonferroni_correct, benjamini_hochberg_correct = map_p_value_correction_method(method)

    # PRUNING OPTIONS
    root_children_prune = request.form.get('root_children_prune') == 'true'
    levels = int(request.form.get('levels', 2))

    linear_branch_prune = request.form.get('linear_branch_prune') == 'true'
    linear_branch_n = int(request.form.get('linear_branch_n', 2))

    high_p_value_prune = request.form.get('high_p_value_prune') == 'true'
    p_value_threshold = float(request.form.get('p_value_threshold', 0.05))

    zero_degree_prune = request.form.get('zero_degree_prune') == 'true'
    
    # Use weighted analysis if weights are present, otherwise use standard analysis
    if weights_dict:
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

    session['correction_method'] = {
        'bonferroni_correct': bonferroni_correct,
        'benjamini_hochberg_correct': benjamini_hochberg_correct
    }

    return render_template('results.html', results=results, graph_json_file=graph_json_file)

@app.route('/graph')
def graph():
    session_id = session.get('session_id')
    graph_file = f'graph_{session_id}.json' if session_id else 'graph.json'
    pruning = session.get('pruning', {})
    correction_method = session.get('correction_method', {})
    classification = session.get('classification', 'structural')
    return render_template('graph.html', 
                         graph_file=graph_file,
                         pruning=pruning, 
                         correction_method=correction_method,
                         classification=classification)



if __name__ == '__main__':
    app.run(debug=True)
