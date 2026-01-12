import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from flask import Flask, render_template, request, redirect, url_for, session
from fishers_calculations import run_enrichment_analysis, run_enrichment_analysis_plain_enrich_pruning_strategy
from visualitations_and_pruning import graph_to_cytospace_json
import re

app = Flask(__name__)
app.secret_key = 'your_secret_key'  # Replace with a secure secret key


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

    # Split on commas OR newlines
    items = re.split(r'[,|\n]', studyset)

    # Clean every item: remove quotes + spaces
    cleaned = [x.strip().replace('"', "") for x in items if x.strip()]

    return cleaned

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
    raw_studyset = session.get('study_set')
    if not raw_studyset:
        return redirect(url_for('submission'))
    # Convert multi-line or comma-separated input into a list
    studyset_list = parse_studyset(raw_studyset)

    looping_prune_method = request.form.get('looping_prune_method')
    if looping_prune_method == 'no_loop_prune':
        # User chose no looping pruning, proceed to other pruning options
        pass
    elif looping_prune_method == 'plain_enrich':
        results, pruned_G = run_enrichment_analysis_plain_enrich_pruning_strategy(studyset_list,
                                                classification=session.get('classification'))
        # Save JSON representation of pruned_G in session for graph visualization
        graph_json_file = 'website/static/data/graph.json'
        graph_to_cytospace_json(pruned_G, graph_json_file, results)

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
    
    results, pruned_G = run_enrichment_analysis(studyset_list,
                                      bonferroni_correct=bonferroni_correct,
                                      benjamini_hochberg_correct=benjamini_hochberg_correct,
                                       root_children_prune=root_children_prune,levels=levels,
                                       linear_branch_prune=linear_branch_prune, n=linear_branch_n,
                                       high_p_value_prune=high_p_value_prune, p_value_threshold=p_value_threshold,
                                       zero_degree_prune=zero_degree_prune,
                                       classification=session.get('classification'))
                                       
    # Save JSON representation of pruned_G in session for graph visualization
    graph_json_file = 'website/static/data/graph.json'
    graph_to_cytospace_json(pruned_G, graph_json_file, results)

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
    pruning = session.get('pruning', {})
    correction_method = session.get('correction_method', {})
    return render_template('graph.html', pruning=pruning, correction_method=correction_method)


if __name__ == '__main__':
    app.run(debug=True)
