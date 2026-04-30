# import convert smiles to chebi from app.py
from flask import session
from website.app import app, convert_smiles_to_chebi

print("Testing convert_smiles_to_chebi function...")
smiles = "CCCCCCCCCCCCCC=CC(C(CO)N)O"

with app.test_request_context("/"):
	# Match your app logic: choose whether to classify with parent terms.
	session["smiles_option"] = "use_parents"
	converted = convert_smiles_to_chebi(smiles)

print(converted)