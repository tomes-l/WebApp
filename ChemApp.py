from flask import Flask, render_template, request, jsonify
from static.python_class.Class_molecule import Molecule 
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import Draw
import base64
from io import BytesIO

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

### Constants ###
ORGANIC = ["H", "B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"]
###############

# Function define for SMILES input
def SMILES_algo(dict_commandes_values):
    dict_return_value = {}

    if not dict_commandes_values['SMILES']:
        dict_return_value['error'] = 'No smiles provided'
        return dict_return_value
    else:
        smiles_string = dict_commandes_values['SMILES']
        mol = Chem.MolFromSmiles(smiles_string)

        if dict_commandes_values.get('mol'):
            dict_return_value['mol'] = Chem.MolToMolBlock(mol)

        if dict_commandes_values.get('img'):
            representation = Chem.MolFromSmiles(smiles_string)
            img = Draw.MolToImage(representation)
    
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue()).decode('utf-8')

            dict_return_value['img'] = img_str
        
        if dict_commandes_values.get('mformula'):
            dict_return_value['mformula'] = CalcMolFormula(mol)

        if dict_commandes_values.get('mweight'):
            dict_return_value['mweight'] = Descriptors.ExactMolWt(mol)

    return dict_return_value

# Function define for mol input
def mol_algo(dict_commandes_values):
    dict_return_value = {}
    molecule = Molecule()

    if dict_commandes_values['file'] == None:
        dict_return_value['error'] = 'No file provided'
        return dict_return_value
    else:
        mol_str = dict_commandes_values['file']
        molecule.load_file(mol_str)
        molecule.initialize_variables()
        molecule.build_graph()
        verify = molecule.verification(ORGANIC)
        if verify == 3:
            dict_return_value['error'] = 'File could not be loaded, bonds block does not cover all atoms'
            return dict_return_value
        elif verify == 2:
            dict_return_value['error'] = 'File could not be loaded, bond blocks are incomplete'
            return dict_return_value
        elif verify == 1:
            dict_return_value['error'] = 'File could not be loaded, contains non-organic elements'
            return dict_return_value
    
    if dict_commandes_values.get('SMILES'):
        SMILES = molecule.to_smiles()
        dict_return_value['smiles'] = SMILES

    if dict_commandes_values.get('img'):
        representation = Chem.MolFromSmiles(SMILES)
        img = Draw.MolToImage(representation)
    
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode('utf-8')

        dict_return_value['img'] = img_str
        
    if dict_commandes_values.get('count_atom') != None:
        count = molecule.count(dict_commandes_values['count_atom'])
        if count > 0:
            dict_return_value['count_atom'] = f"The number of {dict_commandes_values['count_atom']} atom elements in the molecule is {str(count)}."
        else:
            dict_return_value['count_atom'] = f"The atom element {dict_commandes_values['count_atom']} does not exist in the molecule."
    
    if dict_commandes_values.get('room_distance1') and dict_commandes_values.get('room_distance2') != None:
        atom_one = int(dict_commandes_values['room_distance1'])
        atom_two = int(dict_commandes_values['room_distance2'])
        distance = molecule.room_distance(atom_one, atom_two)
        if distance == "Wrong ID":
            dict_return_value['room_distance'] = "Invalid atom identifier provided"
        elif distance == "No XYZ":
            dict_return_value['room_distance'] = "Cannot compute the 3D distance, coordinates are not available"
        elif molecule.coordinates_sum() == 0:
            dict_return_value['room_distance'] = "Cannot compute the 3D distance, coordinates are not available"
        else:
            dict_return_value['room_distance'] = f"The 3D distance between atoms {str(atom_one)} and {str(atom_two)} is {str(distance)}"
    
    if dict_commandes_values.get('bond_distance1') and dict_commandes_values.get('bond_distance2') != None:
        atom_one = int(dict_commandes_values['bond_distance1'])
        atom_two = int(dict_commandes_values['bond_distance2'])
        if molecule.path_distance(atom_one, atom_two) == "False":
            dict_return_value['bond_distance'] = "Invalid atom identifier provided"
        elif molecule.path_distance(atom_one, atom_two) == 0:
            dict_return_value['bond_distance'] = f"The 2D distance between atoms {str(atom_one)} and {str(atom_two)} is {str(0)}"
        else:
            count = molecule.path_distance(atom_one, atom_two)
            dict_return_value['bond_distance'] = f"The 2D distance between atoms {str(atom_one)} and {str(atom_two)} is {str(count)}"
    
    if dict_commandes_values.get('atom_neighbours') != None:
        atom_id = int(dict_commandes_values['atom_neighbours'])
        if molecule.get_neighbour_num(atom_id) == False:
            dict_return_value['atom_neighbours'] = "Invalid atom identifier provided"
        else:
            neighbour_num = molecule.get_neighbour_num(atom_id)
            dict_return_value['atom_neighbours'] = f"Atom {str(atom_id)} has {str(neighbour_num)} neighbours"
    
    if dict_commandes_values.get('ring'):
        results = molecule.find_ring(ORGANIC)
        if results == 0:
            dict_return_value['ring'] = "There is no ring in the molecule"
        elif results == 1:
            dict_return_value['ring'] = f"There is {results} ring in the molecule"
        else:
            dict_return_value['ring'] = f"There are {results} rings in the molecule"
    
    return dict_return_value


@app.route('/submit-form-smiles', methods=['POST'])
def handle_form_smiles():
    data = request.get_json()
    dict_SMILES = {
        'SMILES': data.get('SMILES-input'),
        'mol': data.get('SMILES-to-MolFile'),
        'img': data.get('SMILES-to-representation'),
        'mformula': data.get('Molecular-formula'),
        'mweight': data.get('Molecular-weight')
    }
    print("Dict SMILES:", dict_SMILES)
    response = SMILES_algo(dict_SMILES)
    print("Response:", response)
    return jsonify(response)


@app.route('/submit-form-mol', methods=['POST'])
def handle_form_mol():
    data = request.get_json()
    dict_CM = {
        'file': data.get('file'),
        'SMILES': data.get('Molfile-to-Smiles'),
        'img': data.get('Molfile-to-representation'),
        'count_atom': data.get('count-atom'),
        'bond_distance1': data.get('2D-distance1'),
        'bond_distance2': data.get('2D-distance2'),
        'room_distance1': data.get('3D-distance1'),
        'room_distance2': data.get('3D-distance2'),
        'atom_neighbours': data.get('atom-neighbours'),
        'ring': data.get('Is-ring')
    }
    print("Dict CM: ", dict_CM)
    response = mol_algo(dict_CM)
    print("Reponse: ", response)
    return jsonify(response)

if __name__ == '__main__':
    app.run(debug=True)
