from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import numpy as np
import joblib
import requests
from rdkit import RDLogger
from rdkit.Chem import AllChem, DataStructs
import time


RDLogger.DisableLog('rdApp.*')
clf = joblib.load("tox21_morgan_classifier.pkl")

#smiles obtainer
def smilesname(name):
    try:

        cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
        cid_response = requests.get(cid_url)
        cid_response.raise_for_status()
        cid_data = cid_response.json()

        if 'IdentifierList' not in cid_data or 'CID' not in cid_data['IdentifierList']:
            print("Could not find CID for ",name)
            return None

        cid = cid_data['IdentifierList']['CID'][0]

        smiles_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/ConnectivitySMILES/JSON"
        smiles_response = requests.get(smiles_url)
        smiles_response.raise_for_status()
        smiles_data = smiles_response.json()

        if 'PropertyTable' in smiles_data and 'Properties' in smiles_data['PropertyTable']:
            return smiles_data['PropertyTable']['Properties'][0]['ConnectivitySMILES']
        else:
            print("SMILES not found in response:", smiles_data)
            return None

    except Exception:
        print("Error fetching SMILES:", Exception)
        return None



#1



def morganfin(smiles, radius=2, nBits=1024):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(nBits)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits))

#2



#similarity tester

def tanimoto_similarity(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 is None or mol2 is None:
        return None

    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)

    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    print("Percentage of similarity is:",similarity*100,'%')



#3


while True:

    mol_menu = input(
        "1. Test toxicity levels of a molecule \n2. Obtain SMILES of a molecule \n"
        "3. Generate image of atomic structure of molecule \n4. Compare 2 molecules (Tanimoto Similarity) "
        "\n5. Exit \nEnter your choice (1-5): ")
    time.sleep(1)

    if mol_menu == '1':

        name = input("Enter compound name: ").strip()
        smiles = smilesname(name)

        vec = morganfin(smiles)
        proba = clf.predict_proba([vec])[0]
        mol = Chem.MolFromSmiles(smiles)

        non_toxic_prob = proba[1]
        toxic_prob = proba[0]

        pred = np.argmax(proba)

        print(f"{smiles} → {'✅ Non-toxic' if pred == 1 else '☠️ Toxic'}")
        print(f"🧪 Toxicity Probability: {toxic_prob:.2%}")
        print(f"🧬 Non-Toxic Probability: {non_toxic_prob:.2%}")
        print()
        time.sleep(1)


    elif mol_menu == '2':

        name = input("Enter compound name: ").strip()
        smiles = smilesname(name)

        if smiles:
            print("SMILES for", name, ":", smiles)
        else:
            print("SMILES not found for:", name)
        print()
        time.sleep(1)


    elif mol_menu == '3':

        name = input("Enter compound name: ").strip()
        smiles = smilesname(name)
        mol = Chem.MolFromSmiles(smiles)
        Draw.MolToImage(mol).show()
        print()
        time.sleep(1)



    elif mol_menu == '4':
        name1 = input("Enter first compound's name: ").strip()
        smiles1 = smilesname(name1)
        name2 = input("Enter second compound's name: ").strip()
        smiles2 = smilesname(name2)
        tanimoto_similarity(smiles1, smiles2)
        print()
        time.sleep(1)


    elif mol_menu == '5':
        print("Thanks for trying out ToxiChecker")
        exit()

    else:
        print("Incorrect input, try again ")
        print()
        time.sleep(1)


