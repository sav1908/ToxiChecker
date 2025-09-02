import requests
import pandas as pd
import time

def fetch_approved_drugs(max_molecules=2000):
    url = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
    params = {
        "molecule_properties__aromatic_rings__gt": 0,
        "max_phase": 4,
        "molecule_structures__isnull": False,
        "limit": 1000,
        "offset": 0
    }

    all_smiles = []

    while len(all_smiles) < max_molecules:
        response = requests.get(url, params=params)
        if response.status_code != 200:
            break

        data = response.json()
        molecules = data.get("molecules", [])

        for mol in molecules:
            smi = mol.get("molecule_structures", {}).get("canonical_smiles")
            if smi:
                all_smiles.append(smi)

            if len(all_smiles) >= max_molecules:
                break

        if not data.get("page_meta", {}).get("next"):
            break

        params["offset"] += params["limit"]
        time.sleep(0.5)

    return all_smiles


smiles_list = fetch_approved_drugs(2000)
df = pd.DataFrame({
    "smiles": smiles_list,
    "toxicity": [1] * len(smiles_list)
})

df.to_csv("non_toxic_chembl.csv", index=False)
print(f"✅ Saved {len(df)} non-toxic drugs to 'non_toxic_chembl.csv'")