import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split
import joblib
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

df = pd.read_csv(r"C:\Users\vihit\PycharmProjects\DrugDiscoveryAI\non_toxic_chembl.csv")

target_column = "toxicity"

df = df[['smiles', target_column]].dropna()
df = df[df['smiles'].apply(lambda x: Chem.MolFromSmiles(str(x)) is not None)]

print("Label distribution:")
print(df["toxicity"].value_counts())

def featurize_morgan(smiles, radius=2, nBits=1024):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(nBits)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits))

X = np.array([featurize_morgan(s) for s in df["smiles"]])
y = df["toxicity"].values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.01 , random_state=42, stratify=y)

clf = RandomForestClassifier(n_estimators=300, class_weight="balanced", random_state=42)
clf.fit(X_train, y_train)

y_pred = clf.predict(X_test)
report = classification_report(y_test, y_pred, digits=2, zero_division=0)

joblib.dump(clf, "tox21_morgan_classifier.pkl")

report, len(X_train), len(X_test)
print(report)
print("Training samples:",len(X_train), "Test samples:" ,len(X_test))