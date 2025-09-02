

#1
'''
name = input("Enter compound name: ").strip()
smiles = smilesname(name)

if smiles:
    print("SMILES for",name,":", smiles)
else:
    print("SMILES not found for:", name)

'''

#2
'''
vec = morganfin(smiles)
proba = clf.predict_proba([vec])[0]
mol = Chem.MolFromSmiles(smiles)


non_toxic_prob = proba[1]
toxic_prob = proba[0]

pred = np.argmax(proba)

print(f"{smiles} → {'✅ Non-toxic' if pred == 1 else '☠️ Toxic'}")
print(f"🧪 Toxicity Probability: {toxic_prob:.2%}")
print(f"🧬 Non-Toxic Probability: {non_toxic_prob:.2%}")

image = input("Would you like to see the molecule? ")
if image.lower() in "yes":
    Draw.MolToImage(mol).show()
else:
    exit()

'''

#3

'''sim = input("Would you like to compare 2 compounds? ")
if sim in "yes":
    name1 = input("Enter first compound's name: ").strip()
    smiles1 = smilesname(name1)
    name2 = input("Enter second compound's name: ").strip()
    smiles2 = smilesname(name2)
    print(tanimoto_similarity(smiles1,smiles2))

'''

