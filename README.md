# rnalib

`rnalib` est une bibliothèque Python dédiée à la manipulation et à l'analyse des séquences d'ARN. Elle permet de prédire et de comparer les structures secondaires de l'ARN en utilisant l'algorithme de Nussinov, ainsi que d'exporter les résultats dans différents formats.

## Utilisation

### Création d'une séquence ARN

Pour créer une séquence ARN, utilisez la classe `Rna_seq` :

```python
import Rnalib

rna_seq = Rnalib.Rna_seq("nom_de_la_sequence", "sequence_ARN")
```

### Prédiction de la structure

Pour prédire la structure d'une séquence ARN, utilisez la classe `Predict_structure` :

```python
a = Rnalib.Predict_structure(rna_seq, 3, skipPredAll=True)
```
la variable `skipPredAll` permet déterminer ou non de calculer l'ensemble des structures optimales

### Affichage des informations de prédiction

Pour afficher les informations de prédiction :

```python
print(a.predict_infos)
```

### Exportation de la matrice des scores

Pour exporter la matrice des scores dans un fichier CSV :

```python
a.export_matrixt_to_csv("nom_du_fichier.csv")
```

### Affichage de la structure prédite

Pour afficher la première structure prédite :

```python
print("score:", a.structure.score)
print("fold:", a.structure.fold)
print("format parenthèsé:", a.structure)
```

### Vérification de la structure

Pour vérifier la structure prédite :

```python
print("Vérification de la structure:", a.structure.check_structure())
```

### Comparaison des structures

Pour comparer deux structures :

```python
st = Rnalib.Rna_parser.parse_dotbrackets_file("nom_du_fichier.db", 3)
compacted_tree_s2 = st['sequence_test2_s2'].arbre.compact_tree()
compacted_tree_s3 = st['sequence_test2_s3'].arbre.compact_tree()

if compacted_tree_s2 == compacted_tree_s3:
    print("Les deux structures sont identiques")
else:
    print("Les deux structures ne sont pas identiques")
```

### Prédiction des structures à partir d'un fichier FASTA

Pour prédire les structures à partir d'un fichier FASTA et exporter les résultats dans un fichier dot-brackets :

```python
Rnalib.Rna_parser.fasta_to_db("sequences_human_tRNA.fasta", "nom_du_fichier.db")
```

### Ouverture et affichage des structures à partir d'un fichier CT

Pour ouvrir et afficher les structures à partir d'un fichier CT :

```python
ctstruct = Rnalib.Rna_parser.parse_connect_file("nom_du_fichier.ct")
ctstruct.print_struct()
```
