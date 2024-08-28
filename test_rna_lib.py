#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 23:16:28 2023

@author: mathieu
"""
import Rnalib

#rna_seq=Rnalib.Rna_seq("test1","CGGAUACUUCUUA")
#rna_seq=Rnalib.Rna_seq("test2","CGGAUACUUCUUAGACGA")
#rna_seq=Rnalib.Rna_seq("test2","CCCCUUAAUUUUUUGGGG")
#rna_seq=Rnalib.Rna_seq("test3","GGGGAAACCCAGGUUCGUUUCGGUCAAGACAACCC")
#rna_seq=Rnalib.Rna_seq("test4","GGGGGCGUAGCUCAGAuGGUAGAGCGCUCGCUUgGCgUGUGAGAGGUACCGGGAUCGaUACCCGGCGCCUCCACCA")
#rna_seq=Rnalib.Rna_seq("sequence_test3","CAGAGUAUGAUCACGGUUUCACCUUGGUACAGGGCGUUCCACUGCACUCUG")
rna_seq=Rnalib.Rna_seq("test4","GAAAGAGAAATGTGAAATGAAGAGAACAAGAATTGTTTGAATAGCCAAGGAGACTGCCTGACGACTAACCAACTCTAAGACTATCGATCGATCGAGTAGTCTTTGATATTGTATGGTTTAAAACGCAGGCAAGTCATCCTTGGCTACCACACAAGTTCTCATTCTTCATGTCACCATTTTCTCTTA")
#rna16s_Ecoli="AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTA"
#,skipPredAll=True
#rna_seq=Rnalib.Rna_seq("rna16s_Ecoli",rna16s_Ecoli)

#rna_seq=Rnalib.Rna_seq("tRNA-Ala_hs","AAGGGCTTAGCTTAATTAAAGTGGCTGATTTGCGTTCAGTTGATGCAGAGTGGGGTTTTGCAGTCCTTACCA")
#rna_seq=Rnalib.Rna_seq("test_1","UCUUCUACGCCGGGGUAGCCUCUGCUAACCUACGCUCAUCGCAACCAGUAGCUAGUGAAGCAAAUGAUUAUGUGUGUUUUUGGCACGACUCCUAUCUCCUGGACUGCUGUAAGAUCCACCGUUUUUGGCAUGCCCGGUCGCGACCUCAUCAAAUUGGUGCAUGAGAAUUAUGUCUCGGUUCUUUUUGCAAAGAAAGUCAUUCAAGUCCACUCCAACGGGACGAAGAAGGAUUUCUGCAUCAUUCGCAUGCACCAAUUUCAUGAGGUCGCGACCGGGCAUGCCAAAAACGGUGGAUCUGACAGCAGCCCAGGAGAUAGGAGUCGUGCCAAAAACGCACAGAAUCAUUAGCUUCACUCUUGCUGGUUGCGAUGAUGAGUGAACUCUACCCCCGCGUAGAAGA")

#rna_seq=Rnalib.Rna_seq("test1","GCAAAAAGCUUAAGGGAAAACCUCCAUUCCCC")
#rna_seq=Rnalib.Rna_seq("exemple_4","CAUCCCCAGCGGGUCCCGCGAGCUCACCCCCUAUACCAGA")
#rna_seq=Rnalib.Rna_seq("test1","GGUUUCUCUUUGUC")
#test_seq="gcgaggcuagcgcuacccgugcgccugcguggaacgauucuguggcgagugccggccgaaagcuagguccggauugcacguggagggccgcccgaagggcacucucggacauuaacccgcauucuguaccauggggcgcaaguuggacccuacgaaggagaagcgggggccaggccgaaaggcccggaagcagaagggugccgagacagaacucgucagauucuugccugcaguaagugacgaaaauuccaagaggcugucuagucgugcucgaaagagggcagccaagaggagauugggcucuguugaagccccuaagacaaauaagucuccugaggccaaaccauugccuggaaagcuaccaaaaggagcuguccagacagcugguaagaagggaccccagucccuauuuaaugcuccucgaggcaagaagcgcccagcaccuggcagugaugaggaagaggaggaggaagacucugaagaagaugguauggugaaccacggggaccucuggggcuccgaggacgaugcugauacgguagaugacuauggagcugacuccaacucugaggaugaggaggaaggugaagcguugcugcccauugaaagagcugcucggaagcagaaggcccgggaagcugcugcugggauccaguggagugaagaggagaccgaggacgaggaggaagagaaagaagugaccccugagucaggccccccaaagguggaagaggcagaugggggccugcagaucaauguggaugaggaaccauuugugcugcccccugcuggggagauggagcaggaugcccaggcuccagaccugcaacgaguucacaagcggauccaggauauugugggaauucugcgugauuuuggggcucagcgggaggaagggcggucucguucugaauaccugaaccggcucaagaaggaucuggccauuuacuacuccuauggagacuuccugcuuggcaagcucauggaccucuuc"
#rna_seq=Rnalib.Rna_seq("test1",test_seq)

#a=Rnalib.Predict_structure(rna_seq,3,bases_scores=Rnalib.Scores(GC=1,AU=1,GU=1))
a=Rnalib.Predict_structure(rna_seq,3,skipPredAll=True)


print(a.predict_infos)

print("\n\nMatrice des scores:")
a.print_matrix()
a.export_matrixt_to_csv("test.csv")

print("\nAffichage de la première structure prédite:")
print("score:",a.structure.score)
print("fold:",a.structure.fold)
print("\nformat parenthèsé:")
print(a.rna.seq)
print(a.structure)
print()

a.structure.print_struct(sepsize=2,print_pos=True)

print("\nformat CT: fichier 'test.ct'")
a.structure.structure_to_ct("test.ct")


print("\nVérification de la structure: ",a.structure.check_structure(),"\n")

print("Arbre de la structure:\n")
a.structure.arbre.print_tree()
print("structure issue de l'arbre:",a.structure.arbre.tree_to_dotpar())
print("structure compactée (utilisant dotpar):\t",a.structure.compact_struct())
print("structure compactée (utilisant arbre):\t",a.structure.arbre.compact_tree().tree_to_dotpar(),"\n")
print("Arbre de la structure compactée:\n")
a.structure.arbre.compact_tree().print_tree(print_tuples=False)

print("\n==========================================\n")

print("Liste des structures optimales:")
print("\nil y a {} structure(s) optimale(s)\n".format(a.structures_nbr))
a.print_all_structures()
    
print("\n==========================================\n")

st=Rnalib.Rna_parser.parse_dotbrackets_file("RNA_seq_struct_v2.db",3)

print("Comparaison des structures de sequence_test2_s2:",st['sequence_test2_s2'].dotpar)
print("\net de sequence_test2_s3:",st['sequence_test2_s3'].dotpar,"\n")
compacted_tree_s2=st['sequence_test2_s2'].arbre.compact_tree()
compacted_tree_s3=st['sequence_test2_s3'].arbre.compact_tree()

print("structure compactée de sequence_test2_s2",compacted_tree_s2.tree_to_dotpar())
print("structure compactée de sequence_test2_s3",compacted_tree_s3.tree_to_dotpar())

if compacted_tree_s2==compacted_tree_s3:
    print("\nLes deux structures sont identiques")
else:
    print("\nLes deux structures ne sont pas identiques")
    
print("\n==========================================\n")


fasta_file="sequences_human_tRNA.fasta"
out_db_file="test_fasta.db"
print("Prédiction des structures pour le fichier fasta:",fasta_file,"\n")
print(f"les prédictions sont écrites dans le fichier '{out_db_file}'")
Rnalib.Rna_parser.fasta_to_db(fasta_file,out_db_file)


print("\n==========================================\n")
print("open RNA_seq_struct.db:")
dbfile=Rnalib.Rna_parser.parse_dotbrackets_file("RNA_seq_struct.db")

for id,st in dbfile.items():
    print("id:",id,"- score:",st.score," - structure:",st.dotpar," - compact:",st.compact_struct())

print("\n==========================================\n")
print("open RA7680.ct:")
ctstruct = Rnalib.Rna_parser.parse_connect_file("RA7680.ct")
ctstruct.print_struct()