#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 17:09:47 2023

@author: Mathieu Genete
"""
from .Alphabet import Alphabet
from .Rna_structure import Rna_structure
from .Rna_seq import Rna_seq
from .Predict_structure import Predict_structure

import os

class Rna_parser:
    @staticmethod
    def parse_fasta(filename: str):
        """
        Analyse un fichier FASTA et retourne un dictionnaire avec les IDs de séquence comme clés et leurs descriptions et séquences comme valeurs.

        Args:
            filename (str): Le chemin vers le fichier FASTA.

        Returns:
            dict: Un dictionnaire où les clés sont les IDs de séquence et les valeurs sont des dictionnaires avec les clés 'description' et 'seq'.
        """
        outseq={}
        with open(filename) as fasta:
            id=""
            while True:
                line=fasta.readline().strip()
                if line=="":
                    break
                if line[0]==">":
                    id=line.split(" ")[0][1:]
                    outseq[id]={'description':line[1:],'seq':""}
                elif id!="":
                    outseq[id]['seq']+=line.strip().upper()
        return outseq
    
    @staticmethod
    def fasta_to_db(infasta: str,outdb: str,minloop=3,scores=None):
        """
        Convertit un fichier FASTA au format dot-bracket en prédisant les structures d'ARN.

        Args:
            infasta (str): Le chemin vers le fichier FASTA d'entrée.
            outdb (str): Le chemin vers le fichier dot-bracket de sortie.
            minloop (int, optionnel): La longueur minimale de la boucle pour la prédiction de structure. Par défaut à 3.
            scores (dict, optionnel): Un dictionnaire de scores de bases pour la prédiction de structure. Par défaut à None.

        Returns:
            None
        """
        fasta=Rna_parser.parse_fasta(infasta)
        struct_list=[]
        for id,datas in fasta.items():
            rna_seq=Rna_seq(datas['description'],datas['seq'])
            struct_list.append(Predict_structure(rna_seq,minloop,skipPredAll=True,bases_scores=scores))
            
        with open(outdb,"w") as outdb_file:
            for s in struct_list:
                outdb_file.write(">{}\n{}\n{}\n".format(s.rna.id,s.rna.seq,s.structure.dotpar))
    
    @staticmethod
    def parse_dotbrackets_file(filename: str,minimal_loop_length=3):
        """
        Analyse un fichier de dot-bracket et retourne un dictionnaire de structures d'ARN.

        Args:
            filename (str): Le chemin vers le fichier de dot-bracket.
            minimal_loop_length (int, optionnel): La longueur minimale de la boucle pour la validation des épingles à cheveux. Par défaut à 3.

        Returns:
            dict: Un dictionnaire où les clés sont les IDs de séquence et les valeurs sont des objets de structure d'ARN.

        Raises:
            Exception: Si le fichier n'existe pas ou s'il y a des erreurs de format dans les séquences.
        """
        dotpar_alphabet=Alphabet.dotpar()
        rna_alphabet=Alphabet.rna()
        rnastruct_list={}
        out_error=[]
        if not os.path.exists(filename):
            raise Exception(f"Le fichier '{filename}' n'existe pas")
            
        with open(filename,'r') as dbfile:
            while True:
                line_id=dbfile.readline().strip()
                seq=dbfile.readline().strip().upper()
                dotpar=dbfile.readline().strip()
                
                if line_id=="":
                    break

                #check rna sequence alphabet
                check_rna_alphabet=all([b in rna_alphabet for b in seq])
                
                #check dotpar sequence
                check_dotpar=all([b in dotpar_alphabet for b in dotpar])
                
                if line_id.startswith(">") and check_rna_alphabet and check_dotpar:
                    id=line_id[1:]
                    rna=Rna_seq(id,seq)
                    struct=Rna_structure(rna,dotpar=dotpar)
                    if struct.check_structure() and struct.check_hairpin(minimal_loop_length):
                        rnastruct_list[id]=struct
                    else:
                        out_error.append("\t=> Erreurs dans la séquence {}".format(id))
                else:
                    out_error.append("\t=> Erreurs dans la séquence {}".format(line_id))
                    
        if len(out_error)>0:
            exception_txt="Erreur(s) de format parenthésé:\n{}".format("\n".join(out_error))
            raise Exception(exception_txt)
            
        return rnastruct_list
    
    @staticmethod
    def parse_connect_file(filename:str ,minimal_loop_length=3):
        """
        Analyse un fichier au format connect et retourne un objet de structure d'ARN.

        Args:
            filename (str): Le chemin vers le fichier au format connect.
            minimal_loop_length (int, optionnel): La longueur minimale de la boucle pour la validation de structure. Par défaut à 3.

        Returns:
            Rna_structure: Un objet de structure d'ARN.

        Raises:
            Exception: S'il y a des erreurs de format dans le fichier connect.
        """
        seq=""
        fold=[]
        fold_set=set()
        with open(filename,"r") as inct:
            for line in inct:
                tmp=line.strip().split()
                if len(tmp)==3 and int(tmp[2])>0:
                    pos=int(tmp[0])-1
                    pos_link=int(tmp[2])-1
                    fold.append((pos,pos_link))
                    fold_set.add((min(pos,pos_link),max(pos,pos_link)))
                seq+=tmp[1].upper()
        if not Rna_parser.__check_connect_format(seq,fold,fold_set,minimal_loop_length):
            raise Exception("Erreur dans la verification du format connect")
        rna=Rna_seq(filename,seq)
        return Rna_structure(rna,fold=fold_set)
        
    #================
    #Méthodes privées
    #================ 
    
    def __check_connect_format(seq: str,fold: list,fold_set: set,minimal_loop_length: int):
        """
        Vérifie le format d'un fichier connect.
    
        Args:
            seq (str): La séquence d'ARN.
            fold (list): Liste des paires de positions de liaison.
            fold_set (set): Ensemble des paires de positions de liaison.
            minimal_loop_length (int): Longueur minimale de la boucle.
    
        Returns:
            bool: Retourne True si le format est correct, sinon lève une exception.
    
        Raises:
            Exception: Si des erreurs de format sont détectées.
        """
        out_error=[]
        for b in seq:
            if b not in Alphabet.rna():
                out_error.append("La séquence d'ARN n'est pas au bon format")
                
        open_pos=sorted([v[0] for v in fold])
        close_pos=sorted([v[1] for v in fold])
        
        crossings=Rna_parser.__check_connect_croisements(fold_set)
        
        if len(crossings)>0:
            out_error.append("\t=>il y a des croisements entre les positions")
            out_error=out_error+crossings
        
        if not open_pos==close_pos:
            out_error.append("\t=>les indexes ne correspondent pas entre les positions")
        
        if len([1 for v in fold_set if abs(v[1]-v[0])<minimal_loop_length])>0:
            out_error.append("\t=>une boucle est inférieure à {}".format(minimal_loop_length))
            
        if len(out_error)>0:
            exception_txt="Erreur(s) de format connect:\n{}".format("\n".join(out_error))
            raise Exception(exception_txt)
        
        return True
    
    def __check_connect_croisements(fold_set: set):
        """
        Vérifie les croisements dans les paires de positions de liaison.
    
        Args:
            fold_set (set): Ensemble des paires de positions de liaison.
    
        Returns:
            list: Liste des erreurs de croisements détectées.
        """
        fold=sorted(fold_set)
        crossings=[]
        for n in range(len(fold)-1):
            i,j=fold[n]
            k,l=fold[n+1]
            if not i<k<l<j and not i<j<k<l:
                crossings.append(f"\t\tCroisement - i={i} - j={j} - k={k} - l={l}")
        return crossings