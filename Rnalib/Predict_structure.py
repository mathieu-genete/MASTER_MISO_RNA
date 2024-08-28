#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 17:13:59 2023

@author: Mathieu Genete
"""
import numpy as np
from .Rna_seq import Rna_seq
from .Scores import Scores
from .Rna_structure import Rna_structure
import time

try:
    import pandas as pd
    pandas_mod=True
except:
    pandas_mod=False
    
class Predict_structure:
    """
    Classe pour prédire les structures d'ARN.

    Attributs:
        __rna (Rna_seq): Objet représentant la séquence d'ARN.
        __minimal_loop_length (int): Longueur minimale de la boucle.
        __matrix (None): Matrice pour les calculs de structure.
        __structure (None): Structure prédite de l'ARN.
        __all_structures (list): Liste de toutes les structures prédictes.
        __predict_time (int): Temps de prédiction.
        __skipPredAll (bool): Indicateur pour sauter la prédiction de toutes les structures.
        __use_recurse (bool): Indicateur pour utiliser la récursion.
        __bases_scores (Scores): Scores des bases de l'ARN.
    """
    def __init__(self,rnaSeq,minloop=3,skipPredAll=False,use_recurse=False,bases_scores=None): 
        """
        Initialise une instance de Predict_structure.

        Args:
            rnaSeq (Rna_seq): Objet représentant la séquence d'ARN.
            minloop (int, optionnel): Longueur minimale de la boucle. Par défaut à 3.
            skipPredAll (bool, optionnel): Indicateur pour sauter la prédiction de toutes les structures. Par défaut à False.
            use_recurse (bool, optionnel): Indicateur pour utiliser la récursion. Par défaut à False.
            bases_scores (Scores, optionnel): Scores des bases de l'ARN. Par défaut à None.

        Raises:
            Exception: Si rnaSeq n'est pas un objet Rna_seq ou si bases_scores n'est pas un objet Scores.
        """
        if not isinstance(rnaSeq,Rna_seq):
            raise Exception("'{}' n'est pas un objet Rna_seq".format(rnaSeq))
            
        if bases_scores is None:
            self.__bases_scores=Scores()
        elif isinstance(bases_scores,Scores):
            self.__bases_scores=bases_scores
        else:
            raise Exception("'{}' n'est pas un objet Score()".format(bases_scores))
            
        self.__rna=rnaSeq
        self.__minimal_loop_length = int(minloop)
        self.__matrix=None
        self.__structure=None
        self.__all_structures=[]
        self.__predict_time=0
        self.__skipPredAll=skipPredAll
        self.__use_recurse=use_recurse

        #Fait les premières prédictions
        self.structures_prediction(self.__skipPredAll,self.__use_recurse)
            
    #===================
    #Getters Setters
    #===================
    
    @property
    def rna(self):
        """
        Retourne l'objet Rna_seq.

        Returns:
           Rna_seq: Objet représentant la séquence d'ARN.
       """
        return self.__rna
    
    @property
    def matrix(self):
        """
        Retourne la matrice de calcul de structure.

        Returns:
            None: Matrice pour les calculs de structure.
        """
        return self.__matrix
    
    @property
    def structure(self):
        """
        Retourne la structure prédite de l'ARN.

        Returns:
            None: Structure prédite de l'ARN.
        """
        return self.__structure
    
    @property
    def all_structures(self):
        """
        Retourne la liste de toutes les structures prédictes.

        Returns:
            list: Liste de toutes les structures prédictes.
        """
        return self.__all_structures
    
    @property
    def scores(self):
        """
        Retourne les scores des bases de l'ARN.

        Returns:
            Scores: Scores des bases de l'ARN.
        """
        return self.__bases_scores
    
    @property
    def structures_nbr(self):
        """
        Retourne le nombre de structures optimales.

        Returns:
            int: Nombre de structures optimales.
        """
        return len(self.__all_structures)
    
    @property
    def predict_infos(self):
        """
        Retourne les informations de prédiction sous forme de chaîne de caractères formatée.

        Returns:
            str: Informations de prédiction formatées.
        """
        out_lines=[]
        print_rna="["+self.__rna.seq+"]"
        nbr_bases_show=10
        if len(self.__rna.seq)>nbr_bases_show*2:
            print_rna="[{}...{}]".format(self.__rna.seq[:nbr_bases_show],self.__rna.seq[-nbr_bases_show:])
        show_scores=",".join(["{}={}".format(k,v) for k,v in self.__bases_scores.scores.items()])
        out_lines.append("seqID: {} {} ({} bp)".format(self.__rna.id,print_rna,len(self.__rna.seq)))
        out_lines.append("paramètres: [θ={} - {}]".format(self.__minimal_loop_length,show_scores))
        out_lines.append("score max = {}".format(self.__structure.score))
        if len(self.__all_structures)>0:
            out_lines.append("Nombre de structures optimales: {}".format(len(self.__all_structures)))
        out_lines.append("temps de calcul: {}s".format(self.__predict_time))
        maxlength=max([len(v) for v in out_lines])
        stdout="*"*maxlength+"\n"
        stdout+="\n".join(out_lines)+"\n"
        stdout+="*"*maxlength+"\n"
        return stdout
    
    #===================
    #Méthodes publiques
    #===================
    
    def change_scores(self,scores):
        """
	    Change les scores des bases de l'ARN et refait les prédictions de structures.

	    Args:
		scores (Scores): Nouveaux scores des bases de l'ARN.

	    Raises:
		Exception: Si l'objet scores n'est pas une instance de Scores.
	    """
        if isinstance(scores,Scores):
            self.__bases_scores=scores
            self.structures_prediction(self.__skipPredAll,self.__use_recurse)
        else:
            raise Exception("objet Score() attendu")
            
    def print_all_structures(self,filename=None):
        """
	    Affiche ou enregistre toutes les structures prédites.

	    Args:
		filename (str, optionnel): Nom du fichier pour enregistrer les structures. Par défaut à None.
	    """
        if len(self.__all_structures)>0:
            outstruct="\n{}\n".format(self.__rna.seq)
            for s in self.__all_structures:
                outstruct+="{}  - score: {}\n".format(s.dotpar,s.score)
            if filename is not None:
                with open(filename,"w") as outall:
                    outall.write(self.get_predict_info())
                    outall.write(outstruct)
            else:
                print(outstruct)
                
    def export_matrixt_to_csv(self,outcsv,sep=","):
        """
	    Exporte la matrice des scores vers un fichier CSV.

	    Args:
		outcsv (str): Nom du fichier CSV de sortie.
		sep (str, optionnel): Séparateur pour le fichier CSV. Par défaut à ",".
	    """
        size_rna=len(self.__rna.seq)
        with open(outcsv,"w") as outc:
            outc.write(sep.join([" "]+[v for v in self.__rna.seq])+"\n")
            for i in range(size_rna):
                out_line=sep.join([self.__rna.seq[i]]+[str(v) for v in self.__matrix[i]])
                outc.write("{}\n".format(out_line))
            
    def print_matrix(self):
        """
	    Affiche la matrice des scores.

	    Utilise pandas pour un affichage formaté si disponible, sinon utilise numpy.
	    """
        if pandas_mod:
            bases = [*self.__rna.seq]
            print(pd.DataFrame(np.array(self.__matrix), index = bases, columns = bases))
        else:
            print(np.array(self.__matrix))
    
    def structures_prediction(self,skipPredAll=False,use_recurse=False):
        """
	    Prédit les structures d'ARN en remplissant la matrice des scores et en effectuant le traceback.

	    Args:
		skipPredAll (bool, optionnel): Indicateur pour sauter la prédiction de toutes les structures. Par défaut à False.
		use_recurse (bool, optionnel): Indicateur pour utiliser la récursion. Par défaut à False.
	    """
        start_time=time.time()
        #Remplit la matrice des scores
        self.__matrix=self.__fill_mat(self.__rna.seq,self.__minimal_loop_length)
        
        if use_recurse:
            #Traceback en utilisant la récursivité
            fold=[]
            fold = self.__traceback_rec(self.__matrix,self.__rna.seq,self.__minimal_loop_length,fold,0,len(self.__rna.seq)-1)
        else:
            #Traceback en utilisant une pile
            fold = self.__traceback_stack(self.__matrix,self.__rna.seq,self.__minimal_loop_length)
        
        self.__structure = Rna_structure(self.__rna,fold=fold,scores=self.__bases_scores)
        if not skipPredAll:
            self.__all_structures=[]
            fold_list=[]
            self.__traceback_all(self.__matrix,self.__rna.seq,self.__minimal_loop_length,fold_list)
            for fold in fold_list:
                rna_st=Rna_structure(self.__rna,fold=fold,scores=self.__bases_scores)
                if rna_st.check_structure():
                    self.__all_structures.append(rna_st)
        self.__predict_time=time.time()-start_time
        
    #===================
    #Méthodes privées
    #===================
        
    def __init_matrix(self,s):
        """
	    Initialise une matrice de scores pour une séquence d'ARN.

	    Args:
		s (str): Séquence d'ARN.

	    Returns:
		list: Matrice de scores initialisée avec des zéros.
	    """
        return [[0]*len(s) for i in range(0,len(s))]

    def __pairing(self,pair):
        """
	    Retourne le score d'appariement pour une paire de bases.

	    Args:
		pair (tuple): Paire de bases (base1, base2).

	    Returns:
		int: Score d'appariement pour la paire de bases.
	    """
        if pair in self.__bases_scores.pairs.keys():
            return self.__bases_scores.pairs[pair]
        return 0

    def __fill_mat(self,rna,minimal_loop_length):
        """
	    Remplit la matrice de scores pour une séquence d'ARN.

	    Args:
		rna (str): Séquence d'ARN.
		minimal_loop_length (int): Longueur minimale de la boucle.

	    Returns:
		list: Matrice de scores remplie.
	    """
        M=self.__init_matrix(rna)
        for j in range(1,len(rna)):
            for i in range(j):
                if j - i > minimal_loop_length:
                    c1=M[i][j-1]
                    c2=M[i+1][j-1]+self.__pairing((rna[i],rna[j]))
                    c3_list=[M[i][k-1]+self.__pairing((rna[k],rna[j]))+M[k+1][j-1] for k in range(i+1,j-minimal_loop_length)]
                    if len(c3_list)>0:
                        c3=max(c3_list)
                    else:
                        c3=0
                    M[i][j]=max(c1,c2,c3)
                else:
                    M[i][j]=0
        return M
    
    def __traceback_rec(self,M,rna,minimal_loop_length,fold,i,j):
        """
	    Effectue le traceback récursif pour trouver les appariements optimaux.

	    Args:
		M (list): Matrice de scores.
		rna (str): Séquence d'ARN.
		minimal_loop_length (int): Longueur minimale de la boucle.
		fold (list): Liste des appariements trouvés.
		i (int): Indice de début.
		j (int): Indice de fin.

	    Returns:
		list: Liste des appariements optimaux.
	    """
        rna=rna.upper()    
        if j - i > minimal_loop_length:
            if M[i][j]==M[i][j-1]:
                fold=self.__traceback_rec(M, rna,minimal_loop_length, fold, i, j-1)
            elif M[i][j]==M[i+1][j-1]+self.__pairing((rna[i],rna[j])):
                fold.append((i,j))
                self.__traceback_rec(M, rna,minimal_loop_length, fold, i+1, j-1)
            else:
                for k in range(i+1,j-minimal_loop_length):
                    if M[i][j]==M[i][k-1]+self.__pairing((rna[k],rna[j]))+M[k+1][j-1]:
                        fold.append((k,j))
                        self.__traceback_rec(M, rna,minimal_loop_length, fold, i, k-1)
                        self.__traceback_rec(M, rna,minimal_loop_length, fold, k+1, j-1)
                        break
        return fold
    

    def __traceback_str(self,M,rna,minimal_loop_length,i,j,struct=None):
        """
	    Effectue le traceback récursif pour générer la structure en notation dot-bracket.

	    Args:
		M (list): Matrice de scores.
		rna (str): Séquence d'ARN.
		minimal_loop_length (int): Longueur minimale de la boucle.
		i (int): Indice de début.
		j (int): Indice de fin.
		struct (list, optionnel): Structure en cours de construction. Par défaut à None.

	    Returns:
		str: Structure en notation dot-bracket.
	    """
        rna=rna.upper()
        if struct is None:
            struct=['.']*len(rna)
        if j - i > minimal_loop_length:
            if M[i][j]==M[i][j-1]:
                struct=self.__traceback_str(M, rna,minimal_loop_length,i, j-1,struct)
            elif M[i][j]==M[i+1][j-1]+self.__pairing((rna[i],rna[j])):
                struct[i]="("
                struct[j]=")"
                self.__traceback_str(M, rna,minimal_loop_length,i+1, j-1,struct)
            else:
                for k in range(i+1,j-minimal_loop_length):
                    if M[i][j]==M[i][k-1]+self.__pairing((rna[k],rna[j]))+M[k+1][j-1]:
                        struct[k]="("
                        struct[j]=")"
                        self.__traceback_str(M, rna,minimal_loop_length,i, k-1,struct)
                        self.__traceback_str(M, rna,minimal_loop_length,k+1, j-1,struct)
                        break
        return "".join(struct)

    def __traceback_str2(self,M,rna,minimal_loop_length,i,j):
        """
	    Effectue le traceback récursif pour générer la structure en notation dot-bracket (version alternative).

	    Args:
		M (list): Matrice de scores.
		rna (str): Séquence d'ARN.
		minimal_loop_length (int): Longueur minimale de la boucle.
		i (int): Indice de début.
		j (int): Indice de fin.

	    Returns:
		str: Structure en notation dot-bracket.
	    """
        rna=rna.upper()
        if j - i > minimal_loop_length:
            if M[i][j]==M[i][j-1]:
                return self.__traceback_str2(M, rna,minimal_loop_length,i, j-1)+"."
            elif M[i][j]==M[i+1][j-1]+self.__pairing((rna[i],rna[j])):
                return "("+self.__traceback_str2(M, rna,minimal_loop_length,i+1, j-1)+")"
            else:
                for k in range(i+1,j-minimal_loop_length):
                    if M[i][j]==M[i][k-1]+self.__pairing((rna[k],rna[j]))+M[k+1][j-1]:
                        return self.__traceback_str2(M, rna,minimal_loop_length,i, k-1) +"("+ self.__traceback_str2(M, rna,minimal_loop_length,k+1, j-1)+")"
                        break
        elif j==i:
            return "."
        else:
            return self.__traceback_str2(M, rna,minimal_loop_length,i, j-1)+"."
        
    def __traceback_stack(self,M,rna,minimal_loop_length):
        """
	    Effectue le traceback en utilisant une pile pour trouver les appariements optimaux.

	    Args:
		M (list): Matrice de scores.
		rna (str): Séquence d'ARN.
		minimal_loop_length (int): Longueur minimale de la boucle.

	    Returns:
		list: Liste des appariements optimaux.
	    """
        stack=[(0,len(rna)-1)]
        rna=rna.upper()
        fold=[]
        while len(stack)>0:
            i,j=stack.pop()
            if j - i > minimal_loop_length and j>0:
                if M[i][j]==M[i][j-1]:
                    stack.append((i, j-1))
                elif M[i][j]==M[i+1][j-1]+self.__pairing((rna[i],rna[j])) and self.__pairing((rna[i],rna[j]))>0:
                    fold.append((i,j))
                    stack.append((i+1, j-1))
                else:
                    for k in range(i+1,j-minimal_loop_length):
                        if M[i][j]==M[i][k-1]+self.__pairing((rna[k],rna[j]))+M[k+1][j-1] and self.__pairing((rna[k],rna[j]))>0:
                            fold.append((k,j))
                            stack.append((i, k-1))
                            stack.append((k+1, j-1))
                            break
        return fold
    
    
    def __traceback_all(self,M,rna,minimal_loop_length,output,stack=None,fold=None):
        """
	    Effectue le traceback pour trouver toutes les structures optimales.

	    Args:
		M (list): Matrice de scores.
		rna (str): Séquence d'ARN.
		minimal_loop_length (int): Longueur minimale de la boucle.
		output (list): Liste pour stocker toutes les structures optimales.
		stack (list, optionnel): Pile pour le traceback. Par défaut à None.
		fold (list, optionnel): Liste des appariements trouvés. Par défaut à None.
	    """
        if fold is None:
            stack=[(0,len(rna)-1)]
            fold=[]
        while len(stack)>0:
            i,j=stack.pop()
            c=0
            if j - i > minimal_loop_length and j>0:
                tmpstack2=stack.copy()
                tmpfold2=fold.copy()
                if M[i][j]==M[i][j-1] and M[i][j]==M[i+1][j-1]+self.__pairing((rna[i],rna[j])) and self.__pairing((rna[i],rna[j]))>0:
                    c+=1
                    tmpstack=stack.copy()
                    stack.append((i, j-1))
                    tmpfold=fold.copy()
                    tmpfold.append((i,j))
                    tmpstack.append((i+1, j-1))
                    self.__traceback_all(M,rna,minimal_loop_length,output,tmpstack,tmpfold)
                else:
                    if M[i][j]==M[i][j-1]:
                        c+=1
                        stack.append((i, j-1))
                    elif M[i][j]==M[i+1][j-1]+self.__pairing((rna[i],rna[j])) and self.__pairing((rna[i],rna[j]))>0:
                        c+=1
                        fold.append((i,j))
                        stack.append((i+1, j-1))
    
                for k in range(i+1,j-minimal_loop_length):
                    tmpfold3=tmpfold2.copy()
                    tmpstack3=tmpstack2.copy()
                    if M[i][j]==M[i][k-1]+self.__pairing((rna[k],rna[j]))+M[k+1][j-1] and self.__pairing((rna[k],rna[j]))>0:
                        if c==0:
                            c+=1
                            fold.append((k,j))
                            stack.append((i, k-1))
                            stack.append((k+1, j-1))
                        else:
                            c+=1
                            tmpfold3.append((k,j))
                            tmpstack3.append((i, k-1))
                            tmpstack3.append((k+1, j-1))
                            self.__traceback_all(M,rna,minimal_loop_length,output,tmpstack3,tmpfold3)
                            
        output.append(fold)