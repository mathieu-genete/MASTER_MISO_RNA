#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 17:11:36 2023

@author: Mathieu Genete
"""
import re
from .Tree import Tree
from .Rna_seq import Rna_seq
from .Scores import Scores

class Rna_structure:
    """
    Classe représentant la structure d'un ARN.
    
    Attributs:
    ----------
    __rna : Rna_seq
        Séquence ARN.
    __fold : list
        Liste des paires de bases formant la structure.
    __dotpar : str
        Représentation en notation dot-parenthèse de la structure.
    __scores : Scores
        Scores associés aux paires de bases.
    __arbre : Tree
        Arbre représentant la structure.
    __score : int
        Score total de la structure.
    """
    
    def __init__(self,rnaSeq: Rna_seq,fold=None,scores=None,dotpar=None):
        """
        Initialise une nouvelle instance de Rna_structure.
        
        Paramètres:
        -----------
        rnaSeq : Rna_seq
            Objet représentant la séquence ARN.
        fold : list, optionnel
            Liste des paires de bases formant la structure.
        scores : Scores, optionnel
            Scores associés aux paires de bases.
        dotpar : str, optionnel
            Représentation en notation dot-parenthèse de la structure.
        
        Exceptions:
        -----------
        Exception
            Si rnaSeq n'est pas un objet Rna_seq.
        Exception
            Si scores n'est pas un objet Scores.
        Exception
            Si ni fold ni dotpar, ou si les deux sont fournis.
        """
        if not isinstance(rnaSeq,Rna_seq):
            raise Exception("'{}' n'est pas un objet Rna_seq".format(rnaSeq))
        
        if scores is not None and not isinstance(scores,Scores):
            raise Exception("'{}' n'est pas un objet Score()".format(scores))
            
        self.__rna=rnaSeq
        if (fold is None and dotpar is None) or (fold is not None and dotpar is not None):
            raise Exception("une structure RNA_structure requière soit un fold ou soit un dotpar")
        elif fold is None:
            self.__dotpar=dotpar
            self.__fold=sorted(self.__dot_par_to_bp(dotpar))
            
        elif dotpar is None:
            self.__fold=sorted(fold)
            self.__dotpar=self.__fold_to_dotpar()
        if scores is None:
            self.__scores = Scores()
        else:
            self.__scores=scores
        
        if self.check_structure():
            #self.__arbre_s=self.__arbre_s(self.__dotpar)
            fold_dict={k:v for k,v in sorted(self.__fold)}
            self.__arbre=self.__construct_tree(Tree((-1,-1)),fold_dict,0,len(self.__rna.seq))
        
        rna=self.__rna.seq
        self.__score=0
        if self.__dotpar and self.__dot_par_to_bp(self.__dotpar):
            for i,j in self.__dot_par_to_bp(self.__dotpar):
                pair=(rna[i],rna[j])
                if pair in self.__scores.pairs:
                    self.__score+=self.__scores.pairs[pair]

    #===================
    #Getters Setters
    #===================
    @property
    def rna(self):
        """Retourne la séquence ARN."""
        return self.__rna
    
    @property
    def dotpar(self):
        """Retourne la représentation en notation dot-parenthèse de la structure."""
        return self.__dotpar
    
    @property
    def fold(self):
        """Retourne la liste des paires de bases formant la structure."""
        return self.__fold
        
    @property
    def score(self):
        """Retourne le score total de la structure."""
        return self.__score
    
    @property
    def arbre(self):
        """Retourne l'arbre représentant la structure."""
        return self.__arbre

    #===================
    #Méthodes magiques
    #=================== 
    def __str__(self):
        """Retourne la représentation en notation dot-parenthèse de la structure."""
        return self.__dotpar
    
    #===================
    #Méthodes publiques
    #===================    
    def structure_to_ct(self,filename=None):
        """
        Convertit la structure au format connect (CT).
        
        Paramètres:
        -----------
        filename : str, optionnel
            Nom du fichier où sauvegarder la table CT. Si None, retourne la table CT sous forme de chaîne de caractères.
        
        Retourne:
        ---------
        str
            La table CT sous forme de chaîne de caractères si filename est None.
        """
        dotpar=self.__dotpar
        rna_seq=self.__rna.seq
        outct=self.__dot_par_to_bp(dotpar)
        third_col=[0]*len(dotpar)
        for i,j in outct:
            third_col[i]=j+1
            third_col[j]=i+1
        txtCT=""
        for i in range(0,len(rna_seq)):
            txtCT+="{}\t{}\t{}\n".format(i+1,rna_seq[i],third_col[i])
            
        if filename is not None:
            with open(filename,"w") as foutCT:
                foutCT.write(txtCT)
        else:
            return txtCT
    
    def check_structure(self):
        """
        Vérifie la validité de la structure ARN.
        
        Retourne:
        ---------
        bool
            True si la structure est valide, False sinon.
        """
        #verifie que le nombre de parenthèses ouvrantes et fermantes sont identiques
        if not self.__parens_count(self.__dotpar):
            return False
        
        if not self.__base_pairs_check():
            return False
        
        return True
    
    def check_hairpin(self,minimal_loop_length: int):
        """
        Vérifie que toutes les boucles de la structure respectent une longueur minimale.
        
        Paramètres:
        -----------
        minimal_loop_length : int
            Longueur minimale des boucles.
        
        Retourne:
        ---------
        bool
            True si toutes les boucles respectent la longueur minimale, False sinon.
        """
        fold=self.__fold
        for bp in fold:
            if bp[1] - bp[0] <= minimal_loop_length:
                print("une boucle est inférieure à {}".format(minimal_loop_length))
                return False
        return True
    
    def compact_struct(self):
        """
        Compacte la structure en supprimant les paires de bases adjacentes.
        
        Retourne:
        ---------
        str
            La structure compactée en notation dot-parenthèse.
        """
        struct=list(self.__dotpar)
        fold_dict={i:j for i,j in self.__fold}
        for i in range(0,len(struct)-1):
            if struct[i]=="(" and struct[i+1]=="(":
                if fold_dict[i]-fold_dict[i+1]==1:
                    struct[i]=""
                    struct[fold_dict[i]]=""
        output="".join(struct)
        
        output= re.sub(r"\.+",".",output)
        return output
    
    def dot_par_to_latex(self,dotb=None,print_struct=True,numbers_shift=1):
        """
        Convertit la structure en notation dot-parenthèse en code LaTeX.
        
        Paramètres:
        -----------
        dotb : str, optionnel
            Structure en notation dot-parenthèse. Si None, utilise la structure de l'objet.
        print_struct : bool, optionnel
            Si True, imprime la structure. Sinon, imprime les positions des paires de bases.
        numbers_shift : int, optionnel
            Décalage des numéros de positions.
        """
        if dotb is None:
            dotb=self.__dotpar
        outstr=""
        for i in range(0,len(dotb)):
            rslt=self.__search_link(dotb,i)
            rslt=(rslt[0]+numbers_shift,rslt[1]+numbers_shift)
            if dotb[i]=="(":
                if print_struct:
                    outstr+="[.,green "
                else:
                    outstr+="[{$"+str(rslt)+"$},green "
            elif dotb[i]==".":
                if print_struct:
                    outstr+="[.,blue]"
                else:
                    outstr+="[{$"+str(rslt)+"$},blue]"
            elif dotb[i]==")":
                outstr+="]"
        print("\\begin{forest}")
        if print_struct:
            print("for tree={circle,fill,l sep=30pt}")
        else:
            print("for tree={l sep=30pt}")
        print("["+outstr+"]")
        print("\\end{forest}")
        
    def print_struct(self,sepsize=0,print_pos=False):
        """
        Affiche la structure ARN.
        
        Paramètres:
        -----------
        sepsize : int, optionnel
            Taille de l'espace entre les bases.
        print_pos : bool, optionnel
            Si True, affiche les positions des bases.
        """
        seq=self.__rna.seq
        dotpar=self.__dotpar
        sep=" "*sepsize
        print(f"id: {self.__rna.id}")
        if print_pos:
            print("".join([str(i)+" "*(sepsize-len(str(i))+1) for i in range(len(seq))]))
        print(sep.join([b for b in seq]))
        print(sep.join([p for p in dotpar]))
        print(f"score: {self.__score}")
            
    #===================
    #Méthodes privées
    #===================        
    
    def __dot_par_to_bp(self,struc: str):
        """
        Convertit une structure en notation dot-parenthèse en une liste de paires de bases.
        
        Paramètres:
        -----------
        struc : str
            Structure en notation dot-parenthèse.
        
        Retourne:
        ---------
        tuple
            Liste triée des paires de bases sous forme de tuples.
        """
        open_parens = []
        bps = []
        if self.__parens_count(struc):        
            for i, x in enumerate(struc):
                if x == '(':
                    open_parens.append(i)
                elif x == ')':
                    if len(open_parens) > 0:
                        bps.append((open_parens.pop(), i))

        return tuple(sorted(bps))
    
    def __fold_to_dotpar(self):
        """
        Convertit une liste de paires de bases en notation dot-parenthèse.
        
        Retourne:
        ---------
        str
            Structure en notation dot-parenthèse.
        """
        rna=self.__rna.seq
        fold=self.__fold
        out_dot=['.']*len(rna)
        for i,j in fold:
            out_dot[i]="("
            out_dot[j]=")"
        dotpar_txt="".join(out_dot)
        return dotpar_txt
            
    def __parens_count(self,dotpar: str):
        """
        Vérifie que le nombre de parenthèses ouvrantes et fermantes est identique.
        
        Paramètres:
        -----------
        dotpar : str
            Structure en notation dot-parenthèse.
        
        Retourne:
        ---------
        bool
            True si le nombre de parenthèses ouvrantes et fermantes est identique, False sinon.
        """
        if not dotpar.count('(') == dotpar.count(')'):
            print("Erreur dans la structure: {}".format(dotpar))
            return False
        return True
    
    def __base_pairs_check(self):
        """
        Vérifie que toutes les paires de bases sont valides.
        
        Retourne:
        ---------
        bool
            True si toutes les paires de bases sont valides, False sinon.
        """
        fold=self.__fold
        seq=self.__rna.seq
        for bp in fold:
            bp_str = (seq[bp[0]] + seq[bp[1]]).upper()
            if bp_str not in self.__scores.allowed_bp:
                print("Paire de base non valide: {}".format(bp_str))
                return False
        return True
    
    def __search_link(self,rnastruct: str,i: int):
        """
        Trouve la paire de base correspondante pour une parenthèse ouvrante.
        
        Paramètres:
        -----------
        rnastruct : str
            Structure en notation dot-parenthèse.
        i : int
            Position de la parenthèse ouvrante.
        
        Retourne:
        ---------
        tuple
            Indices de la paire de base correspondante.
        """
        j=i+1
        while j<=len(rnastruct):
            if j-i>=1 and rnastruct[i]=="(" and rnastruct[i:j].count("(")==rnastruct[i:j].count(")") :
                return(i,j-1)
            j+=1
        return (i,i)

    def __arbre_s(self,structure: str,start=0,end=0,tree=None):
        """
        Construit un arbre représentant la structure ARN.
        
        Paramètres:
        -----------
        structure : str
            Structure en notation dot-parenthèse.
        start : int, optionnel
            Position de départ pour la construction de l'arbre.
        end : int, optionnel
            Position de fin pour la construction de l'arbre.
        tree : Tree, optionnel
            Arbre à compléter.
        
        Retourne:
        ---------
        Tree
            Arbre représentant la structure ARN.
        """
        if tree is None:
            start=0
            end=len(structure)
            tree=Tree((-1,-1))
        while start<end:
            start+=1
            car=structure[start-1]
            if car=='.':
                tree.add_child(Tree((start-1,start-1)))
            if car=='(':
                sub_struct=self.__search_link(structure,start-1)
                subtree=Tree(sub_struct)
                tree.add_child(subtree)
                self.__arbre_s(structure,sub_struct[0]+1,sub_struct[1],subtree)
                self.__arbre_s(structure,sub_struct[1],end,tree)
                break
        return tree
    
    def __construct_tree(self,tree: Tree,fold_dict: dict,start=0,end=0):
        """
        Construit un arbre à partir d'un dictionnaire de paires de bases (folds).
        
        Paramètres:
        -----------
        tree : Tree
            Arbre à compléter.
        fold_dict : dict
            Dictionnaire des paires de bases (fold).
        start : int, optionnel
            Position de départ pour la construction de l'arbre.
        end : int, optionnel
            Position de fin pour la construction de l'arbre.
        
        Retourne:
        ---------
        Tree
            Arbre représentant la structure ARN.
        """
        while start<end:
            if start not in fold_dict.keys():
                tree.add_child(Tree((start,start)))
                start+=1
            else:
                sub_struct=(start,fold_dict[start])
                subtree=Tree(sub_struct)
                tree.add_child(subtree)
                self.__construct_tree(subtree,fold_dict,sub_struct[0]+1,sub_struct[1])
                start=sub_struct[1]+1
        return tree