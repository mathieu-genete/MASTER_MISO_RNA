#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 14:10:04 2023

@author: Mathieu Genete
"""
import copy

class Tree:
    """
    Classe représentant un arbre avec des nœuds pouvant avoir plusieurs enfants.
    """
    def __init__(self,value):
        """
        Initialise un nœud de l'arbre avec une valeur donnée.
        
        :param value: La valeur du nœud.
        """
        self.__valeur=value
        self.__childs=[]
        self.__parent=self
        
    def add_child(self,tree: 'Tree'):
        """
        Ajoute un enfant à ce nœud.
        
        :param tree: L'arbre enfant à ajouter.
        :return: L'arbre enfant ajouté.
        """
        self.__childs.append(tree)
        tree.set_parent(self)
        return tree

    #===================
    #Getters Setters
    #===================
    
    @property
    def valeur(self):
        """
        Retourne la valeur du nœud.
        
        :return: La valeur du nœud.
        """
        return self.__valeur
    
    @property
    def childs(self):
        """
        Retourne la liste des enfants du nœud.
        
        :return: La liste des enfants.
        """
        return self.__childs
    
    @property
    def parent(self):
        """
        Retourne le parent du nœud.
        
        :return: Le parent du nœud.
        """
        return self.__parent
    
    @property
    def is_leaf(self):
        """
        Vérifie si le nœud est une feuille (n'a pas d'enfants).
           
        :return: True si le nœud est une feuille, sinon False.
        """
        return len(self.__childs)==0
    
    @property
    def is_root(self):
        """
        Vérifie si le nœud est la racine (n'a pas de parent).
        
        :return: True si le nœud est la racine, sinon False.
        """
        return self.__parent==self
        
    def set_valeur(self,value):
        """
        Définit la valeur du nœud.
        
        :param value: La nouvelle valeur du nœud.
        """
        self.__valeur=value
        
    def set_childs(self,value: list):
        """
        Définit la liste des enfants du nœud.
        
        :param value: La nouvelle liste des enfants.
        """
        self.__childs=value

    def set_parent(self,value: 'Tree'):
        """
        Définit le parent du nœud.
        
        :param value: Le nouveau parent du nœud.
        """
        self.__parent=value

    #===================
    #Méthodes magiques
    #=================== 
    
    def __eq__(self, other: 'Tree'):
        """
        Vérifie si deux arbres ont la même architecture en comparant leur représentation en dotpar.
        
        :param other: L'autre arbre à comparer.
        :return: True si les arbres ont la même architecture, sinon False.
        """
        return self.tree_to_dotpar() == other.tree_to_dotpar()
    
    #===================
    #Méthodes publiques
    #===================
    def compact_tree(self):
        """
        Retourne une version compacte de l'arbre.
        
        :return: L'arbre compacté.
        """
        compacted_tree=copy.deepcopy(self)
        self.__vert_compact_tree(compacted_tree)    
        self.__horiz_compact_tree(compacted_tree)
        return compacted_tree

    def print_tree(self,print_tuples=True):
        """
        Affiche l'arbre.
        
        :param print_tuples: Indique si les tuples doivent être imprimés.
        """
        print(self.__print_tree(print_tuples=print_tuples))
        
    def tree_to_dotpar(self):
        """
        Convertit l'arbre en une représentation dot-bracket.
        
        :return: La représentation dot-bracket de l'arbre.
        """
        out_dotpar=[]
        stack=[(self,0)]
        while len(stack)>0:
            tree,position=stack.pop()
            for child in tree.childs:
                if len(child.childs)>0:
                    out_dotpar = out_dotpar[:position]+["(",")"]+out_dotpar[position:]
                    position+=2
                    stack.append((child,position-1))
                else:
                    out_dotpar = out_dotpar[:position]+["."]+out_dotpar[position:]
                    position+=1
        return "".join(out_dotpar)

    #===================
    #Méthodes privées
    #===================
    def __print_tree(self,tree=None, markerStr="+- ", levelMarkers=[],outstr="",print_tuples=True):
        """
        Affiche l'arbre avec une représentation graphique utilisant des marqueurs pour indiquer les niveaux.
        
        :param tree: L'arbre à imprimer. Si None, utilise l'arbre actuel.
        :param markerStr: La chaîne de caractères utilisée pour marquer les nœuds.
        :param levelMarkers: Les marqueurs de niveau pour l'indentation.
        :param outstr: La chaîne de caractères de sortie.
        :param print_tuples: Indique si les valeurs des nœuds doivent être imprimées sous forme de tuples.
        :return: La chaîne de caractères représentant l'arbre.
        """
        if tree is None:
            tree=self
        emptyStr = " "*len(markerStr)
        connectionStr = "|" + emptyStr[:-1]
        level = len(levelMarkers)
        mapper = lambda draw: connectionStr if draw else emptyStr
        markers = "".join(map(mapper, levelMarkers[:-1]))
        markers += markerStr if level > 0 else ""
        if print_tuples:
            valeur=str(tree.valeur)
        else:
            if tree.is_root:
                valeur="R"
            elif len(tree.childs)>0:
                valeur="N"
            else:
                valeur="L"
                
        outstr+=markers+valeur+"\n"
        for i, child in enumerate(tree.childs):
            isLast = i == len(tree.childs) - 1
            outstr=self.__print_tree(child, markerStr, [*levelMarkers, not isLast],outstr,print_tuples)
        return outstr
    
    def __print_simple_tree(self,tree=None, level=0,outstr=""):
        """
        Affiche l'arbre de manière simple avec une indentation pour chaque niveau.
        
        :param tree: L'arbre à imprimer. Si None, utilise l'arbre actuel.
        :param level: Le niveau actuel de l'arbre.
        :param outstr: La chaîne de caractères de sortie.
        :return: La chaîne de caractères représentant l'arbre.
        """
        if tree is None:
            tree=self
        outstr+="  " * level +" "+ str(tree.valeur) +"\n"
        for child in tree.childs:
            outstr=self.__print_simple_tree(child, level + 1,outstr)
        return outstr
        
    def __vert_compact_tree(self,tree: 'Tree'):
        """
        Compacte verticalement l'arbre en fusionnant les nœuds avec un seul enfant.
        
        :param tree: L'arbre à compacter.
        """
        for v in tree.childs:
            if len(v.childs)==1:
                v.set_valeur(v.childs[0].valeur)
                v.set_childs(v.childs[0].childs)
                for c in v.childs:
                    c.set_parent(v)
                self.__vert_compact_tree(v.parent)
            self.__vert_compact_tree(v)
            
    def __horiz_compact_tree(self,tree: 'Tree'):
        """
        Compacte horizontalement l'arbre en supprimant les feuilles consécutives.
        
        :param tree: L'arbre à compacter.
        """
        for i in range(1,len(tree.childs)):
            if tree.childs[i-1].is_leaf and tree.childs[i].is_leaf:
                tree.childs[i-1]=None
        tree.set_childs([j for j in tree.childs if j is not None])
        
        for v in tree.childs:
            self.__horiz_compact_tree(v)