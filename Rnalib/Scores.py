#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 17:16:16 2023

@author: Mathieu Genete
"""

class Scores:
    """
    Classe pour gérer les scores entre différentes paires.
    
    Attributs:
    ----------
    __pairs : dict
        Dictionnaire contenant les scores pour chaque paire.
    __allowed_bp : list
        Liste des paires autorisées sous forme de chaînes de caractères.
    """
    def __init__(self,GC=3,AU=2,GU=1):
        """
        Initialise la classe Scores avec des valeurs par défaut pour GC, AU et GU.
        
        Paramètres:
        -----------
        GC : int, optionnel
            Score pour la paire GC (par défaut 3).
        AU : int, optionnel
            Score pour la paire AU (par défaut 2).
        GU : int, optionnel
            Score pour la paire GU (par défaut 1).
        
        Exceptions:
        -----------
        Exception
            Si les scores ne sont pas des entiers.
        """
        if not all([type(GC)==int,type(AU)==int,type(GU)==int]):
            raise Exception("Les scores doivent être des entiers")
        self.__pairs=None
        outscore={}    
        if GC<0: GC=0
        if AU<0: AU=0
        if AU<0: AU=0
        for b,s in {"GC":GC,"AU":AU,"GU":GU}.items():
            if s>0:
                outscore[(b[0],b[1])]=s
                outscore[(b[1],b[0])]=s
        self.__pairs=outscore
        self.__allowed_bp=["".join(k) for k in outscore.keys()]

    #===================
    #Getters Setters
    #===================
    
    @property
    def pairs(self):
        """
        Retourne le dictionnaire des paires et leurs scores.
        
        Retourne:
        ---------
        dict
            Dictionnaire des paires et leurs scores.
        """
        return self.__pairs
    
    @property
    def allowed_bp(self):
        """
        Retourne la liste des paires autorisées.
        
        Retourne:
        ---------
        list
            Liste des paires autorisées sous forme de chaînes de caractères.
        """
        return self.__allowed_bp
    
    @property
    def scores(self):
        """
        Retourne un dictionnaire des scores avec les paires triées.
        
        Retourne:
        ---------
        dict
            Dictionnaire des scores avec les paires triées.
        """
        return {"".join(sorted(v)):self.__pairs[v] for v in self.__pairs.keys()}

    #===================
    #Méthodes magiques
    #=================== 
    def __str__(self):
        """
        Retourne un dictionnaire contenant les scores pour chaque paire.
        
        Retourne:
        ---------
        dict
            Dictionnaire contenant les scores pour chaque paire.
        """
        return self.__pairs