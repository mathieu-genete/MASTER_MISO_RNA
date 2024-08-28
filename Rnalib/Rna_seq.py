#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 17:08:39 2023

@author: Mathieu Genete
"""
from .Alphabet import Alphabet

class Rna_seq:
    """
    Classe représentant une séquence d'ARN.

    Attributs:
    ----------
    __id : str
        Identifiant de la séquence.
    __seq : str
        Séquence d'ARN en majuscules.
    """
    def __init__(self,seqid,inrna):
        """
        Initialise une nouvelle instance de la classe Rna_seq.

        Paramètres:
        -----------
        seqid : str
            Identifiant de la séquence.
        inrna : str
            Séquence d'ARN ou d'ADNc.
        """
        self.__id=seqid
        self.__seq=inrna.upper()
        self.__check_rna()

    #===================
    #Méthodes magiques
    #=================== 
    
    def __str__(self):
        """
        Retourne la séquence d'ARN sous forme de chaîne de caractères.

        Returns:
        ---------
        str
            La séquence d'ARN.
        """
        return self.__seq
    
    #===================
    #Getters Setters
    #===================
    
    @property
    def seq(self):
        """
        Retourne la séquence d'ARN.

        Returns:
        ---------
        str
            La séquence d'ARN.
        """
        return self.__seq
        
    @property
    def id(self):
        """
        Retourne l'identifiant de la séquence.

        Returns:
        ---------
        str
            L'identifiant de la séquence.
        """
        return self.__id
    
    #===================
    #Méthodes privées
    #===================
    def __check_rna(self):
        """
        Vérifie et convertit la séquence en ARN si nécessaire.

        Remplace les thymidines (T) par des uraciles (U) et vérifie que la séquence
        ne contient que des bases valides pour l'ARN.

        Raise:
        -----
        Exception
            Si la séquence contient des bases non valides pour l'ARN.
        """
        self.__seq = self.seq.replace("T","U")
        for b in self.seq:
            if b not in Alphabet.rna():
                raise Exception("La séquence n'est pas un ARN ou ADNc")