#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 17:07:57 2023

@author: Mathieu Genete
"""

class Alphabet:
    """
    Classe contenant des méthodes statiques pour retourner des alphabets spécifiques.
    """

    @staticmethod
    def rna():
        """
        Retourne l'alphabet de l'ARN.

        Returns:
            str: Une chaîne de caractères représentant les bases de l'ARN (A, C, G, U).
        """
        return "ACGU"

    @staticmethod
    def dna():
        """
        Retourne l'alphabet de l'ADN.

        Returns:
            str: Une chaîne de caractères représentant les bases de l'ADN (A, G, C, T).
        """
        return "AGCT"

    @staticmethod
    def dotpar():
        """
        Retourne les caractères utilisés pour la notation en parenthèses.

        Returns:
            str: Une chaîne de caractères représentant la notation en parenthèses (., ).
        """
        return "(.)"
