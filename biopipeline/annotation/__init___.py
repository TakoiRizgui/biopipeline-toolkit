"""
Module d'annotation génomique

Ce module contient des outils pour analyser des génomes annotés :
- EnzymeFinder : Identifier et cataloguer les enzymes industrielles
- ProkkaRunner : Automatiser l'annotation avec Prokka (à venir)
"""

from .enzyme_finder import EnzymeFinder, quick_enzyme_analysis

__all__ = ['EnzymeFinder', 'quick_enzyme_analysis']
