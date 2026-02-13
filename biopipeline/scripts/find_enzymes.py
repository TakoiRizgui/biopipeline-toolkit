#!/usr/bin/env python3
"""
Script rapide pour identifier les enzymes dans un génome annoté

Usage:
    python scripts/find_enzymes.py genome.gbk --output enzyme_results/
    
Auteur : Takoi Rizgui
Projet : Mastère Bioinformatique - Annotation Génomique
"""

import argparse
import sys
from pathlib import Path

# Import du module
sys.path.insert(0, str(Path(__file__).parent.parent))
from biopipeline.annotation import quick_enzyme_analysis


def main():
    parser = argparse.ArgumentParser(
        description='Identifier les enzymes industrielles dans un génome annoté',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemples :
  python scripts/find_enzymes.py prokka_output.gbk
  python scripts/find_enzymes.py prokka_output.gbk --output results/
  
Le script génère automatiquement :
  - Catalogue CSV de toutes les enzymes
  - Fichier FASTA pour AlphaFold
  - Graphiques (camembert, histogrammes)
  - Rapport HTML interactif
  
Familles d'enzymes identifiées :
  • Lipases / Estérases
  • Protéases / Peptidases
  • Cellulases / Glucosidases
  • Laccases / Oxydoréductases
  • Amylases
  • Peroxydases
  • Xylanases
  • Chitinases

Auteur: Takoi Rizgui
Projet: Mastère Bioinformatique
        """
    )
    
    parser.add_argument('genbank', help='Fichier GenBank annoté (sortie Prokka)')
    parser.add_argument('--output', '-o', default='./enzyme_results',
                       help='Dossier de sortie (défaut: ./enzyme_results)')
    
    args = parser.parse_args()
    
    try:
        quick_enzyme_analysis(args.genbank, args.output)
    except FileNotFoundError as e:
        print(f"\n❌ ERREUR : {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ ERREUR INATTENDUE : {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
