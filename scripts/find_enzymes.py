#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from biopipeline.annotation import quick_enzyme_analysis

def main():
    parser = argparse.ArgumentParser(description='Identifier les enzymes industrielles dans un genome annote')
    parser.add_argument('genbank', help='Fichier GenBank annote (sortie Prokka)')
    parser.add_argument('--output', '-o', default='./enzyme_results', help='Dossier de sortie')
    
    args = parser.parse_args()
    
    try:
        quick_enzyme_analysis(args.genbank, args.output)
    except FileNotFoundError as e:
        print(f"\nERREUR : {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\nERREUR INATTENDUE : {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()