from Bio import SeqIO
import pandas as pd
from pathlib import Path

class EnzymeFinder:
    ENZYME_FAMILIES = {
        'Lipases': {'keywords': ['lipase', 'esterase'], 'color': '#FF6B6B'},
        'Proteases': {'keywords': ['protease', 'peptidase'], 'color': '#4ECDC4'},
        'Cellulases': {'keywords': ['cellulase', 'glucosidase'], 'color': '#95E1D3'},
    }
    
    def __init__(self, genbank_file):
        self.genbank_file = Path(genbank_file)
        if not self.genbank_file.exists():
            raise FileNotFoundError(f"Fichier non trouve: {genbank_file}")
        self.records = list(SeqIO.parse(str(self.genbank_file), "genbank"))
        self.enzymes = None
    
    def find_all_enzymes(self):
        enzymes_list = []
        for record in self.records:
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                product = feature.qualifiers.get('product', [''])[0]
                family = self._classify_enzyme(product)
                if family:
                    enzymes_list.append({
                        'locus_tag': feature.qualifiers.get('locus_tag', [''])[0],
                        'product': product,
                        'family': family,
                        'length': len(feature.qualifiers.get('translation', [''])[0]),
                        'sequence': feature.qualifiers.get('translation', [''])[0]
                    })
        self.enzymes = pd.DataFrame(enzymes_list)
        return self.enzymes
    
    def _classify_enzyme(self, product):
        product_lower = product.lower()
        for family, info in self.ENZYME_FAMILIES.items():
            for keyword in info['keywords']:
                if keyword in product_lower:
                    return family
        return None
    
    def export_to_csv(self, output_file):
        if self.enzymes is not None and len(self.enzymes) > 0:
            export_df = self.enzymes.drop('sequence', axis=1)
            export_df.to_csv(output_file, index=False)
            print(f"Catalogue exporte : {output_file}")
    
    def __str__(self):
        if self.enzymes is None:
            return "EnzymeFinder (non analyse)"
        return f"\nTotal enzymes : {len(self.enzymes)}\n"

def quick_enzyme_analysis(genbank_file, output_dir="./enzyme_results/"):
    from pathlib import Path
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print(f"\nAnalyse des enzymes : {genbank_file}")
    finder = EnzymeFinder(genbank_file)
    enzymes = finder.find_all_enzymes()
    
    if len(enzymes) == 0:
        print("Aucune enzyme trouvee")
        return
    
    print(finder)
    base_name = Path(genbank_file).stem
    finder.export_to_csv(output_path / f"{base_name}_enzymes_catalog.csv")
    print(f"\nAnalyse terminee ! Resultats dans {output_dir}")