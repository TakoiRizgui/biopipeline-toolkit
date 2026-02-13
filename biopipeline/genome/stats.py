"""
Module pour calculer les statistiques génomiques
"""

from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import List, Dict, Tuple


class GenomeStats:
    """
    Classe pour analyser la qualité d'un génome assemblé
    
    Attributes:
        fasta_file (str): Chemin vers le fichier FASTA
        sequences (list): Liste des séquences du génome
    
    Example:
        >>> stats = GenomeStats("assembly.fasta")
        >>> n50 = stats.calculate_n50()
        >>> report = stats.generate_report()
        >>> report.to_csv("genome_stats.csv")
    """
    
    def __init__(self, fasta_file: str):
        """
        Initialiser avec un fichier FASTA
        
        Args:
            fasta_file: Chemin vers le fichier FASTA d'assemblage
        """
        self.fasta_file = Path(fasta_file)
        if not self.fasta_file.exists():
            raise FileNotFoundError(f"Fichier non trouvé: {fasta_file}")
        
        self.sequences = list(SeqIO.parse(str(self.fasta_file), "fasta"))
        if not self.sequences:
            raise ValueError(f"Aucune séquence trouvée dans {fasta_file}")
    
    def calculate_n50(self) -> int:
        """
        Calculer le N50 (longueur médiane pondérée des contigs)
        
        Returns:
            int: Valeur du N50 en bases
        """
        lengths = sorted([len(seq) for seq in self.sequences], reverse=True)
        total_length = sum(lengths)
        cumulative = 0
        
        for length in lengths:
            cumulative += length
            if cumulative >= total_length / 2:
                return length
        
        return 0
    
    def calculate_n90(self) -> int:
        """
        Calculer le N90
        
        Returns:
            int: Valeur du N90 en bases
        """
        lengths = sorted([len(seq) for seq in self.sequences], reverse=True)
        total_length = sum(lengths)
        cumulative = 0
        
        for length in lengths:
            cumulative += length
            if cumulative >= total_length * 0.9:
                return length
        
        return 0
    
    def gc_content(self) -> float:
        """
        Calculer le contenu GC global du génome
        
        Returns:
            float: Pourcentage GC (0-100)
        """
        total_gc = 0
        total_bases = 0
        
        for seq in self.sequences:
            sequence = str(seq.seq).upper()
            total_gc += sequence.count('G') + sequence.count('C')
            total_bases += len(sequence)
        
        if total_bases == 0:
            return 0.0
        
        return (total_gc / total_bases) * 100
    
    def gc_content_per_sequence(self) -> pd.DataFrame:
        """
        Calculer le GC% pour chaque séquence
        
        Returns:
            DataFrame avec colonnes: sequence_id, length, gc_percent
        """
        data = []
        
        for seq in self.sequences:
            sequence = str(seq.seq).upper()
            gc = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
            data.append({
                'sequence_id': seq.id,
                'length': len(seq),
                'gc_percent': gc
            })
        
        return pd.DataFrame(data)
    
    def get_basic_stats(self) -> Dict:
        """
        Obtenir toutes les statistiques de base
        
        Returns:
            dict: Dictionnaire avec toutes les stats
        """
        lengths = [len(seq) for seq in self.sequences]
        
        return {
            'genome_file': self.fasta_file.name,
            'total_sequences': len(self.sequences),
            'total_length': sum(lengths),
            'n50': self.calculate_n50(),
            'n90': self.calculate_n90(),
            'gc_percent': round(self.gc_content(), 2),
            'longest_contig': max(lengths) if lengths else 0,
            'shortest_contig': min(lengths) if lengths else 0,
            'mean_length': round(sum(lengths) / len(lengths), 2) if lengths else 0,
            'median_length': sorted(lengths)[len(lengths)//2] if lengths else 0
        }
    
    def generate_report(self, output_csv: str = None) -> pd.DataFrame:
        """
        Générer un rapport complet sous forme de DataFrame
        
        Args:
            output_csv: Chemin optionnel pour sauvegarder en CSV
        
        Returns:
            DataFrame avec les statistiques
        """
        stats = self.get_basic_stats()
        df = pd.DataFrame([stats])
        
        if output_csv:
            df.to_csv(output_csv, index=False)
            print(f"✅ Rapport sauvegardé: {output_csv}")
        
        return df
    
    def plot_length_distribution(self, output_file: str = None, min_length: int = 0):
        """
        Créer un histogramme de la distribution des longueurs de contigs
        
        Args:
            output_file: Chemin pour sauvegarder le graphique
            min_length: Longueur minimale pour filtrer les petits contigs
        """
        lengths = [len(seq) for seq in self.sequences if len(seq) >= min_length]
        
        plt.figure(figsize=(10, 6))
        plt.hist(lengths, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
        plt.xlabel('Longueur des contigs (pb)', fontsize=12)
        plt.ylabel('Nombre de contigs', fontsize=12)
        plt.title(f'Distribution des longueurs de contigs\n{self.fasta_file.name}', 
                  fontsize=14, fontweight='bold')
        plt.grid(axis='y', alpha=0.3)
        
        # Ajouter statistiques sur le graphique
        stats = self.get_basic_stats()
        textstr = f"N50: {stats['n50']:,} pb\nTotal: {stats['total_sequences']} contigs"
        plt.text(0.98, 0.97, textstr, transform=plt.gca().transAxes,
                fontsize=10, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"✅ Graphique sauvegardé: {output_file}")
        else:
            plt.show()
    
    def plot_gc_distribution(self, output_file: str = None):
        """
        Créer un graphique du contenu GC par contig
        
        Args:
            output_file: Chemin pour sauvegarder le graphique
        """
        gc_data = self.gc_content_per_sequence()
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # Histogramme GC%
        ax1.hist(gc_data['gc_percent'], bins=30, color='coral', 
                edgecolor='black', alpha=0.7)
        ax1.axvline(gc_data['gc_percent'].mean(), color='red', 
                   linestyle='--', linewidth=2, label=f'Moyenne: {gc_data["gc_percent"].mean():.1f}%')
        ax1.set_xlabel('Contenu GC (%)', fontsize=12)
        ax1.set_ylabel('Nombre de contigs', fontsize=12)
        ax1.set_title('Distribution du contenu GC', fontsize=14, fontweight='bold')
        ax1.legend()
        ax1.grid(axis='y', alpha=0.3)
        
        # Scatter GC% vs longueur
        ax2.scatter(gc_data['length'], gc_data['gc_percent'], 
                   alpha=0.5, color='steelblue', s=20)
        ax2.set_xlabel('Longueur du contig (pb)', fontsize=12)
        ax2.set_ylabel('Contenu GC (%)', fontsize=12)
        ax2.set_title('GC% vs Longueur', fontsize=14, fontweight='bold')
        ax2.set_xscale('log')
        ax2.grid(alpha=0.3)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"✅ Graphique sauvegardé: {output_file}")
        else:
            plt.show()
    
    def __str__(self) -> str:
        """Représentation textuelle des stats"""
        stats = self.get_basic_stats()
        return f"""
Statistiques Génomiques - {stats['genome_file']}
{'='*50}
Nombre de séquences : {stats['total_sequences']:,}
Longueur totale     : {stats['total_length']:,} pb
N50                 : {stats['n50']:,} pb
N90                 : {stats['n90']:,} pb
GC%                 : {stats['gc_percent']}%
Plus long contig    : {stats['longest_contig']:,} pb
Plus court contig   : {stats['shortest_contig']:,} pb
Longueur moyenne    : {stats['mean_length']:,} pb
        """


# Fonction helper pour analyse rapide
def quick_stats(fasta_file: str, output_dir: str = "./"):
    """
    Analyse rapide d'un génome et génération de tous les rapports
    
    Args:
        fasta_file: Fichier FASTA à analyser
        output_dir: Dossier pour sauvegarder les résultats
    
    Example:
        >>> quick_stats("assembly.fasta", "./results/")
    """
    from pathlib import Path
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Analyser
    stats = GenomeStats(fasta_file)
    
    # Rapport texte
    print(stats)
    
    # Rapport CSV
    base_name = Path(fasta_file).stem
    stats.generate_report(output_path / f"{base_name}_stats.csv")
    
    # Graphiques
    stats.plot_length_distribution(output_path / f"{base_name}_length_dist.png")
    stats.plot_gc_distribution(output_path / f"{base_name}_gc_dist.png")
    
    print(f"\n✅ Analyse terminée ! Résultats dans {output_dir}")


if __name__ == "__main__":
    # Exemple d'utilisation
    import sys
    
    if len(sys.argv) > 1:
        fasta_file = sys.argv[1]
        output_dir = sys.argv[2] if len(sys.argv) > 2 else "./"
        quick_stats(fasta_file, output_dir)
    else:
        print("Usage: python stats.py <genome.fasta> [output_dir]")
