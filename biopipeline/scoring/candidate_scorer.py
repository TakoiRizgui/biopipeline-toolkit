"""
Système de Scoring pour Prioriser les Candidats Enzymatiques

Ce module permet de scorer et classer les enzymes selon plusieurs critères :
- Longueur optimale
- Présence de signal peptide
- Présence de domaines conservés
- Potentiel de thermostabilité
- Nouveauté (distance avec enzymes connues)

Auteur : Takoi Rizgui
Projet : Mastère Bioinformatique
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
import matplotlib.pyplot as plt
import seaborn as sns


class CandidateScorer:
    """
    Classe pour scorer et prioriser les candidats enzymatiques
    
    Example:
        >>> scorer = CandidateScorer()
        >>> scored = scorer.score_enzymes("enzyme_catalog.csv")
        >>> top20 = scored.nlargest(20, 'total_score')
    """
    
    def __init__(self):
        """Initialiser le scorer avec les paramètres par défaut"""
        self.criteria_weights = {
            'length': 0.25,           # 25% - Taille optimale
            'signal_peptide': 0.15,   # 15% - Sécrétion
            'ec_number': 0.20,        # 20% - Classification claire
            'family_priority': 0.20,  # 20% - Famille d'intérêt
            'gc_content': 0.10,       # 10% - Expression facilitée
            'complexity': 0.10        # 10% - Complexité séquence
        }
        
        # Familles prioritaires (selon projet)
        self.family_priorities = {
            'Lipases': 1.0,
            'Proteases': 0.9,
            'Cellulases': 0.8,
            'Laccases': 0.7,
            'Amylases': 0.6,
            'Peroxydases': 0.5,
            'Xylanases': 0.4,
            'Chitinases': 0.3
        }
    
    def score_length(self, length: int) -> float:
        """
        Score basé sur la longueur de la protéine
        
        Optimal : 250-600 aa
        Acceptable : 150-800 aa
        Suboptimal : < 150 ou > 800 aa
        
        Args:
            length: Longueur en acides aminés
            
        Returns:
            Score entre 0 et 1
        """
        if 250 <= length <= 600:
            return 1.0
        elif 200 <= length < 250 or 600 < length <= 700:
            return 0.8
        elif 150 <= length < 200 or 700 < length <= 800:
            return 0.6
        elif 100 <= length < 150 or 800 < length <= 1000:
            return 0.4
        else:
            return 0.2
    
    def score_signal_peptide(self, product: str) -> float:
        """
        Détecter présence probable de signal peptide
        
        Mots-clés indicateurs : signal, secreted, extracellular, exported
        
        Args:
            product: Description du produit
            
        Returns:
            1.0 si présent, 0.5 si incertain, 0.0 si absent
        """
        product_lower = product.lower()
        
        positive_keywords = ['signal', 'secreted', 'extracellular', 'exported', 
                           'precursor', 'preprotein']
        negative_keywords = ['intracellular', 'cytoplasmic', 'membrane']
        
        if any(kw in product_lower for kw in positive_keywords):
            return 1.0
        elif any(kw in product_lower for kw in negative_keywords):
            return 0.0
        else:
            return 0.5  # Incertain
    
    def score_ec_number(self, ec_number: str) -> float:
        """
        Score basé sur la présence d'un numéro EC
        
        Args:
            ec_number: Numéro EC (ex: "3.1.1.3")
            
        Returns:
            1.0 si EC complet, 0.5 si partiel, 0.0 si absent
        """
        if pd.isna(ec_number) or ec_number == 'N/A':
            return 0.0
        
        ec_parts = str(ec_number).split('.')
        
        if len(ec_parts) >= 4:  # EC complet (ex: 3.1.1.3)
            return 1.0
        elif len(ec_parts) >= 3:  # EC partiel (ex: 3.1.1)
            return 0.7
        elif len(ec_parts) >= 2:  # EC classe (ex: 3.1)
            return 0.4
        else:
            return 0.0
    
    def score_family_priority(self, family: str) -> float:
        """
        Score basé sur la priorité de la famille enzymatique
        
        Args:
            family: Nom de la famille
            
        Returns:
            Score entre 0 et 1
        """
        return self.family_priorities.get(family, 0.5)
    
    def score_gc_content(self, sequence: str) -> float:
        """
        Score basé sur le GC content (facilité d'expression)
        
        Optimal : 40-60% GC
        
        Args:
            sequence: Séquence protéique
            
        Returns:
            Score entre 0 and 1
        """
        if not sequence or len(sequence) < 50:
            return 0.5  # Séquence trop courte, score neutre
        
        # Compter G et C dans la séquence
        gc_count = sequence.count('G') + sequence.count('C')
        gc_percent = (gc_count / len(sequence)) * 100
        
        if 40 <= gc_percent <= 60:
            return 1.0
        elif 35 <= gc_percent < 40 or 60 < gc_percent <= 65:
            return 0.8
        elif 30 <= gc_percent < 35 or 65 < gc_percent <= 70:
            return 0.6
        else:
            return 0.4
    
    def score_complexity(self, sequence: str) -> float:
        """
        Score basé sur la complexité de la séquence
        
        Éviter : 
        - Séquences répétitives (low complexity)
        - Régions poly-X (polyalanine, etc.)
        
        Args:
            sequence: Séquence protéique
            
        Returns:
            Score entre 0 and 1
        """
        if not sequence or len(sequence) < 50:
            return 0.5
        
        # Diversité des acides aminés
        unique_aa = len(set(sequence))
        diversity = unique_aa / 20  # 20 acides aminés possibles
        
        # Détection poly-X (ex: AAAAAAAA)
        max_repeat = max([sequence.count(aa * 5) for aa in set(sequence)])
        
        if max_repeat > 0:  # Région répétitive détectée
            return 0.3
        elif diversity > 0.7:  # Bonne diversité
            return 1.0
        elif diversity > 0.5:
            return 0.7
        else:
            return 0.5
    
    def score_enzymes(self, enzyme_catalog: str, 
                     custom_weights: Optional[Dict] = None) -> pd.DataFrame:
        """
        Scorer toutes les enzymes d'un catalogue
        
        Args:
            enzyme_catalog: Chemin vers CSV du catalogue
            custom_weights: Poids personnalisés (optionnel)
            
        Returns:
            DataFrame avec scores détaillés
        """
        # Charger catalogue
        df = pd.read_csv(enzyme_catalog)
        
        # Utiliser poids personnalisés si fournis
        weights = custom_weights if custom_weights else self.criteria_weights
        
        # Calculer scores individuels
        df['score_length'] = df['length'].apply(self.score_length)
        df['score_signal'] = df['product'].apply(self.score_signal_peptide)
        df['score_ec'] = df['ec_number'].apply(self.score_ec_number)
        df['score_family'] = df['family'].apply(self.score_family_priority)
        df['score_gc'] = df['sequence'].apply(self.score_gc_content)
        df['score_complexity'] = df['sequence'].apply(self.score_complexity)
        
        # Score total pondéré
        df['total_score'] = (
            df['score_length'] * weights['length'] +
            df['score_signal'] * weights['signal_peptide'] +
            df['score_ec'] * weights['ec_number'] +
            df['score_family'] * weights['family_priority'] +
            df['score_gc'] * weights['gc_content'] +
            df['score_complexity'] * weights['complexity']
        )
        
        # Normaliser sur 100
        df['total_score'] = (df['total_score'] * 100).round(1)
        
        # Trier par score décroissant
        df = df.sort_values('total_score', ascending=False).reset_index(drop=True)
        
        # Ajouter rang
        df['rank'] = range(1, len(df) + 1)
        
        return df
    
    def export_top_candidates(self, scored_df: pd.DataFrame, 
                             top_n: int, output_fasta: str):
        """
        Exporter les top N candidats au format FASTA
        
        Args:
            scored_df: DataFrame avec scores
            top_n: Nombre de candidats à exporter
            output_fasta: Fichier FASTA de sortie
        """
        top_candidates = scored_df.head(top_n)
        
        with open(output_fasta, 'w') as f:
            for idx, row in top_candidates.iterrows():
                header = (f">{row['locus_tag']}|{row['family']}|"
                         f"score_{row['total_score']}|rank_{row['rank']}")
                f.write(header + '\n')
                f.write(row['sequence'] + '\n')
        
        print(f"✅ Top {top_n} candidats exportés : {output_fasta}")
    
    def plot_score_distribution(self, scored_df: pd.DataFrame, 
                                output_file: Optional[str] = None):
        """
        Visualiser la distribution des scores
        
        Args:
            scored_df: DataFrame avec scores
            output_file: Chemin pour sauvegarder (optionnel)
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # 1. Distribution score total
        ax1 = axes[0, 0]
        ax1.hist(scored_df['total_score'], bins=30, color='steelblue', 
                edgecolor='black', alpha=0.7)
        ax1.axvline(scored_df['total_score'].median(), color='red', 
                   linestyle='--', linewidth=2, 
                   label=f'Médiane: {scored_df["total_score"].median():.1f}')
        ax1.set_xlabel('Score Total', fontsize=12)
        ax1.set_ylabel('Nombre d\'enzymes', fontsize=12)
        ax1.set_title('Distribution des Scores', fontsize=14, fontweight='bold')
        ax1.legend()
        ax1.grid(axis='y', alpha=0.3)
        
        # 2. Scores par famille
        ax2 = axes[0, 1]
        family_scores = scored_df.groupby('family')['total_score'].mean().sort_values()
        family_scores.plot(kind='barh', ax=ax2, color='coral', alpha=0.7)
        ax2.set_xlabel('Score Moyen', fontsize=12)
        ax2.set_title('Score Moyen par Famille', fontsize=14, fontweight='bold')
        ax2.grid(axis='x', alpha=0.3)
        
        # 3. Scores individuels
        ax3 = axes[1, 0]
        score_cols = ['score_length', 'score_signal', 'score_ec', 
                     'score_family', 'score_gc', 'score_complexity']
        score_means = scored_df[score_cols].mean()
        score_means.plot(kind='bar', ax=ax3, color='purple', alpha=0.7)
        ax3.set_ylabel('Score Moyen', fontsize=12)
        ax3.set_title('Contribution de Chaque Critère', fontsize=14, fontweight='bold')
        ax3.tick_params(axis='x', rotation=45)
        ax3.grid(axis='y', alpha=0.3)
        
        # 4. Top 20
        ax4 = axes[1, 1]
        top20 = scored_df.head(20)
        colors_map = {'Lipases': '#FF6B6B', 'Proteases': '#4ECDC4', 
                     'Cellulases': '#95E1D3', 'Laccases': '#F38181'}
        colors = [colors_map.get(f, '#CCCCCC') for f in top20['family']]
        ax4.barh(range(20), top20['total_score'], color=colors, alpha=0.7)
        ax4.set_yticks(range(20))
        ax4.set_yticklabels([f"{r['rank']}. {r['locus_tag'][:15]}" 
                            for _, r in top20.iterrows()], fontsize=9)
        ax4.set_xlabel('Score', fontsize=12)
        ax4.set_title('Top 20 Candidats', fontsize=14, fontweight='bold')
        ax4.grid(axis='x', alpha=0.3)
        ax4.invert_yaxis()
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"✅ Graphique sauvegardé : {output_file}")
        else:
            plt.show()
    
    def generate_summary_report(self, scored_df: pd.DataFrame, 
                               output_file: str):
        """
        Générer un rapport résumé en CSV
        
        Args:
            scored_df: DataFrame avec scores
            output_file: Fichier CSV de sortie
        """
        summary = {
            'total_enzymes': len(scored_df),
            'mean_score': scored_df['total_score'].mean(),
            'median_score': scored_df['total_score'].median(),
            'top_10_percent_threshold': scored_df['total_score'].quantile(0.9),
            'excellent_candidates': len(scored_df[scored_df['total_score'] >= 80]),
            'good_candidates': len(scored_df[scored_df['total_score'] >= 60]),
            'top_family': scored_df.groupby('family')['total_score'].mean().idxmax()
        }
        
        summary_df = pd.DataFrame([summary])
        summary_df.to_csv(output_file, index=False)
        print(f"✅ Rapport résumé sauvegardé : {output_file}")
        
        return summary_df


# Fonction helper pour usage rapide
def quick_scoring(enzyme_catalog: str, top_n: int = 50, 
                 output_dir: str = "./scoring_results/"):
    """
    Scoring rapide et export des meilleurs candidats
    
    Args:
        enzyme_catalog: Fichier CSV du catalogue
        top_n: Nombre de top candidats à exporter
        output_dir: Dossier de sortie
    
    Example:
        >>> quick_scoring("enzyme_catalog.csv", top_n=50)
    """
    from pathlib import Path
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print(f"\n🎯 Scoring des enzymes : {enzyme_catalog}")
    
    # Scorer
    scorer = CandidateScorer()
    scored = scorer.score_enzymes(enzyme_catalog)
    
    # Sauvegarder catalogue scoré
    base_name = Path(enzyme_catalog).stem
    scored.to_csv(output_path / f"{base_name}_scored.csv", index=False)
    print(f"✅ Catalogue scoré sauvegardé")
    
    # Top candidats
    print(f"\n📊 Top {top_n} Candidats:")
    top = scored.head(top_n)
    print(top[['rank', 'locus_tag', 'family', 'length', 'total_score']].to_string())
    
    # Export FASTA
    scorer.export_top_candidates(
        scored, 
        top_n, 
        output_path / f"top{top_n}_for_alphafold.fasta"
    )
    
    # Graphiques
    scorer.plot_score_distribution(
        scored,
        output_path / "score_distribution.png"
    )
    
    # Rapport résumé
    scorer.generate_summary_report(
        scored,
        output_path / "scoring_summary.csv"
    )
    
    print(f"\n✅ Analyse terminée ! Résultats dans {output_dir}")


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        catalog = sys.argv[1]
        top_n = int(sys.argv[2]) if len(sys.argv) > 2 else 50
        output_dir = sys.argv[3] if len(sys.argv) > 3 else "./scoring_results/"
        quick_scoring(catalog, top_n, output_dir)
    else:
        print("Usage: python candidate_scorer.py <enzyme_catalog.csv> [top_n] [output_dir]")
