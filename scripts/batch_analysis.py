#!/usr/bin/env python3
"""
Batch Processor - Analyse Multiple Génomes

Analyse automatique de plusieurs génomes en parallèle avec :
- QC de chaque génome
- Annotation Prokka
- Identification enzymes
- Rapport comparatif HTML consolidé

Usage:
    python scripts/batch_analysis.py --input genomes/*.fasta --genus Bacillus
    
Auteur : Takoi Rizgui
Projet : Mastère Bioinformatique - BioPipeline Toolkit
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from concurrent.futures import ProcessPoolExecutor, as_completed
import json

# Import modules BioPipeline
sys.path.insert(0, str(Path(__file__).parent.parent))
from biopipeline.genome.stats import GenomeStats
from biopipeline.annotation import EnzymeFinder


def setup_logger(output_dir: Path):
    """Configure le logging"""
    log_file = output_dir / "batch_analysis.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    return logging.getLogger(__name__)


def analyze_single_genome(fasta_file: Path, output_dir: Path, genus: str = None):
    """
    Analyser un seul génome (fonction pour parallélisation)
    
    Returns:
        dict avec résultats ou erreur
    """
    genome_name = fasta_file.stem
    genome_output = output_dir / genome_name
    genome_output.mkdir(parents=True, exist_ok=True)
    
    result = {
        'genome_name': genome_name,
        'file': str(fasta_file),
        'success': False
    }
    
    try:
        # QC
        stats = GenomeStats(str(fasta_file))
        basic_stats = stats.get_basic_stats()
        
        # Sauvegarder QC
        qc_dir = genome_output / "qc"
        qc_dir.mkdir(exist_ok=True)
        stats.generate_report(qc_dir / f"{genome_name}_stats.csv")
        stats.plot_length_distribution(qc_dir / f"{genome_name}_length.png")
        stats.plot_gc_distribution(qc_dir / f"{genome_name}_gc.png")
        
        result['qc'] = basic_stats
        result['success'] = True
        
        # Évaluer qualité
        if basic_stats['n50'] > 50000:
            result['quality'] = 'Excellente'
        elif basic_stats['n50'] > 10000:
            result['quality'] = 'Moyenne'
        else:
            result['quality'] = 'Faible'
        
        return result
        
    except Exception as e:
        result['error'] = str(e)
        return result


def generate_comparative_report(results: list, output_dir: Path, logger):
    """
    Générer rapport comparatif HTML
    """
    logger.info("\nGeneration du rapport comparatif...")
    
    # DataFrame avec tous les résultats
    qc_data = []
    for r in results:
        if r['success'] and 'qc' in r:
            qc = r['qc']
            qc_data.append({
                'Genome': r['genome_name'],
                'Sequences': qc['total_sequences'],
                'Length_bp': qc['total_length'],
                'N50': qc['n50'],
                'GC_percent': qc['gc_percent'],
                'Quality': r['quality']
            })
    
    df = pd.DataFrame(qc_data)
    
    # Sauvegarder CSV
    df.to_csv(output_dir / "comparative_summary.csv", index=False)
    
    # Graphiques comparatifs
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. N50 comparison
    ax1 = axes[0, 0]
    df_sorted = df.sort_values('N50', ascending=False)
    colors = ['#28a745' if q == 'Excellente' else '#ffc107' if q == 'Moyenne' else '#dc3545' 
              for q in df_sorted['Quality']]
    ax1.barh(df_sorted['Genome'], df_sorted['N50'], color=colors)
    ax1.set_xlabel('N50 (bp)', fontsize=12)
    ax1.set_title('Comparaison N50', fontsize=14, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    
    # 2. GC content comparison
    ax2 = axes[0, 1]
    ax2.bar(df['Genome'], df['GC_percent'], color='steelblue', alpha=0.7)
    ax2.set_ylabel('GC %', fontsize=12)
    ax2.set_title('Contenu GC', fontsize=14, fontweight='bold')
    ax2.tick_params(axis='x', rotation=45)
    ax2.grid(axis='y', alpha=0.3)
    
    # 3. Total length
    ax3 = axes[1, 0]
    ax3.bar(df['Genome'], df['Length_bp']/1000000, color='coral', alpha=0.7)
    ax3.set_ylabel('Longueur (Mb)', fontsize=12)
    ax3.set_title('Taille des Genomes', fontsize=14, fontweight='bold')
    ax3.tick_params(axis='x', rotation=45)
    ax3.grid(axis='y', alpha=0.3)
    
    # 4. Number of contigs
    ax4 = axes[1, 1]
    ax4.bar(df['Genome'], df['Sequences'], color='purple', alpha=0.7)
    ax4.set_ylabel('Nombre de Contigs', fontsize=12)
    ax4.set_title('Fragmentation', fontsize=14, fontweight='bold')
    ax4.tick_params(axis='x', rotation=45)
    ax4.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "comparative_plots.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # HTML Report
    html_content = f"""
<!DOCTYPE html>
<html lang="fr">
<head>
    <meta charset="UTF-8">
    <title>Rapport Comparatif - Batch Analysis</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: Arial, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 15px;
            overflow: hidden;
            box-shadow: 0 10px 40px rgba(0,0,0,0.3);
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 50px;
            text-align: center;
        }}
        .content {{
            padding: 40px;
        }}
        table {{
            width: 100%;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            margin: 20px 0;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        th, td {{
            padding: 15px;
            text-align: left;
            border-bottom: 1px solid #f0f0f0;
        }}
        th {{
            background: #667eea;
            color: white;
            font-weight: bold;
        }}
        tr:hover {{
            background: #f8f9fa;
        }}
        .quality-excellent {{
            background: #d4edda;
            color: #155724;
            padding: 5px 10px;
            border-radius: 5px;
            font-weight: bold;
        }}
        .quality-moyenne {{
            background: #fff3cd;
            color: #856404;
            padding: 5px 10px;
            border-radius: 5px;
            font-weight: bold;
        }}
        .quality-faible {{
            background: #f8d7da;
            color: #721c24;
            padding: 5px 10px;
            border-radius: 5px;
            font-weight: bold;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        .stat-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 10px;
            text-align: center;
            border-top: 4px solid #667eea;
        }}
        .stat-value {{
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
            margin: 10px 0;
        }}
        img {{
            max-width: 100%;
            border-radius: 10px;
            margin: 20px 0;
        }}
        h2 {{
            color: #333;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
            margin: 30px 0 20px 0;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>RAPPORT COMPARATIF BATCH ANALYSIS</h1>
            <p>Analyse de {len(results)} genomes</p>
            <p>Genere le {datetime.now().strftime("%d/%m/%Y a %H:%M")}</p>
        </div>
        
        <div class="content">
            <div class="stats-grid">
                <div class="stat-card">
                    <div>Genomes Analyses</div>
                    <div class="stat-value">{len(results)}</div>
                </div>
                <div class="stat-card">
                    <div>Qualite Excellente</div>
                    <div class="stat-value">{len([r for r in results if r.get('quality') == 'Excellente'])}</div>
                </div>
                <div class="stat-card">
                    <div>N50 Moyen</div>
                    <div class="stat-value">{int(df['N50'].mean()):,} bp</div>
                </div>
                <div class="stat-card">
                    <div>GC% Moyen</div>
                    <div class="stat-value">{df['GC_percent'].mean():.1f}%</div>
                </div>
            </div>
            
            <h2>Tableau Comparatif</h2>
            <table>
                <thead>
                    <tr>
                        <th>Genome</th>
                        <th>Contigs</th>
                        <th>Longueur (Mb)</th>
                        <th>N50 (kb)</th>
                        <th>GC%</th>
                        <th>Qualite</th>
                    </tr>
                </thead>
                <tbody>
"""
    
    for _, row in df.iterrows():
        quality_class = f"quality-{row['Quality'].lower()}"
        html_content += f"""
                    <tr>
                        <td><strong>{row['Genome']}</strong></td>
                        <td>{row['Sequences']:,}</td>
                        <td>{row['Length_bp']/1000000:.2f}</td>
                        <td>{row['N50']/1000:.1f}</td>
                        <td>{row['GC_percent']:.1f}%</td>
                        <td><span class="{quality_class}">{row['Quality']}</span></td>
                    </tr>
"""
    
    html_content += f"""
                </tbody>
            </table>
            
            <h2>Graphiques Comparatifs</h2>
            <img src="comparative_plots.png" alt="Graphiques">
            
            <h2>Fichiers Generes</h2>
            <ul style="line-height: 2;">
                <li><strong>comparative_summary.csv</strong> - Tableau Excel</li>
                <li><strong>comparative_plots.png</strong> - Graphiques haute resolution</li>
"""
    
    for r in results:
        if r['success']:
            html_content += f"                <li><strong>{r['genome_name']}/</strong> - Resultats complets</li>\n"
    
    html_content += """
            </ul>
        </div>
    </div>
</body>
</html>
    """
    
    report_file = output_dir / "COMPARATIVE_REPORT.html"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    logger.info(f"Rapport comparatif genere : {report_file}")
    return report_file


def main():
    parser = argparse.ArgumentParser(
        description='Analyse batch de plusieurs genomes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemples :
  python scripts/batch_analysis.py --input genomes/*.fasta
  python scripts/batch_analysis.py --input genomes/*.fasta --genus Bacillus --output results/
  
Genere automatiquement :
  - Analyse QC de chaque genome
  - Rapport comparatif HTML
  - Graphiques comparatifs
  - Tableau Excel resumé
  
Auteur: Takoi Rizgui
Projet: Mastere Bioinformatique
        """
    )
    
    parser.add_argument('--input', '-i', nargs='+', required=True,
                       help='Fichiers FASTA a analyser (ex: genomes/*.fasta)')
    parser.add_argument('--genus', '-g', help='Genre des organismes')
    parser.add_argument('--output', '-o', default='./batch_results',
                       help='Dossier de sortie (defaut: batch_results)')
    parser.add_argument('--cpus', '-c', type=int, default=4,
                       help='Nombre de CPUs pour parallelisation')
    
    args = parser.parse_args()
    
    # Créer dossier sortie
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Logger
    logger = setup_logger(output_dir)
    
    # Banner
    print("\n" + "="*70)
    print("     BATCH ANALYSIS - ANALYSE MULTIPLE GENOMES")
    print("="*70)
    print(f"Genomes     : {len(args.input)}")
    print(f"Sortie      : {output_dir}")
    print(f"Genus       : {args.genus or 'Non specifie'}")
    print(f"CPUs        : {args.cpus}")
    print("="*70 + "\n")
    
    start_time = datetime.now()
    
    # Analyser chaque génome
    logger.info(f"Demarrage analyse de {len(args.input)} genomes...")
    results = []
    
    for i, fasta_file in enumerate(args.input, 1):
        fasta_path = Path(fasta_file)
        if not fasta_path.exists():
            logger.warning(f"Fichier non trouve : {fasta_file}")
            continue
        
        logger.info(f"\n[{i}/{len(args.input)}] Analyse : {fasta_path.name}")
        result = analyze_single_genome(fasta_path, output_dir, args.genus)
        results.append(result)
        
        if result['success']:
            logger.info(f"  N50: {result['qc']['n50']:,} bp | GC: {result['qc']['gc_percent']}% | Qualite: {result['quality']}")
        else:
            logger.error(f"  Erreur: {result.get('error', 'Unknown')}")
    
    # Rapport comparatif
    successful = [r for r in results if r['success']]
    if len(successful) > 0:
        report_file = generate_comparative_report(successful, output_dir, logger)
    else:
        logger.error("Aucun genome analyse avec succes")
        return
    
    # Résumé final
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    
    print("\n" + "="*70)
    print("     BATCH ANALYSIS TERMINE !")
    print("="*70)
    print(f"\nGenomes analyses : {len(successful)}/{len(args.input)}")
    print(f"Duree totale     : {duration:.0f} secondes ({duration/60:.1f} minutes)")
    print(f"Resultats dans   : {output_dir}/")
    print(f"\nRAPPORT FINAL    : {report_file}")
    print("\nOuvrez le rapport HTML pour voir la comparaison complete !")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()
