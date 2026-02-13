#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from datetime import datetime
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))
from biopipeline.genome.stats import GenomeStats

def setup_logger(output_dir, genome_name):
    log_file = output_dir / f"{genome_name}_analysis.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def create_html_report(stats_df, genome_name, output_dir):
    stats = stats_df.iloc[0].to_dict()
    
    html_content = f"""<!DOCTYPE html>
<html lang="fr">
<head>
    <meta charset="UTF-8">
    <title>Rapport - {genome_name}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: #f4f4f4;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            border-radius: 10px;
            text-align: center;
            margin-bottom: 30px;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        .stat-card {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .stat-label {{
            color: #666;
            font-size: 0.9em;
            text-transform: uppercase;
        }}
        .stat-value {{
            color: #667eea;
            font-size: 2em;
            font-weight: bold;
            margin-top: 10px;
        }}
        .image-container {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            margin: 20px 0;
            text-align: center;
        }}
        .image-container img {{
            max-width: 100%;
            border-radius: 5px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 Rapport d'Analyse Génomique</h1>
        <p>{genome_name}</p>
        <p>Généré le {datetime.now().strftime("%d/%m/%Y à %H:%M")}</p>
    </div>
    
    <h2>📊 Statistiques Clés</h2>
    <div class="stats-grid">
        <div class="stat-card">
            <div class="stat-label">Nombre de Contigs</div>
            <div class="stat-value">{stats['total_sequences']:,}</div>
        </div>
        <div class="stat-card">
            <div class="stat-label">Longueur Totale</div>
            <div class="stat-value">{stats['total_length']:,} pb</div>
        </div>
        <div class="stat-card">
            <div class="stat-label">N50</div>
            <div class="stat-value">{stats['n50']:,} pb</div>
        </div>
        <div class="stat-card">
            <div class="stat-label">Contenu GC</div>
            <div class="stat-value">{stats['gc_percent']}%</div>
        </div>
    </div>
    
    <h2>📈 Visualisations</h2>
    <div class="image-container">
        <h3>Distribution des Longueurs</h3>
        <img src="{genome_name}_length_dist.png" alt="Longueurs">
    </div>
    <div class="image-container">
        <h3>Analyse du Contenu GC</h3>
        <img src="{genome_name}_gc_dist.png" alt="GC">
    </div>
    
    <div style="background: white; padding: 20px; border-radius: 8px; margin-top: 30px; text-align: center;">
        <p><strong>BioPipeline Toolkit v0.1.0</strong></p>
        <p>Takoi Rizgui - Mastère Bioinformatique</p>
    </div>
</body>
</html>"""
    
    html_file = output_dir / f"{genome_name}_report.html"
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    return html_file

def main():
    parser = argparse.ArgumentParser(description='Analyse rapide d\'un génome assemblé')
    parser.add_argument('fasta', help='Fichier FASTA du génome')
    parser.add_argument('--output', '-o', default='./results', help='Dossier de sortie')
    parser.add_argument('--min-length', '-m', type=int, default=0, help='Longueur minimale des contigs')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    genome_name = Path(args.fasta).stem
    logger = setup_logger(output_dir, genome_name)
    
    logger.info("="*60)
    logger.info("ANALYSE GÉNOMIQUE - BioPipeline Toolkit")
    logger.info("="*60)
    logger.info(f"Fichier : {args.fasta}")
    logger.info(f"Sortie : {output_dir}")
    
    try:
        logger.info("\n📊 Chargement du génome...")
        stats = GenomeStats(args.fasta)
        logger.info(f"✅ {len(stats.sequences)} séquences chargées")
        
        logger.info("\n📈 Calcul des statistiques...")
        stats_df = stats.generate_report(output_dir / f"{genome_name}_stats.csv")
        logger.info("✅ Statistiques sauvegardées")
        
        logger.info("\n" + str(stats))
        
        logger.info("\n🎨 Génération des graphiques...")
        stats.plot_length_distribution(
            output_dir / f"{genome_name}_length_dist.png",
            min_length=args.min_length
        )
        stats.plot_gc_distribution(output_dir / f"{genome_name}_gc_dist.png")
        logger.info("✅ Graphiques générés")
        
        logger.info("\n📄 Génération du rapport HTML...")
        html_file = create_html_report(stats_df, genome_name, output_dir)
        logger.info(f"✅ Rapport : {html_file}")
        
        logger.info("\n" + "="*60)
        logger.info("✅ ANALYSE TERMINÉE")
        logger.info("="*60)
        logger.info(f"\n📁 Fichiers générés :")
        logger.info(f"  • {genome_name}_stats.csv")
        logger.info(f"  • {genome_name}_length_dist.png")
        logger.info(f"  • {genome_name}_gc_dist.png")
        logger.info(f"  • {genome_name}_report.html  ← Ouvre ce fichier !")
        logger.info(f"  • {genome_name}_analysis.log")
        
    except Exception as e:
        logger.error(f"❌ ERREUR : {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()