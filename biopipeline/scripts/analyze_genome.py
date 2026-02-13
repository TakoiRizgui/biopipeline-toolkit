#!/usr/bin/env python3
"""
Script d'analyse rapide d'un génome assemblé

Usage:
    python scripts/analyze_genome.py genome.fasta --output results/
    
Génère automatiquement :
    - Statistiques CSV
    - Graphiques de qualité (PNG)
    - Rapport HTML
    - Fichier log

Auteur : Takoi Rizgui
Projet : Mastère Bioinformatique - Annotation Génomique
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime
import logging

# Import du module
sys.path.insert(0, str(Path(__file__).parent.parent))
from biopipeline.genome.stats import GenomeStats


def setup_logger(output_dir: Path, genome_name: str) -> logging.Logger:
    """Configure le système de logging"""
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


def create_html_report(stats_df, genome_name: str, output_dir: Path):
    """Créer un rapport HTML professionnel"""
    
    stats = stats_df.iloc[0].to_dict()
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="fr">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Rapport d'Analyse - {genome_name}</title>
        <style>
            * {{
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }}
            
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                line-height: 1.6;
                color: #333;
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                padding: 20px;
            }}
            
            .container {{
                max-width: 1200px;
                margin: 0 auto;
                background: white;
                border-radius: 15px;
                box-shadow: 0 10px 40px rgba(0,0,0,0.2);
                overflow: hidden;
            }}
            
            .header {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 40px;
                text-align: center;
            }}
            
            .header h1 {{
                font-size: 2.5em;
                margin-bottom: 10px;
            }}
            
            .header .subtitle {{
                font-size: 1.1em;
                opacity: 0.9;
            }}
            
            .content {{
                padding: 40px;
            }}
            
            .info-section {{
                background: #f8f9fa;
                padding: 20px;
                border-radius: 10px;
                margin-bottom: 30px;
                border-left: 5px solid #667eea;
            }}
            
            .stats-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
                gap: 20px;
                margin: 30px 0;
            }}
            
            .stat-card {{
                background: white;
                padding: 25px;
                border-radius: 10px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
                border-top: 3px solid #667eea;
                transition: transform 0.3s ease;
            }}
            
            .stat-card:hover {{
                transform: translateY(-5px);
                box-shadow: 0 5px 20px rgba(0,0,0,0.15);
            }}
            
            .stat-label {{
                color: #666;
                font-size: 0.9em;
                text-transform: uppercase;
                letter-spacing: 1px;
                margin-bottom: 10px;
            }}
            
            .stat-value {{
                color: #667eea;
                font-size: 2em;
                font-weight: bold;
            }}
            
            .images-section {{
                margin: 40px 0;
            }}
            
            .image-container {{
                margin: 20px 0;
                text-align: center;
            }}
            
            .image-container img {{
                max-width: 100%;
                border-radius: 10px;
                box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            }}
            
            .image-title {{
                font-size: 1.2em;
                margin: 15px 0;
                color: #333;
                font-weight: 600;
            }}
            
            .footer {{
                background: #f8f9fa;
                padding: 30px;
                text-align: center;
                color: #666;
                border-top: 1px solid #ddd;
            }}
            
            .badge {{
                display: inline-block;
                padding: 5px 15px;
                background: #667eea;
                color: white;
                border-radius: 20px;
                font-size: 0.9em;
                margin: 5px;
            }}
            
            .quality-indicator {{
                padding: 10px 20px;
                border-radius: 5px;
                display: inline-block;
                font-weight: bold;
                margin: 10px 0;
            }}
            
            .quality-good {{
                background: #d4edda;
                color: #155724;
            }}
            
            .quality-warning {{
                background: #fff3cd;
                color: #856404;
            }}
            
            .quality-poor {{
                background: #f8d7da;
                color: #721c24;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>🧬 Rapport d'Analyse Génomique</h1>
                <div class="subtitle">
                    {genome_name}<br>
                    Généré le {datetime.now().strftime("%d/%m/%Y à %H:%M")}
                </div>
            </div>
            
            <div class="content">
                <div class="info-section">
                    <h2>📋 Informations Générales</h2>
                    <p style="margin-top: 10px;">
                        <span class="badge">Fichier : {stats['genome_file']}</span>
                        <span class="badge">Analyse : BioPipeline Toolkit v0.1</span>
                    </p>
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
                    
                    <div class="stat-card">
                        <div class="stat-label">Plus Long Contig</div>
                        <div class="stat-value">{stats['longest_contig']:,} pb</div>
                    </div>
                    
                    <div class="stat-card">
                        <div class="stat-label">Longueur Moyenne</div>
                        <div class="stat-value">{stats['mean_length']:,.0f} pb</div>
                    </div>
                </div>
                
                <div class="info-section">
                    <h2>✅ Évaluation de la Qualité</h2>
                    <div style="margin-top: 15px;">
                        {'<div class="quality-indicator quality-good">✅ Excellent : N50 > 50 kb</div>' if stats['n50'] > 50000 else 
                         '<div class="quality-indicator quality-warning">⚠️ Moyen : N50 entre 10-50 kb</div>' if stats['n50'] > 10000 else
                         '<div class="quality-indicator quality-poor">❌ Fragmented : N50 < 10 kb</div>'}
                        
                        <p style="margin-top: 10px; color: #666;">
                            Le N50 de {stats['n50']:,} pb indique 
                            {'une bonne qualité d\'assemblage' if stats['n50'] > 50000 else 
                             'une qualité d\'assemblage moyenne' if stats['n50'] > 10000 else
                             'un assemblage fragmenté'}.
                        </p>
                    </div>
                </div>
                
                <div class="images-section">
                    <h2>📈 Visualisations</h2>
                    
                    <div class="image-container">
                        <h3 class="image-title">Distribution des Longueurs de Contigs</h3>
                        <img src="{genome_name}_length_dist.png" alt="Distribution des longueurs">
                    </div>
                    
                    <div class="image-container">
                        <h3 class="image-title">Analyse du Contenu GC</h3>
                        <img src="{genome_name}_gc_dist.png" alt="Distribution GC">
                    </div>
                </div>
                
                <div class="info-section">
                    <h2>ℹ️ Notes Méthodologiques</h2>
                    <ul style="margin-top: 10px; line-height: 2;">
                        <li><strong>N50</strong> : Longueur médiane pondérée (50% des bases dans contigs ≥ N50)</li>
                        <li><strong>Contenu GC</strong> : Pourcentage de guanine + cytosine dans le génome</li>
                        <li><strong>Analyse</strong> : Effectuée avec BioPython et BioPipeline Toolkit</li>
                    </ul>
                </div>
            </div>
            
            <div class="footer">
                <p><strong>BioPipeline Toolkit v0.1.0</strong></p>
                <p>Développé par Takoi Rizgui - Mastère Bioinformatique</p>
                <p>Projet Tuniso-Italien - Annotation Génomique & Enzymes Industrielles</p>
                <p style="margin-top: 10px; font-size: 0.9em; color: #999;">
                    🔗 <a href="https://github.com/TakoiRizgui/biopipeline-toolkit" style="color: #667eea;">
                        github.com/TakoiRizgui/biopipeline-toolkit
                    </a>
                </p>
            </div>
        </div>
    </body>
    </html>
    """
    
    html_file = output_dir / f"{genome_name}_report.html"
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    return html_file


def main():
    """Fonction principale"""
    
    parser = argparse.ArgumentParser(
        description='Analyse rapide d\'un génome assemblé',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemples d'utilisation:
  python scripts/analyze_genome.py genome.fasta
  python scripts/analyze_genome.py genome.fasta --output results/
  python scripts/analyze_genome.py genome.fasta --output results/ --min-length 1000

Auteur: Takoi Rizgui
Projet: Mastère Bioinformatique - Annotation Génomique
        """
    )
    
    parser.add_argument('fasta', help='Fichier FASTA du génome assemblé')
    parser.add_argument('--output', '-o', default='./results',
                       help='Dossier de sortie (défaut: ./results)')
    parser.add_argument('--min-length', '-m', type=int, default=0,
                       help='Longueur minimale des contigs pour les graphiques (défaut: 0)')
    
    args = parser.parse_args()
    
    # Créer dossier de sortie
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Nom de base du génome
    genome_name = Path(args.fasta).stem
    
    # Logger
    logger = setup_logger(output_dir, genome_name)
    
    logger.info("="*60)
    logger.info("ANALYSE GÉNOMIQUE - BioPipeline Toolkit")
    logger.info("="*60)
    logger.info(f"Fichier d'entrée : {args.fasta}")
    logger.info(f"Dossier de sortie : {output_dir}")
    logger.info(f"Longueur minimale : {args.min_length} pb")
    logger.info("")
    
    try:
        # 1. Charger et analyser le génome
        logger.info("📊 Chargement du génome...")
        stats = GenomeStats(args.fasta)
        logger.info(f"✅ {len(stats.sequences)} séquences chargées")
        
        # 2. Générer statistiques
        logger.info("\n📈 Calcul des statistiques...")
        stats_df = stats.generate_report(
            output_dir / f"{genome_name}_stats.csv"
        )
        logger.info("✅ Statistiques calculées et sauvegardées")
        
        # 3. Afficher résumé
        logger.info("\n" + "="*60)
        logger.info("RÉSUMÉ DES STATISTIQUES")
        logger.info("="*60)
        logger.info(str(stats))
        
        # 4. Générer graphiques
        logger.info("🎨 Génération des graphiques...")
        stats.plot_length_distribution(
            output_dir / f"{genome_name}_length_dist.png",
            min_length=args.min_length
        )
        stats.plot_gc_distribution(
            output_dir / f"{genome_name}_gc_dist.png"
        )
        logger.info("✅ Graphiques générés")
        
        # 5. Générer rapport HTML
        logger.info("\n📄 Génération du rapport HTML...")
        html_file = create_html_report(stats_df, genome_name, output_dir)
        logger.info(f"✅ Rapport HTML : {html_file}")
        
        # 6. Récapitulatif final
        logger.info("\n" + "="*60)
        logger.info("✅ ANALYSE TERMINÉE AVEC SUCCÈS")
        logger.info("="*60)
        logger.info(f"\n📁 Fichiers générés dans : {output_dir}/")
        logger.info(f"  • {genome_name}_stats.csv")
        logger.info(f"  • {genome_name}_length_dist.png")
        logger.info(f"  • {genome_name}_gc_dist.png")
        logger.info(f"  • {genome_name}_report.html")
        logger.info(f"  • {genome_name}_analysis.log")
        logger.info("")
        logger.info(f"🌐 Ouvrez le rapport HTML dans votre navigateur :")
        logger.info(f"   {html_file.absolute()}")
        logger.info("")
        
    except FileNotFoundError as e:
        logger.error(f"❌ ERREUR : {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"❌ ERREUR INATTENDUE : {e}")
        logger.exception("Détails de l'erreur :")
        sys.exit(1)


if __name__ == "__main__":
    main()
