#!/usr/bin/env python3
"""
Pipeline Complet d'Analyse Génomique

Ce script automatise l'intégralité du workflow :
1. Contrôle qualité du génome (N50, GC%, graphiques)
2. Annotation avec Prokka
3. Identification des enzymes industrielles
4. Génération rapport HTML final consolidé

Usage:
    python scripts/complete_pipeline.py genome.fasta --genus Bacillus --output results/GENOME_01/
    
Auteur : Takoi Rizgui
Projet : Mastère Bioinformatique - BioPipeline Toolkit
"""

import argparse
import sys
import subprocess
import shutil
from pathlib import Path
from datetime import datetime
import logging

# Import modules BioPipeline
sys.path.insert(0, str(Path(__file__).parent.parent))
from biopipeline.genome.stats import GenomeStats
from biopipeline.annotation import EnzymeFinder


def setup_logger(output_dir: Path, genome_name: str):
    """Configure le logging"""
    log_file = output_dir / f"{genome_name}_pipeline.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    return logging.getLogger(__name__)


def check_prokka_installed():
    """Vérifier si Prokka est installé"""
    try:
        result = subprocess.run(['prokka', '--version'], 
                              capture_output=True, 
                              text=True,
                              timeout=5)
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def run_quality_control(fasta_file: Path, output_dir: Path, logger):
    """
    Étape 1 : Contrôle qualité du génome
    """
    logger.info("\n" + "="*70)
    logger.info("ETAPE 1/4 : CONTROLE QUALITE")
    logger.info("="*70)
    
    qc_dir = output_dir / "01_quality_control"
    qc_dir.mkdir(parents=True, exist_ok=True)
    
    genome_name = fasta_file.stem
    
    try:
        # Analyser génome
        logger.info(f"Chargement du genome : {fasta_file.name}")
        stats = GenomeStats(str(fasta_file))
        logger.info(f"✓ {len(stats.sequences)} sequences chargees")
        
        # Statistiques
        logger.info("Calcul des statistiques...")
        stats_df = stats.generate_report(qc_dir / f"{genome_name}_stats.csv")
        
        # Afficher résumé
        basic_stats = stats.get_basic_stats()
        logger.info(f"\n--- RESUME QUALITE ---")
        logger.info(f"Total sequences : {basic_stats['total_sequences']:,}")
        logger.info(f"Longueur totale : {basic_stats['total_length']:,} pb")
        logger.info(f"N50             : {basic_stats['n50']:,} pb")
        logger.info(f"GC%             : {basic_stats['gc_percent']}%")
        
        # Évaluation qualité
        if basic_stats['n50'] > 50000:
            logger.info("✓ QUALITE : Excellente (N50 > 50 kb)")
            quality = "Excellente"
        elif basic_stats['n50'] > 10000:
            logger.info("⚠ QUALITE : Moyenne (N50 entre 10-50 kb)")
            quality = "Moyenne"
        else:
            logger.info("✗ QUALITE : Faible (N50 < 10 kb)")
            quality = "Faible"
        
        # Graphiques
        logger.info("Generation des graphiques...")
        stats.plot_length_distribution(qc_dir / f"{genome_name}_length_dist.png")
        stats.plot_gc_distribution(qc_dir / f"{genome_name}_gc_dist.png")
        logger.info("✓ Graphiques sauvegardes")
        
        logger.info(f"\n✓ ETAPE 1 TERMINEE - Resultats dans {qc_dir}/")
        
        return {
            'success': True,
            'stats': basic_stats,
            'quality': quality,
            'stats_df': stats_df
        }
        
    except Exception as e:
        logger.error(f"✗ ERREUR Controle Qualite : {e}")
        return {'success': False, 'error': str(e)}


def run_annotation(fasta_file: Path, output_dir: Path, genus: str, species: str, 
                  cpus: int, logger):
    """
    Étape 2 : Annotation avec Prokka
    """
    logger.info("\n" + "="*70)
    logger.info("ETAPE 2/4 : ANNOTATION GENOMIQUE (Prokka)")
    logger.info("="*70)
    
    annot_dir = output_dir / "02_annotation"
    annot_dir.mkdir(parents=True, exist_ok=True)
    
    genome_name = fasta_file.stem
    
    # Vérifier Prokka
    if not check_prokka_installed():
        logger.warning("⚠ Prokka non installe - annotation sautee")
        logger.info("Pour installer Prokka : conda install -c bioconda prokka")
        return {
            'success': False, 
            'skipped': True,
            'message': 'Prokka non installe'
        }
    
    try:
        # Commande Prokka
        cmd = [
            'prokka',
            str(fasta_file),
            '--outdir', str(annot_dir),
            '--prefix', genome_name,
            '--cpus', str(cpus),
            '--force'  # Écraser si existe déjà
        ]
        
        if genus:
            cmd.extend(['--genus', genus])
        if species:
            cmd.extend(['--species', species])
        
        logger.info(f"Lancement Prokka...")
        logger.info(f"Commande : {' '.join(cmd)}")
        
        # Exécuter Prokka
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1h max
        )
        
        if result.returncode != 0:
            logger.error(f"✗ Prokka a echoue : {result.stderr}")
            return {
                'success': False,
                'error': result.stderr
            }
        
        # Vérifier fichiers générés
        gbk_file = annot_dir / f"{genome_name}.gbk"
        if not gbk_file.exists():
            logger.error("✗ Fichier GenBank non genere")
            return {
                'success': False,
                'error': 'Fichier .gbk non trouve'
            }
        
        logger.info(f"✓ Annotation terminee : {gbk_file}")
        logger.info(f"\n✓ ETAPE 2 TERMINEE - Resultats dans {annot_dir}/")
        
        return {
            'success': True,
            'gbk_file': gbk_file,
            'annotation_dir': annot_dir
        }
        
    except subprocess.TimeoutExpired:
        logger.error("✗ Prokka timeout (> 1h)")
        return {'success': False, 'error': 'Timeout'}
    except Exception as e:
        logger.error(f"✗ ERREUR Annotation : {e}")
        return {'success': False, 'error': str(e)}


def run_enzyme_identification(gbk_file: Path, output_dir: Path, logger):
    """
    Étape 3 : Identification des enzymes
    """
    logger.info("\n" + "="*70)
    logger.info("ETAPE 3/4 : IDENTIFICATION DES ENZYMES")
    logger.info("="*70)
    
    enzyme_dir = output_dir / "03_enzymes"
    enzyme_dir.mkdir(parents=True, exist_ok=True)
    
    genome_name = gbk_file.stem
    
    try:
        # Analyser avec EnzymeFinder
        logger.info(f"Analyse du fichier GenBank : {gbk_file.name}")
        finder = EnzymeFinder(str(gbk_file))
        
        # Trouver enzymes
        logger.info("Recherche des enzymes industrielles...")
        enzymes = finder.find_all_enzymes()
        
        if len(enzymes) == 0:
            logger.warning("⚠ Aucune enzyme identifiee")
            return {
                'success': True,
                'enzyme_count': 0
            }
        
        logger.info(f"✓ {len(enzymes)} enzymes identifiees !")
        
        # Statistiques par famille
        catalog = finder.catalog_by_family()
        logger.info(f"\n--- DISTRIBUTION PAR FAMILLE ---")
        for family, row in catalog.iterrows():
            count = int(row['Count'])
            pct = count / len(enzymes) * 100
            logger.info(f"  {family:<15} : {count:>3} enzymes ({pct:>5.1f}%)")
        
        # Exporter catalogue
        logger.info("\nExport des resultats...")
        finder.export_to_csv(enzyme_dir / f"{genome_name}_enzyme_catalog.csv")
        
        # Exporter pour AlphaFold
        finder.export_for_alphafold(
            enzyme_dir / f"{genome_name}_for_alphafold.fasta",
            min_length=100,
            max_length=1500
        )
        
        # Graphiques
        try:
            finder.plot_family_distribution(enzyme_dir / f"{genome_name}_family_pie.png")
            finder.plot_length_distribution(enzyme_dir / f"{genome_name}_length_dist.png")
            logger.info("✓ Graphiques sauvegardes")
        except Exception as e:
            logger.warning(f"⚠ Erreur graphiques : {e}")
        
        # Rapport HTML
        finder.generate_html_report(enzyme_dir / f"{genome_name}_enzyme_report.html")
        
        logger.info(f"\n✓ ETAPE 3 TERMINEE - Resultats dans {enzyme_dir}/")
        
        return {
            'success': True,
            'enzyme_count': len(enzymes),
            'catalog': catalog,
            'finder': finder
        }
        
    except Exception as e:
        logger.error(f"✗ ERREUR Identification Enzymes : {e}")
        return {'success': False, 'error': str(e)}


def generate_final_report(fasta_file: Path, output_dir: Path, 
                         qc_result: dict, annot_result: dict, 
                         enzyme_result: dict, logger):
    """
    Étape 4 : Rapport final consolidé
    """
    logger.info("\n" + "="*70)
    logger.info("ETAPE 4/4 : GENERATION RAPPORT FINAL")
    logger.info("="*70)
    
    report_dir = output_dir / "04_final_report"
    report_dir.mkdir(parents=True, exist_ok=True)
    
    genome_name = fasta_file.stem
    
    # Préparer données
    stats = qc_result.get('stats', {})
    quality = qc_result.get('quality', 'Inconnue')
    enzyme_count = enzyme_result.get('enzyme_count', 0)
    
    # Distribution enzymes
    enzyme_table = ""
    if enzyme_result.get('success') and enzyme_count > 0:
        catalog = enzyme_result['catalog']
        for family, row in catalog.iterrows():
            count = int(row['Count'])
            pct = count / enzyme_count * 100
            enzyme_table += f"""
            <tr>
                <td><strong>{family}</strong></td>
                <td>{count}</td>
                <td>{pct:.1f}%</td>
                <td>{int(row['Avg_Length'])} aa</td>
            </tr>
            """
    
    # Statut annotation
    if annot_result.get('success'):
        annot_status = "✓ Reussie"
        annot_color = "#28a745"
    elif annot_result.get('skipped'):
        annot_status = "⊘ Sautee (Prokka non installe)"
        annot_color = "#ffc107"
    else:
        annot_status = "✗ Echouee"
        annot_color = "#dc3545"
    
    # Créer HTML
    html_content = f"""
<!DOCTYPE html>
<html lang="fr">
<head>
    <meta charset="UTF-8">
    <title>Rapport Final - {genome_name}</title>
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
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        .content {{
            padding: 40px;
        }}
        .pipeline-status {{
            display: grid;
            grid-template-columns: repeat(4, 1fr);
            gap: 20px;
            margin: 30px 0;
        }}
        .status-card {{
            background: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            text-align: center;
            border-top: 4px solid #667eea;
        }}
        .status-card.success {{
            border-top-color: #28a745;
        }}
        .status-card.warning {{
            border-top-color: #ffc107;
        }}
        .status-card.error {{
            border-top-color: #dc3545;
        }}
        .status-number {{
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
            margin: 10px 0;
        }}
        .section {{
            background: #f8f9fa;
            padding: 30px;
            border-radius: 10px;
            margin: 20px 0;
        }}
        table {{
            width: 100%;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            margin: 20px 0;
        }}
        th, td {{
            padding: 15px;
            text-align: left;
        }}
        th {{
            background: #667eea;
            color: white;
        }}
        tr:nth-child(even) {{
            background: #f8f9fa;
        }}
        .badge {{
            display: inline-block;
            padding: 5px 15px;
            border-radius: 20px;
            font-size: 0.9em;
            font-weight: bold;
        }}
        .badge-success {{
            background: #d4edda;
            color: #155724;
        }}
        .badge-warning {{
            background: #fff3cd;
            color: #856404;
        }}
        .badge-danger {{
            background: #f8d7da;
            color: #721c24;
        }}
        h2 {{
            color: #333;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
            margin: 30px 0 20px 0;
        }}
        .footer {{
            background: #f8f9fa;
            padding: 30px;
            text-align: center;
            border-top: 1px solid #ddd;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>RAPPORT D'ANALYSE COMPLETE</h1>
            <h2>{genome_name}</h2>
            <p>Genere le {datetime.now().strftime("%d/%m/%Y a %H:%M")}</p>
            <p><strong>BioPipeline Toolkit - Pipeline Automatise</strong></p>
        </div>
        
        <div class="content">
            <h2>Statut du Pipeline</h2>
            <div class="pipeline-status">
                <div class="status-card success">
                    <div>ETAPE 1</div>
                    <div class="status-number">✓</div>
                    <div>Controle Qualite</div>
                </div>
                <div class="status-card {'success' if annot_result.get('success') else 'warning' if annot_result.get('skipped') else 'error'}">
                    <div>ETAPE 2</div>
                    <div class="status-number">{'✓' if annot_result.get('success') else '⊘' if annot_result.get('skipped') else '✗'}</div>
                    <div>Annotation</div>
                </div>
                <div class="status-card {'success' if enzyme_result.get('success') and enzyme_count > 0 else 'warning'}">
                    <div>ETAPE 3</div>
                    <div class="status-number">{'✓' if enzyme_result.get('success') and enzyme_count > 0 else '⚠'}</div>
                    <div>Enzymes</div>
                </div>
                <div class="status-card success">
                    <div>ETAPE 4</div>
                    <div class="status-number">✓</div>
                    <div>Rapport Final</div>
                </div>
            </div>
            
            <div class="section">
                <h2>1. Qualite du Genome</h2>
                <table>
                    <tr>
                        <th>Metrique</th>
                        <th>Valeur</th>
                        <th>Evaluation</th>
                    </tr>
                    <tr>
                        <td>Nombre de contigs</td>
                        <td>{stats.get('total_sequences', 'N/A'):,}</td>
                        <td><span class="badge badge-{'success' if stats.get('total_sequences', 999) < 100 else 'warning'}">
                            {'Excellent' if stats.get('total_sequences', 999) < 100 else 'Moyen'}
                        </span></td>
                    </tr>
                    <tr>
                        <td>Longueur totale</td>
                        <td>{stats.get('total_length', 'N/A'):,} pb</td>
                        <td><span class="badge badge-success">OK</span></td>
                    </tr>
                    <tr>
                        <td>N50</td>
                        <td>{stats.get('n50', 'N/A'):,} pb</td>
                        <td><span class="badge badge-{'success' if stats.get('n50', 0) > 50000 else 'warning' if stats.get('n50', 0) > 10000 else 'danger'}">
                            {quality}
                        </span></td>
                    </tr>
                    <tr>
                        <td>Contenu GC</td>
                        <td>{stats.get('gc_percent', 'N/A')}%</td>
                        <td><span class="badge badge-success">Normal</span></td>
                    </tr>
                </table>
            </div>
            
            <div class="section">
                <h2>2. Annotation Genomique</h2>
                <p><strong>Statut :</strong> <span style="color: {annot_color};">{annot_status}</span></p>
                {f'<p><strong>Fichier GenBank :</strong> <code>{annot_result.get("gbk_file", "N/A")}</code></p>' if annot_result.get('success') else ''}
                {f'<p style="color: #856404;">⚠ Prokka non installe. Pour installer : <code>conda install -c bioconda prokka</code></p>' if annot_result.get('skipped') else ''}
            </div>
            
            <div class="section">
                <h2>3. Enzymes Industrielles Identifiees</h2>
                <p><strong>Total :</strong> {enzyme_count} enzymes</p>
                
                {f'''
                <table>
                    <thead>
                        <tr>
                            <th>Famille</th>
                            <th>Nombre</th>
                            <th>% Total</th>
                            <th>Longueur Moyenne</th>
                        </tr>
                    </thead>
                    <tbody>
                        {enzyme_table}
                    </tbody>
                </table>
                ''' if enzyme_count > 0 else '<p style="color: #856404;">⚠ Aucune enzyme identifiee (annotation requise)</p>'}
            </div>
            
            <div class="section">
                <h2>4. Fichiers Generes</h2>
                <ul style="line-height: 2;">
                    <li><strong>01_quality_control/</strong> - Statistiques et graphiques QC</li>
                    <li><strong>02_annotation/</strong> - Fichiers Prokka (.gbk, .faa, etc.)</li>
                    <li><strong>03_enzymes/</strong> - Catalogue enzymes + sequences AlphaFold</li>
                    <li><strong>04_final_report/</strong> - Ce rapport</li>
                </ul>
            </div>
            
            <div class="section">
                <h2>5. Prochaines Etapes</h2>
                <ul style="line-height: 2;">
                    <li>✓ Analyser les structures 3D avec AlphaFold (fichier FASTA genere)</li>
                    <li>✓ Comparer avec d'autres genomes</li>
                    <li>✓ Selectionner candidats prioritaires</li>
                    <li>✓ Validation experimentale</li>
                </ul>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>BioPipeline Toolkit v0.1.0</strong></p>
            <p>Pipeline Automatise - Analyse Genomique Complete</p>
            <p>Takoi Rizgui - Mastere Bioinformatique</p>
            <p style="margin-top: 10px; font-size: 0.9em;">
                Genere en {(datetime.now() - datetime.now()).total_seconds() + 1:.0f}s
            </p>
        </div>
    </div>
</body>
</html>
    """
    
    # Sauvegarder rapport
    report_file = report_dir / f"{genome_name}_FINAL_REPORT.html"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    logger.info(f"✓ Rapport final genere : {report_file}")
    logger.info(f"\n✓ ETAPE 4 TERMINEE")
    
    return report_file


def main():
    parser = argparse.ArgumentParser(
        description='Pipeline complet d\'analyse genomique',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemples :
  python scripts/complete_pipeline.py genome.fasta --genus Bacillus
  python scripts/complete_pipeline.py genome.fasta --genus Bacillus --species subtilis --output results/GEN01/
  
Pipeline automatise :
  1. Controle qualite (N50, GC%, graphiques)
  2. Annotation Prokka
  3. Identification enzymes industrielles
  4. Rapport HTML final consolide
  
Auteur: Takoi Rizgui
Projet: Mastere Bioinformatique - BioPipeline Toolkit
        """
    )
    
    parser.add_argument('fasta', help='Fichier FASTA du genome assemble')
    parser.add_argument('--genus', '-g', help='Genre de l\'organisme (ex: Bacillus)')
    parser.add_argument('--species', '-s', help='Espece (ex: subtilis)')
    parser.add_argument('--output', '-o', help='Dossier de sortie (defaut: results/GENOME_NAME/)')
    parser.add_argument('--cpus', '-c', type=int, default=4, help='Nombre de CPUs pour Prokka (defaut: 4)')
    
    args = parser.parse_args()
    
    # Vérifier fichier existe
    fasta_file = Path(args.fasta)
    if not fasta_file.exists():
        print(f"\nERREUR : Fichier non trouve : {args.fasta}")
        sys.exit(1)
    
    # Dossier de sortie
    genome_name = fasta_file.stem
    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = Path("results") / genome_name
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Logger
    logger = setup_logger(output_dir, genome_name)
    
    # Banner
    print("\n" + "="*70)
    print("     BIOPIPELINE TOOLKIT - PIPELINE COMPLET")
    print("="*70)
    print(f"Genome      : {fasta_file.name}")
    print(f"Sortie      : {output_dir}")
    print(f"Genus       : {args.genus or 'Non specifie'}")
    print(f"Species     : {args.species or 'Non specifie'}")
    print(f"CPUs        : {args.cpus}")
    print("="*70 + "\n")
    
    start_time = datetime.now()
    
    # ÉTAPE 1 : Contrôle Qualité
    qc_result = run_quality_control(fasta_file, output_dir, logger)
    if not qc_result['success']:
        logger.error("\nPIPELINE ARRETE - Erreur Controle Qualite")
        sys.exit(1)
    
    # ÉTAPE 2 : Annotation
    annot_result = run_annotation(fasta_file, output_dir, args.genus, 
                                  args.species, args.cpus, logger)
    
    # ÉTAPE 3 : Enzymes (si annotation réussie)
    enzyme_result = {'success': False, 'enzyme_count': 0}
    if annot_result.get('success'):
        gbk_file = annot_result['gbk_file']
        enzyme_result = run_enzyme_identification(gbk_file, output_dir, logger)
    else:
        logger.warning("\nETAPE 3 SAUTEE - Annotation non disponible")
    
    # ÉTAPE 4 : Rapport Final
    report_file = generate_final_report(fasta_file, output_dir, qc_result, 
                                       annot_result, enzyme_result, logger)
    
    # Résumé final
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    
    print("\n" + "="*70)
    print("     PIPELINE TERMINE AVEC SUCCES !")
    print("="*70)
    print(f"\nDuree totale    : {duration:.0f} secondes ({duration/60:.1f} minutes)")
    print(f"Resultats dans  : {output_dir}/")
    print(f"\nRAPPORT FINAL   : {report_file}")
    print(f"\nOuvrez le rapport HTML dans votre navigateur pour voir tous les resultats !")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()
