# 🧬 BioPipeline Toolkit

**Toolkit Python pour l'annotation génomique et la découverte d'enzymes industrielles**

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![BioPython](https://img.shields.io/badge/BioPython-1.79+-green.svg)](https://biopython.org/)

Développé dans le cadre d'un **Mastère de Recherche en Bioinformatique** (Projet Tuniso-Italien).

---

## 🎯 Objectif

Automatiser et standardiser les workflows d'annotation génomique pour accélérer la découverte d'enzymes industrielles à partir de génomes microbiens (bactéries et champignons).

## ✨ Fonctionnalités Actuelles

### Module Genome
- ✅ **Analyse qualité** : N50, N90, GC%, statistiques complètes
- ✅ **Visualisations** : Histogrammes, scatter plots, rapports HTML
- ✅ **Export** : CSV, PNG, HTML

### Scripts
- ✅ **Analyse rapide** : Tout-en-un pour analyse complète d'un génome

## 🚀 Installation

### Prérequis
- Python 3.8 ou supérieur
- pip

### Installation rapide

```bash
# Cloner le repository
git clone https://github.com/TakoiRizgui/biopipeline-toolkit.git
cd biopipeline-toolkit

# Installer en mode développement
pip install -e .
```

## 📖 Utilisation

### 1. Script d'analyse rapide (Recommandé)

```bash
python scripts/analyze_genome.py genome.fasta --output results/
```

**Génère automatiquement :**
- ✅ Statistiques CSV
- ✅ Graphiques PNG (distribution longueurs, GC%)
- ✅ Rapport HTML professionnel
- ✅ Fichier log

**Options :**
```bash
python scripts/analyze_genome.py genome.fasta \
    --output results/ \
    --min-length 1000  # Filtrer contigs < 1000 pb
```

### 2. Utilisation en Python

```python
from biopipeline.genome.stats import GenomeStats

# Analyser un génome
stats = GenomeStats("my_assembly.fasta")

# Obtenir statistiques
print(stats)  # Affiche résumé

# Générer rapport
report = stats.generate_report("genome_stats.csv")

# Créer graphiques
stats.plot_length_distribution("length_dist.png")
stats.plot_gc_distribution("gc_dist.png")
```

### 3. Analyse rapide en une ligne

```python
from biopipeline.genome.stats import quick_stats

quick_stats("genome.fasta", output_dir="./results/")
```

## 📊 Exemple de Sortie

### Rapport texte
```
Statistiques Génomiques - my_genome.fasta
==================================================
Nombre de séquences : 247
Longueur totale     : 4,832,451 pb
N50                 : 54,321 pb
N90                 : 12,345 pb
GC%                 : 52.3%
Plus long contig    : 245,678 pb
Plus court contig   : 1,234 pb
Longueur moyenne    : 19,563 pb
```

### Fichiers générés
```
results/
├── my_genome_stats.csv           # Statistiques
├── my_genome_length_dist.png     # Graphique longueurs
├── my_genome_gc_dist.png         # Graphique GC%
├── my_genome_report.html         # Rapport professionnel
└── my_genome_analysis.log        # Journal d'exécution
```

## 🗂️ Structure du Projet

```
biopipeline-toolkit/
├── README.md
├── LICENSE
├── requirements.txt
├── setup.py
│
├── biopipeline/              # Package principal
│   ├── __init__.py
│   ├── genome/               # Module analyse génomique
│   │   ├── __init__.py
│   │   └── stats.py          # ✅ Statistiques & visualisations
│   │
│   ├── annotation/           # 🚧 En développement
│   ├── analysis/             # 🚧 En développement
│   ├── structure/            # 🚧 En développement
│   ├── ml/                   # 🚧 En développement
│   └── utils/                # 🚧 En développement
│
├── scripts/                  # Scripts standalone
│   └── analyze_genome.py     # ✅ Analyse rapide
│
├── notebooks/                # Jupyter notebooks
├── data/                     # Données exemple
├── tests/                    # Tests unitaires
└── docs/                     # Documentation
```

## 🧪 Tests

```bash
# Installer pytest
pip install pytest

# Lancer les tests
pytest tests/

# Avec couverture
pytest tests/ --cov=biopipeline
```

## 🚧 Fonctionnalités à Venir

### Phase 2 (En développement)
- 🔄 **Module Annotation** : Wrapper Prokka, batch processing
- 🔍 **Module Enzyme Finder** : Identification CAZymes, lipases, protéases
- 📐 **Module Structure** : Batch AlphaFold, analyse sites actifs
- 🤖 **Module ML** : Prédiction thermostabilité, activité catalytique

### Phase 3 (Planifié)
- 📊 Dashboard interactif (Streamlit)
- 🔗 API REST
- 📦 Package PyPI
- 📚 Documentation complète (ReadTheDocs)

## 🎓 Contexte Académique

Ce projet a été développé dans le cadre d'un **Mastère de Recherche en Bioinformatique** :

**Titre du projet :** *Exploration du Potentiel Enzymatique de Microorganismes Tunisiens par Annotation Génomique et Prédiction Assistée par Intelligence Artificielle*

**Objectifs :**
- Annoter 5-10 génomes bactériens d'environnements extrêmes tunisiens
- Identifier 100-500 enzymes d'intérêt industriel (lipases, protéases, cellulases, etc.)
- Prédire structures 3D avec AlphaFold2
- Développer modèles ML pour prédire propriétés enzymatiques
- Publications : 2 articles dans revues internationales

**Collaboration :** Projet Tuniso-Italien - Biotechnologie Durable

## 📄 Citation

Si vous utilisez cet outil dans vos travaux, merci de citer :

```bibtex
@software{rizgui2026biopipeline,
  author = {Rizgui, Takoi},
  title = {BioPipeline Toolkit: Automated Genomic Annotation for Industrial Enzyme Discovery},
  year = {2026},
  url = {https://github.com/TakoiRizgui/biopipeline-toolkit}
}
```

## 👩‍🔬 Auteur

**Takoi Rizgui**
- 🎓 Mastère Recherche Bioinformatique - Projet Tuniso-Italien (2026)
- 🎓 Master Big Data, Data Science & IA - Horizon School of Digital Technologies (2024-2025)
- 🔬 Ex-Technicienne Laboratoire Médical (5 ans d'expérience)
- 🌍 Tunis, Tunisie

**Parcours :** Reconversion de la biologie médicale vers la data science et la bioinformatique, combinant expertise métier avec compétences techniques en IA/ML.

## 🤝 Contribution

Les contributions sont bienvenues ! N'hésitez pas à :
- 🐛 Signaler des bugs (Issues)
- 💡 Proposer des fonctionnalités
- 🔧 Soumettre des Pull Requests

## 📜 License

Ce projet est sous licence MIT - voir le fichier [LICENSE](LICENSE) pour plus de détails.

## 🙏 Remerciements

- **PubChem (NIH)** pour les données scientifiques
- **BioPython** pour les outils bioinformatiques
- **Projet Tuniso-Italien** pour le financement et l'encadrement
- **Horizon School of Digital Technologies** pour la formation en Data Science & IA
- La communauté bioinformatique open-source

---

## 📞 Contact

- 🐙 GitHub : [@TakoiRizgui](https://github.com/TakoiRizgui)
- 💼 LinkedIn : [Profil](https://linkedin.com/in/takoi-rizgui)
- 📧 Email : [Via GitHub]

---

**Fait avec ❤️ pour la communauté bioinformatique**

*Combining medical laboratory expertise with artificial intelligence for sustainable biotechnology*
