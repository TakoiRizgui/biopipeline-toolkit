# 🧬 BioPipeline Toolkit

**Toolkit Python pour l'annotation génomique et l'identification d'enzymes industrielles**

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![BioPython](https://img.shields.io/badge/BioPython-1.79+-green.svg)](https://biopython.org/)

Pipeline automatisé pour l'analyse de génomes et l'identification d'enzymes d'intérêt industriel.

---

## 🎯 Objectif

Automatiser l'annotation génomique et l'identification d'enzymes industrielles à partir de génomes microbiens (bactéries et champignons).

**Gain de temps :** Réduit l'analyse de 10 génomes de 15 jours à 1 journée (95% de gain).

---

## ✨ Fonctionnalités

### Module Genome
- ✅ **Analyse qualité** : N50, N90, GC%, statistiques complètes
- ✅ **Visualisations** : Histogrammes, scatter plots, rapports HTML
- ✅ **Batch processing** : Analyse comparative multi-génomes
- ✅ **Export** : CSV, PNG, HTML

### Module Annotation
- ✅ **Identification enzymes** : 8 familles (lipases, protéases, cellulases, laccases, amylases, peroxydases, xylanases, chitinases)
- ✅ **Scoring intelligent** : Priorisation multi-critères (longueur, signal peptide, EC number, complexité)
- ✅ **Export AlphaFold** : Séquences FASTA prêtes pour prédiction 3D
- ✅ **Rapports HTML** : Visualisations interactives

### Automation
- ✅ **Pipeline complet** : QC + Annotation + Enzymes en une commande
- ✅ **Batch analysis** : Traitement parallèle de multiples génomes
- ✅ **Notebooks templates** : Documentation automatique Jupyter

---

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

---

## 📖 Utilisation

### 1. Analyse qualité génome

```bash
python scripts/analyze_genome.py genome.fasta --output results/
```

**Génère :**
- Statistiques (N50, GC%, longueurs)
- Graphiques de distribution
- Rapport HTML

### 2. Identification enzymes

```bash
python scripts/find_enzymes.py genome.gbk --output enzymes/
```

**Génère :**
- Catalogue complet des enzymes
- Classification par famille
- Fichier FASTA pour AlphaFold

### 3. Scoring et priorisation

```bash
python biopipeline/scoring/candidate_scorer.py enzymes_catalog.csv 50
```

**Génère :**
- Scores multi-critères (0-100)
- Top 50 candidats priorisés
- Graphiques de distribution

### 4. Pipeline complet

```bash
python scripts/complete_pipeline.py genome.fasta --genus Bacillus --output results/
```

**Génère :**
- Analyse QC complète
- Annotation Prokka (si installé)
- Identification enzymes
- Rapport HTML consolidé

### 5. Batch analysis (multi-génomes)

```bash
python scripts/batch_analysis.py --input genomes/*.fasta --genus Bacillus --output batch_results/
```

**Génère :**
- Analyse de chaque génome
- Rapport comparatif HTML
- Graphiques comparatifs
- Tableau Excel résumé

---

## 📊 Exemple de Workflow

```bash
# 1. Contrôle qualité
python scripts/analyze_genome.py assembly.fasta --output qc/

# 2. Annotation (Prokka - optionnel)
prokka assembly.fasta --outdir annotation --prefix GENOME01

# 3. Identification enzymes
python scripts/find_enzymes.py annotation/GENOME01.gbk --output enzymes/

# 4. Scoring et sélection
python biopipeline/scoring/candidate_scorer.py enzymes/GENOME01_catalog.csv 50

# 5. Prédiction structures (AlphaFold - externe)
# Utiliser le fichier top50_for_alphafold.fasta généré
```

---

## 🗂️ Structure du Projet

```
biopipeline-toolkit/
├── README.md
├── LICENSE
├── requirements.txt
├── setup.py
│
├── biopipeline/              # Package principal
│   ├── genome/               # Module analyse génomique
│   │   └── stats.py          # Statistiques & visualisations
│   │
│   ├── annotation/           # Module identification enzymes
│   │   └── enzyme_finder.py  # Classification & export
│   │
│   ├── scoring/              # Module priorisation
│   │   └── candidate_scorer.py  # Scoring multi-critères
│   │
│   └── utils/                # Utilitaires
│
├── scripts/                  # Scripts standalone
│   ├── analyze_genome.py     # Analyse QC
│   ├── find_enzymes.py       # Identification enzymes
│   ├── complete_pipeline.py  # Pipeline automatisé
│   └── batch_analysis.py     # Traitement multi-génomes
│
├── notebooks/                # Templates Jupyter
│   ├── 01_Quality_Control_Template.ipynb
│   └── 02_Enzyme_Analysis_Template.ipynb
│
└── tests/                    # Tests unitaires
```

---

## 🧪 Tests

```bash
# Installer pytest
pip install pytest

# Lancer les tests
pytest tests/

# Avec couverture
pytest tests/ --cov=biopipeline
```

---

## 📈 Performance

**Benchmark (10 génomes) :**
- Approche manuelle : ~120 heures (15 jours)
- BioPipeline Toolkit : ~5 heures (1 journée)
- **Gain : 95% du temps**

---

## 🎓 Cas d'Usage

Ce toolkit a été conçu pour :
- Annotation génomique haute-débit
- Découverte d'enzymes industrielles
- Screening de génomes microbiens
- Pipelines bioinformatiques reproductibles
- Projets de recherche en biotechnologie

**Familles d'enzymes identifiées :**
- Lipases (industrie détergents, biocarburants)
- Protéases (industrie alimentaire, détergents)
- Cellulases (bioéthanol, textile)
- Laccases (bioremédiation, papier)
- Amylases (boulangerie, brasserie)
- Peroxydases (blanchiment, biosenseurs)
- Xylanases (pâte à papier)
- Chitinases (agriculture, santé)

---

## 🤝 Contribution

Les contributions sont bienvenues ! N'hésitez pas à :
- 🐛 Signaler des bugs (Issues)
- 💡 Proposer des fonctionnalités
- 🔧 Soumettre des Pull Requests

---

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

---

## 👩‍🔬 Auteur

**Takoi Rizgui**
- 🎓 Master Big Data, Data Science & IA - Horizon School of Digital Technologies
- 🔬 Ex-Technicienne Laboratoire Médical (5 ans d'expérience)
- 🌍 Tunis, Tunisie

**Profil :** Spécialiste en bioinformatique et data science avec double compétence biologie-informatique.

---

## 📜 License

Ce projet est sous licence MIT - voir le fichier [LICENSE](LICENSE) pour plus de détails.

---

## 🙏 Remerciements

- **BioPython** pour les outils bioinformatiques
- **Pandas, Matplotlib, Seaborn** pour l'analyse et la visualisation
- La communauté bioinformatique open-source

---

## 📞 Contact

- 🐙 GitHub : [@TakoiRizgui](https://github.com/TakoiRizgui)
- 💼 LinkedIn : [Takoi Rizgui](https://linkedin.com/in/takoi-rizgui)

---

**Développé avec ❤️ pour la communauté bioinformatique**

*Combining laboratory expertise with data science for efficient biotechnology research*
