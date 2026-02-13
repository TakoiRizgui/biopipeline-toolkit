from setuptools import setup, find_packages 
 
setup( 
    name="biopipeline-toolkit", 
    version="0.1.0", 
    packages=find_packages(), 
    install_requires=[ 
        "biopython>=1.79", 
        "pandas>=1.3.0", 
        "numpy>=1.21.0", 
        "matplotlib>=3.4.0", 
        "seaborn>=0.11.0", 
        "scikit-learn>=0.24.0", 
    ], 
) 
