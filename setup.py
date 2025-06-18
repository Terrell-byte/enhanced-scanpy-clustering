from setuptools import setup, find_packages

setup(
    name="rna_clustering_algorithms",
    version="0.1.4",
    description="Modular clustering extension for scanpy",
    author=["Daniel Sutton", "Laurits Madsen", "Sebastian Svendsen", "Valdemar Fuglsang"],
    packages=find_packages(),
    install_requires=[
        "scanpy>=1.9.0",
        "numpy>=1.20.0",
        "anndata>=0.8.0",
        "scipy>=1.8.0",
        "scikit-learn>=1.7.0",
    ],
    python_requires=">=3.8",
) 