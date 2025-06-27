from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="riborez",
    version="1.0.0",
    author="Bjorn Shockey",
    author_email="bs128@rice.edu",
    description="A tool for ribosomal analysis and data management",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bbshockey/RiboRez_V1",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=[
        # No external dependencies needed - NCBI Datasets CLI is installed automatically
    ],
    entry_points={
        "console_scripts": [
            "riborez=riborez.cli:main",
        ],
    },
    keywords="bioinformatics, ribosomal, analysis, ncbi, genomes",
    project_urls={
        "Bug Reports": "https://github.com/bbshockey/RiboRez_V1/issues",
        "Source": "https://github.com/bbshockey/RiboRez_V1",
    },
) 