# PNMT SNP Analysis: Parkinson's Disease

Systematic computational analysis identifying deleterious single nucleotide polymorphisms in the PNMT gene associated with Parkinson's disease.

## Overview

This project analyzed 286 missense SNPs in the PNMT gene using multiple bioinformatics tools to identify high-risk variants that could serve as therapeutic targets or diagnostic markers for Parkinson's disease.

## Key Findings

- **3 high-risk deleterious SNPs identified:** L217Q, G79D, Y35D
- Molecular dynamics simulations revealed significant structural changes
- Published in peer-reviewed journal (Taylor & Francis)

## Methods

### Variant Prioritization
- **Tools used:** PROVEAN, PolyPhen-2, PredictSNP, DynaMut, CUPSAT
- **Consensus approach:** Variants predicted deleterious by ≥80% of tools

### Structural Analysis
- **Molecular dynamics:** GROMACS simulations (100ns)
- **Metrics analyzed:** RMSD, RMSF, radius of gyration, SASA
- **Structure modeling:** Homology modeling and validation

## Repository Structure

```
PNMT-SNP-Analysis/
├── data/
│   ├── raw/                    # Raw SNP data
│   └── processed/              # Filtered variants
├── scripts/
│   ├── 01_variant_filtering.py
│   ├── 02_pathogenicity_prediction.py
│   ├── 03_structural_analysis.py
│   └── 04_visualization.py
├── results/
│   ├── figures/
│   └── tables/
├── notebooks/
│   └── analysis.ipynb
└── README.md
```

## Technologies

- **Python:** Biopython, Pandas, NumPy, Matplotlib
- **Molecular Dynamics:** GROMACS
- **Structural Biology:** PyMOL, SWISS-MODEL
- **Statistics:** SciPy, statsmodels

## Publication

Das, P. et al. (2021). *In-silico Analysis of Deleterious Single Nucleotide Polymorphism in PNMT gene associated with Parkinson's Disease.* Taylor & Francis.

## Future Directions

- Wet-lab validation of top candidates
- Integration with transcriptomic data
- Extension to other neurodegenerative diseases

## Contact

Papiya Das - papiya2810das@gmail.com

## License

This project is licensed under the MIT License.
