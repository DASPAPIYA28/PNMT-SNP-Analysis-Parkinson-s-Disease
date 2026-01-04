# PNMT SNP Analysis: Parkinson's Disease

**Published Research:** In-silico Analysis of Deleterious Single Nucleotide Polymorphism in PNMT gene associated with Parkinson's Disease  
**Publication:** Taylor & Francis (2022)  
**Authors:** Maansi, Papiya Das, Sineagha V, Sriraksha Prakash

## Overview

This project systematically analyzed 286 missense SNPs in the Phenylethanolamine N-methyltransferase (PNMT) gene to identify deleterious variants associated with Parkinson's disease. Using 17 computational prediction tools and molecular dynamics simulations, we identified three high-risk SNPs that significantly alter protein structure and function.

## Key Findings

- **3 high-risk deleterious SNPs identified:** L217Q (rs5641), G79D (rs371769078), Y35D (rs1187184066)
- Molecular dynamics simulations (100ns) revealed significant structural destabilization
- All three variants show reduced protein stability compared to wild-type
- Published peer-reviewed findings with experimental validation pathway

## Methodology

### Computational Pipeline

**1. Sequence-Based Prediction (10 tools)**
- PROVEAN, Mutation Assessor, PolyPhen-2
- Meta-SNP, PredictSNP1, SuSPect
- PMut, PhD-SNP, SNPs&GO, SNAP2

**2. Structure-Based Prediction (7 tools)**
- I-Mutant 2.0, DUET, MUPRO
- iStable, DynaMut, Align-GVGD, CUPSAT

**3. Conservation Analysis (2 tools)**
- Multiple Sequence Alignment (MSA): Clustal Omega, MAFFT, MUSCLE, T-Coffee
- ConSurf for evolutionary conservation scoring

**4. Molecular Dynamics Simulations**
- Software: GROMACS 2019.4
- Force field: GROMOS 54A7
- Simulation time: 100 nanoseconds
- Analysis: RMSD, RMSF, Rg, SASA, hydrogen bonding

### Analysis Parameters

```
Initial SNPs analyzed: 286 missense variants
First filter (sequence tools): 30 deleterious
Second filter (structure tools): 16 significant
Final filter (conservation + MDS): 3 high-risk
```

## Results Summary

### Structural Stability Analysis

| Mutation | RMSD (nm) | RMSF (nm) | Rg (nm) | SASA (nm²) | H-bonds | Overall Impact |
|----------|-----------|-----------|---------|------------|---------|----------------|
| Wild-type| 0.279     | 0.131     | 1.731   | 122.735    | 194     | Baseline       |
| L217Q    | 0.309     | 0.120     | 1.742   | 123.940    | 199     | Destabilizing  |
| G79D     | 0.300     | 0.134     | 1.747   | 124.796    | 195     | Most deleterious|
| Y35D     | 0.274     | 0.139     | 1.741   | 123.445    | 199     | Destabilizing  |

### Relative Deleterious Impact

G79D > Y35D > L217Q (based on cumulative MDS analysis)

## Repository Structure

```
PNMT-SNP-Analysis/
├── data/
│   ├── raw/
│   │   └── pnmt_snps_dbSNP.csv          # 286 missense SNPs from dbSNP
│   └── processed/
│       ├── filtered_30_snps.csv          # First filter results
│       ├── filtered_16_snps.csv          # Second filter results
│       └── final_3_high_risk.csv         # Final validated SNPs
├── scripts/
│   ├── 01_data_collection.py             # Retrieve SNPs from NCBI
│   ├── 02_sequence_analysis.py           # Run PROVEAN, PolyPhen-2, etc.
│   ├── 03_structure_analysis.py          # Run I-Mutant, DUET, etc.
│   ├── 04_conservation_analysis.py       # MSA and ConSurf analysis
│   ├── 05_prepare_gromacs.sh             # GROMACS setup script
│   ├── 06_molecular_dynamics.py          # MDS analysis
│   └── 07_visualization.py               # Generate plots
├── results/
│   ├── figures/
│   │   ├── rmsd_analysis.png
│   │   ├── rmsf_analysis.png
│   │   ├── radius_gyration.png
│   │   └── sasa_analysis.png
│   ├── tables/
│   │   ├── sequence_predictions.csv
│   │   ├── structure_predictions.csv
│   │   └── mds_results.csv
│   └── protein_structures/
│       ├── wild_type_1HNN.pdb
│       ├── L217Q_mutant.pdb
│       ├── G79D_mutant.pdb
│       └── Y35D_mutant.pdb
├── notebooks/
│   └── complete_analysis.ipynb           # Full analysis workflow
├── docs/
│   ├── methodology.md
│   └── supplementary_data.pdf
├── requirements.txt
└── README.md
```

## Technologies Used

**Programming & Analysis:**
- Python 3.8+ (Biopython, Pandas, NumPy, Matplotlib, Seaborn)
- R (for statistical analysis)
- Bash scripting (GROMACS workflows)

**Molecular Dynamics:**
- GROMACS 2019.4
- GROMOS 54A7 force field
- VMD for trajectory visualization

**Structural Biology:**
- PyMOL (structure visualization)
- UCSF Chimera 1.15
- SWISS-MODEL (homology modeling)

**Bioinformatics Tools:**
- Multiple sequence alignment tools
- Conservation analysis servers
- Protein stability prediction servers

## Installation

```bash
# Clone repository
git clone https://github.com/DASPAPIYA28/PNMT-SNP-Analysis.git
cd PNMT-SNP-Analysis

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install GROMACS (for MDS analysis)
# Follow instructions at: http://www.gromacs.org/
```

## Usage

### Complete Analysis Pipeline

```python
from scripts import data_collection, sequence_analysis, molecular_dynamics

# Step 1: Collect SNP data
snps = data_collection.fetch_pnmt_snps()

# Step 2: Run sequence predictions
sequence_results = sequence_analysis.analyze_all_tools(snps)

# Step 3: Filter high-risk variants
high_risk = sequence_analysis.filter_deleterious(sequence_results)

# Step 4: Run molecular dynamics
mds_results = molecular_dynamics.run_gromacs_analysis(high_risk)

# Step 5: Generate report
molecular_dynamics.generate_report(mds_results)
```

### Individual Analysis Examples

```python
# PROVEAN analysis
from scripts.sequence_analysis import run_provean
results = run_provean('data/raw/pnmt_snps.csv')

# RMSD calculation
from scripts.molecular_dynamics import calculate_rmsd
rmsd = calculate_rmsd('results/protein_structures/wild_type_1HNN.pdb',
                      'results/protein_structures/G79D_mutant.pdb')
```

## Key Results

### G79D Mutation (Most Deleterious)

**Structural Changes:**
- Mutant residue (Aspartic acid) is larger than wild-type (Glycine)
- Introduces negative charge in buried residue
- Loss of flexibility required for protein function
- Located in S-adenosyl-L-methionine binding region
- RMSD increased by 7.5% vs wild-type

**Functional Impact:**
- Disrupts local backbone conformation
- Affects ligand binding capacity
- Located near highly conserved position
- MDS shows significant destabilization

### Y35D Mutation

**Structural Changes:**
- Mutant residue smaller than wild-type (creates empty space)
- Introduces negative charge (wild-type is neutral)
- Loss of hydrophobic interactions
- Direct interaction with SKF ligand affected

**Functional Impact:**
- Ligand binding disruption likely
- RMSF peaks at residue 102 (SAH binding site)
- Affects multiple ligand binding sites

### L217Q Mutation

**Structural Changes:**
- Mutant residue (Glutamine) larger than wild-type (Leucine)
- May not fit in protein core
- Loss of hydrophobic interactions
- Affects α-helix stability

**Functional Impact:**
- Local stability alterations affect neighboring residue ligand contacts
- Increased average hydrogen bonds suggest compensatory mechanisms

## Clinical Significance

### Parkinson's Disease Connection

PNMT catalyzes 2N-methylation of β-carbolines, creating 2N-methylated β-carbolinium cations. These are structural analogs of MPP+ (neurotoxin causing Parkinsonism) and may contribute to:

1. Selective loss of PNMT-positive neurons in Parkinson's patients
2. Mitochondrial dysfunction in dopaminergic cells
3. Neuropathology in substantia nigra and locus ceruleus

### Protein Interaction Network

PNMT interacts with:
- Dopamine beta-hydroxylase (DBH) - Score: 0.990
- Catechol O-methyltransferase (COMT) - Score: 0.960
- Amine oxidase (MAOA/MAOB) - Score: 0.947/0.946
- Tyrosine 3-monooxygenase (TH) - Score: 0.890

## Publication

Das, P., Maansi, Sineagha V., & Sriraksha Prakash. (2022). *In-silico Analysis of Deleterious Single Nucleotide Polymorphism in PNMT gene associated with Parkinson's Disease.* Taylor & Francis.

**DOI:** [Available upon request]  
**PubMed ID:** [Available upon request]

## Future Directions

1. **Wet-lab validation** of top three candidate SNPs
2. **Integration with transcriptomic data** from Parkinson's patients
3. **Extension to other neurodegenerative diseases** (Alzheimer's, ALS)
4. **Functional assays** to measure enzymatic activity of mutant proteins
5. **Population genetics studies** to determine allele frequencies

## Citation

If you use this code or methodology, please cite:

```bibtex
@article{das2022pnmt,
  title={In-silico Analysis of Deleterious Single Nucleotide Polymorphism in PNMT gene associated with Parkinson's Disease},
  author={Das, Papiya and Maansi and Sineagha V. and Sriraksha Prakash},
  journal={Taylor \& Francis},
  year={2022}
}
```

## Contact

**Papiya Das**  
MSc Biotechnology, University of Leeds  
Email: papiya2810das@gmail.com  
LinkedIn: [linkedin.com/in/p--d](https://linkedin.com/in/p--d)  
GitHub: [@DASPAPIYA28](https://github.com/DASPAPIYA28)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Dr. T P Krishna Murthy (Project Guide)
- Department of Biotechnology, M.S. Ramaiah Institute of Technology
- Visvesvaraya Technological University, Belgaum
- GROMACS development team
- All bioinformatics tool developers whose servers we used

---

**Note:** This repository contains the computational analysis pipeline. Raw data and detailed results are available upon reasonable request for research purposes.
