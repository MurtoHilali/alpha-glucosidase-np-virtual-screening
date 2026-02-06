# Docking and virtual screening data for natural product alternatives to alpha-glucosidase inhibitor acarbose

This repository contains the **code, analysis notebooks, and workflows** used for the discovery of potential **natural product inhibitors of alpha-glucosidase**, with an emphasis on alternatives to acarbose.

The project performs large-scale virtual screening and docking against multiple alpha-glucosidase targets using **GNINA**, followed by statistical analysis, decoy validation, and redocking controls. Due to size constraints, all heavy docking outputs and large datasets are hosted externally on Zenodo.

---
![](/alpha-glucosidase-np-virtual-screening/aadd-workflow.drawio.png)

## Project Overview

Acarbose, an alpha-glucosidase inhibitor, is used to manage type 2 diabetes, but with limited effectivenes. Natural products provide a chemically diverse search space for identifying alternative inhibitors with improved profiles.

This project explores natural product candidates by:
- Docking large ligand libraries against multiple alpha-glucosidase structures
- Performing active–decoy virtual screening
- Validating docking protocols via redocking
- Ranking and filtering candidates based on CNN-based scoring


## Targets

Docking and screening were performed against the following alpha-glucosidase structures:

- [**2QMJ**](https://www.rcsb.org/structure/2QMJ)
- [**3BAJ**](https://www.rcsb.org/structure/3BAJ)
- [**5NN8**](https://www.rcsb.org/structure/5NN8)

For each target, cleaned PDB and PDBQT structures are provided.


## Repository Structure

```
├── 00_data
│ ├── high_similarity_candidates.csv (external; see Data Availability)
│ └── prep (data prep notebooks)
├── 01_proteins
│ ├── *_orig.pdb
│ ├── *_cleaned.pdb
│ └── *_cleaned.pdbqt
├── 02_ligands (not tracked; hosted externally)
├── 03_results (not tracked; hosted externally)
├── 04_controls
│ ├── 041_decoys
│ ├── 042_redock
│ ├── *_active_ids.txt
│ ├── *_active_smi.txt
│ ├── analysis (control analysis notebooks)
│ ├── create_decoys.py
│ └── decoy_docking.sh
├── 05_tools
│ ├── convert_to_pdbqt.sh
│ ├── create_ligands.py
│ └── gnina_log_parser.py
└── README.md

```

---

### Dependencies

This project uses several key external tools and datasets to support reproducible docking and screening workflows:

* **LUDe (LIDEB’s Useful Decoys)** — a decoy-generation toolkit that produces property-matched molecules for benchmarking and validation. This was used to generate decoy sets for control experiments and enrichment analysis. ([lideb-lude-v2.streamlit.app][1])

* **GNINA** — used to perform docking and rescoring of natural product ligands against alpha-glucosidase targets. ([GNINA][2])

* **COCONUT (Collection of Open Natural Products)** — natural products from COCONUT served as the primary screening library in this study, offering structural diversity well-suited for probing novel inhibitor space. ([COCONUT database][3])

These tools and resources are open access and form the foundation of the docking, decoy generation, and natural product screening pipelines used in this repository.

[1]: https://lideb-lude-v2.streamlit.app/ "Decoys: LIDEB's Useful Decoys"
[2]: https://github.com/gnina/gnina "GNINA"
[3]: https://coconut.naturalproducts.net/ "COCONUT database"


## Methods Summary

- **Docking engine:** GNINA  
- **Scoring:** CNN_VS (GNINA)  
- **Validation:** Active–decoy screening and redocking RMSD analysis  
- **Ligand sources:** Natural product libraries (e.g., COCONUT-derived subsets)  
- **Outputs:** Parquet-formatted score tables, docked poses for actives, and ranked candidate lists  

Detailed methods and analysis steps are documented in the notebooks under `00_data/prep/` and `04_controls/analysis/`.


## Data Availability

Due to GitHub file size limits, all large datasets and docking outputs are hosted on **Zenodo**.

**Zenodo DOI:**  
`10.5281/zenodo.18500700`

The Zenodo archive includes:
- Virtual screening score tables (`*_vs.parquet`)
- Active–decoy datasets (`*_decoys.parquet`)
- Redocking validation results (`*_redock.parquet`)
- Docked poses for active compounds
- High-similarity candidate lists

---

## Reproducibility

This repository is designed so that:
- All **analysis** can be reproduced using the Zenodo datasets
- Docking pipelines can be rerun using the provided scripts and configurations
- Large intermediate outputs are intentionally excluded from version control

Paths, scripts, and directory conventions assume execution on an HPC cluster environment.

---

## License

- **Code:** MIT License  
- **Data:** CC-BY 4.0 (see Zenodo record)

---

## Citation

If you use this code or data, please cite the Zenodo dataset:


```
Murto Hilali. Docking and virtual screening data for natural product
alternatives to alpha-glucosidase inhibitor acarbose. Zenodo. https://doi.org/10.5281/zenodo.18500700
```

---

## Notes

This repository is under active development and may be extended to additional targets, ligand libraries, or scoring methodologies in future releases.
