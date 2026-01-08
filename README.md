<div align="center">

<img src="https://raw.githubusercontent.com/qimmuno/metaclonotypist/main/metaclonotypist.svg" width=700>

[![Latest release](https://img.shields.io/pypi/v/metaclonotypist.svg)](https://pypi.python.org/pypi/metaclonotypist)
[![License](https://img.shields.io/pypi/l/metaclonotypist.svg)](https://github.com/qimmuno/metaclonotypist/blob/master/LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1101/2025.04.12.648537-pink)](https://doi.org/10.1101/2025.04.12.648537)
[![Zenodo](https://zenodo.org/badge/972076505.svg)](https://doi.org/10.5281/zenodo.17977728)

Metaclonotypist is a modular pipeline for the discovery of TCR metaclones powered by the [pyrepseq](github.com/andim/pyrepseq) package for repertoire sequencing analysis.

</div>

## Features

- Automated identification of T cell metaclones from repertoire sequencing data
- HLA-association analysis with robust false discovery rate control
- A highly modular pipeline combining speed with accuracy. This is achieved by combining the [Symdel algorithm](https://arxiv.org/abs/2403.09010) for fast edit distance sequence neighbor candidate identification with refinement by more complex similarity metrics. Metaclonotypist supports([TCRdist](https://doi.org/10.1038/nature22383) as its default similarity metric, as well as (currently experimental) [SCEPTR](https://doi.org/10.1016/j.cels.2024.12.006)) filtering. Metaclonotypist also supports different graph-based clustering algorithms, including the default Leiden clustering.

## Requirements

- Python 3.9 or later
- Install dependencies via pip (see `pyproject.toml`)
- Metaclonotypist has been tested to work with pyrepseq v1.5.1, pandas v2.2.1, numpy v1.26.4, scipy v1.12.0, statsmodels v0.14.1,  matplotlib v3.8.3, seaborn v0.13.2
- Linux or macOS recommended (Windows untested)

## Installation

```bash
pip install metaclonotypist
```

Note that installation might take a couple of minutes, if dependencies need to be installed.

## Usage

### Basic run on example data

To run the CLI with example data:

```bash
git clone https://github.com/qimmuno/metaclonotypist.git
cd metaclonotypist
pip install -e .
bash examples/run_cli_example.sh
```

The example data is small in size so the analysis should run in <10s. The analysis is based on a dataset (in `examples/data`) of the 30 top-most expanded clones at the site of a tuberculin-skin test from 150 individuals with associated HLA metadata.

### Outputs

This will create (if successful) the following outputs in the folder `examples/out`:
- `volcano_plot*.png`: a volcano plot of cluster-HLA associations
- `cluster_associations*.csv`: a table of significant cluster-HLA associations
- `clustering*.csv`: a corresponding table reporting the TCRs associated with all identified metaclones
- `stats*.csv`: a table of summary statistics and parameter values

`*` is a string reporting parameter settings used during the analysis.

<details>
<summary>Click to view full output documentation</summary>

#### Output file: `volcano_plot*.png`

This PNG file displays a volcano plot summarizing the results of the cluster-HLA association analysis. Each point on the plot represents a specific cluster-HLA allele combination, with the x-axis showing the log-transformed odds ratio and the y-axis showing the negative log10 p-value for the association. Cluster-HLA pairs with strong associations appear further from the origin. Colour indicates associations judged to be statistically significant following false discovery rate control. Infinite odds ratios (arising from perfect separation) are plotted at a large fixed value to ensure they are visible on the plot.

As a control the second panel shows the same analysis on data where the donor metadata was shuffled.

#### Output file: `cluster_associations.csv`

This CSV file contains the results of the cluster-HLA association analysis. The columns are:

- **cluster**: Identifier for the TCR cluster (metaclone).
- **hla**: HLA allele tested for association.
- **count_cluster**: Number of individuals within the cluster that have the specified HLA allele.
- **total_cluster**: Total number of individuals in the cluster.
- **count_other**: Number of individuals not part of the cluster with the specified HLA allele.
- **total_other**: Total number of individuals not part of the cluster.
- **pvalue**: P-value from the statistical test assessing the association between the cluster and the HLA allele.
- **odds_ratio**: Odds ratio quantifying the strength of the association between the cluster and the HLA allele.

#### Output file: `clustering*.csv`

This CSV file contains the mapping of TCRs to identified metaclones (clusters). The columns are:

- **index**: Unique identifier for each TCR sequence in the input data.
- **cluster**: Identifier for the metaclone (cluster) to which the TCR has been assigned.

Each row represents a TCR and the cluster it belongs to, allowing users to trace which TCRs are grouped together as metaclones.

#### Output file: `stats.csv`

This CSV file provides summary statistics and parameter values used in the analysis. The columns are:

- **parameter**: Name of the parameter or statistic.
- **value**: The corresponding value for the parameter or statistic.

Each row reports a specific parameter setting or summary metric, allowing users to track the configuration and results of their analysis.
</details>

### Run on custom data

To run on your own dataset:

```bash
metaclonotypist --tcrpath path/to/tcr.csv --hlapath path/to/hla.csv --output-dir my_results/
```

Refer to `examples/data/` for input file format.

<details>
<summary>Click to view full input documentation</summary>

#### Input file: `tcrdata.csv`

This CSV file contains TCR sequence data for each sample. The columns are:

- **TRBV**: The TCR beta variable gene segment (e.g., TRBV20-1).
- **TRBJ**: The TCR beta joining gene segment (e.g., TRBJ2-7).
- **CDR3B**: The amino acid sequence of the TCR beta chain CDR3 region.
- **Sample.ID**: Identifier for the sample or donor from which the TCR was derived.
- **clonal_count**: The number of times this TCR sequence was observed in the sample (clone count).

Each row represents a unique TCR sequence observed in a particular sample, along with its gene usage and abundance.

For alpha chain analysis please supply the argument `--chain alpha` to metaclonotypist, and replace `B` with `A` in the above, e.g. `TRBV` -> `TRAV`.

#### Input file: `metadata.csv`

This CSV file contains metadata for each sample (typically HLA genotypes). The columns are:

- **Sample.ID**: Identifier for the sample or donor (must match the `Sample.ID` in the TCR data).
- **HLA columns**: Each subsequent column represents a specific HLA allele for a given gene and copy (e.g., `DPA1.1`, `DPA1.2`, `B.1`, `B.2`, `DPB1.1`, `DPB1.2`, `A.1`, `A.2`, `C.1`, `C.2`, `DRB1.1`, `DRB1.2`, `DQAB1`, `DQAB2`, `DQAB3`, `DQAB4`). These columns record the HLA alleles present in each individual for the corresponding gene and copy.

Each row corresponds to a single donor, listing their HLA alleles for each locus. The HLA columns may vary depending on the typing resolution and available data, but should be consistent across all samples. Where an individual is homozygous at a particular locus, the same allele name can be repeated twice or one of the alleles can be left blank.

Note: This metadata file could also contain other donor characteristics that might be differentially associated with metaclone presence in the repertoire. We have only used the pipeline so far to test for HLA association, but it is very much possible using this same setup to test for other associations (e.g., disease status).
</details>

### Advanced usage

Run `metaclonotypist --help` for full usage instructions:


<details>
<summary>Click to view full help output</summary>

```text
usage: metaclonotypist [-h] --tcrpath TCRPATH --hlapath HLAPATH -o OUTPUT_DIR [--chain {alpha,beta}] [--tcrdistmethod {tcrdist,sceptr}] [--mincount MINCOUNT] [--maxtcrdist MAXTCRDIST]
                       [--clustering {leiden,multilevel}] [--hlatest {fisher,agresti-caffo}] [--mindonors MINDONORS] [--maxedits MAXEDITS] [--version]

options:
  -h, --help            show this help message and exit
  --tcrpath TCRPATH     Path to input TCR data (CSV file)
  --hlapath HLAPATH     Path to input HLA metadata (CSV file)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Path to the output directory
  --chain {alpha,beta}  chain to use (default: beta)
  --tcrdistmethod {tcrdist,sceptr}
                        TCR distance method (default: tcrdist)
  --mincount MINCOUNT   Minimum count for clones (default: None, no filtering)
  --maxtcrdist MAXTCRDIST
                        Maximum TCR distance (default: 15)
  --clustering {leiden,multilevel}
                        Clustering algorithm (default: leiden)
  --hlatest {fisher,agresti-caffo}
                        Statistical test method for HLA association (default: fisher)
  --mindonors MINDONORS
                        Minimum number of donors for HLA filtering (default: 4)
  --maxedits MAXEDITS   Maximum edits for TCR distance (default: 2)
  --version             Show the version of Metaclonotypist
```
</details> 

## Citing Metaclonotypist
Please cite our [preprint](https://doi.org/10.1101/2025.04.12.648537).

### BibTex
```bibtex
@article{turner_tst_2025,
	title = {Evolution of {T} cell responses in the tuberculin skin test reveals generalisable Mtb-reactive {T} cell metaclones},
	doi = {10.1101/2025.04.12.648537},
	journal = {biorXiv preprint},
	author = {Turner, Carolin T and Tiffeau-Mayer, Andreas and Rosenheim, Joshua and Chandran, Aneesh and Saxena, Rishika and Zhang, Ping and Jiang, Jana and Berkeley, Michelle and Pang, Flora and Uddin, Imran and Nageswaran, Gayathri and Byrne, Suzanne and Karthikeyan, Akshay and Smidt, Werner and Ogongo, Paul and Byng-Maddick, Rachel and Capocci, Santino and Lipman, Marc and Kunst, Heike and Lozewicz, Stefan and Rasmussen, Veron and Pollara, Gabriele and Knight, Julian C and Leslie, Alasdair and Chain, Benny M and Noursadeghi, Mahdad},
	year = {2025},
}
```

## License

Metaclonotypist is released under the MIT License.

## Contributing

Contributions, bug reports, and feature requests are welcome! Please open an issue or pull request on [GitHub](https://github.com/qimmuno/metaclonotypist).

