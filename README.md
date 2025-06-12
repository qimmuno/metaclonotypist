[![License](https://img.shields.io/pypi/l/metaclonotypist.svg)](https://github.com/qimmuno/metaclonotypist/blob/master/LICENSE)
[![Latest release](https://img.shields.io/pypi/v/metaclonotypist.svg)](https://pypi.python.org/pypi/metaclonotypist)

# Metaclonotypist

Metaclonotypist is a flexible, modular pipeline for the discovery of TCR metaclones. It is powered by the [pyrepseq](github.com/andim/pyrepseq) package for repertoire sequencing analysis.

## Features

- Automated identification of T cell metaclones from repertoire sequencing data
- HLA-association analysis with robust false discovery rate control
- A highly modular pipeline combining speed with accuracy, achieved respectively by combining the [symdel algorithm](https://arxiv.org/abs/2403.09010))  with advanced [TCR similarity metrics](https://doi.org/10.1016/j.cels.2024.12.006) and clustering algorithms.

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

### Outputs

This will create (if successful) the following outputs in the folder `examples/out`:
- a volcano plot of cluster-HLA associations
- a table of significant cluster-HLA associations
- a corresponding table reporting the TCRs associated with all identified metaclones
- a table of summary statistics and parameter values

The example data is small in size so the analysis should run in <10s. The analysis is based on a dataset (in `examples/data`) of the 30 top-most expanded clones at the site of a tuberculin-skin test from 150 individuals with associated HLA metadata.

### Run on custom data

To run on your own dataset:

```bash
metaclonotypist --tcrpath path/to/tcr.csv --hlapath path/to/hla.csv --output-dir my_results/
```

Refer to `examples/data/` for input file format.

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

