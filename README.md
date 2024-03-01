# Reconstructing Signed Relations from Interaction Data: TUTORIAL

This repository contains the R code to reproduce some of the results in the study "Reconstructing signed relations from interaction data," published in Nature Scientific Reports. The φ-method introduced in this paper infers weighted signed relations—positive or negative social ties—from unsigned interaction data, such as proximity or communication logs.

## Paper Reference

For a comprehensive understanding of the methodology, findings, and implications of this study, please refer to the [full paper](https://www.nature.com/articles/s41598-023-47822-1).

## Datasets

The analysis utilizes different community datasets to validate the φ-method. These datasets include unsigned interaction data, which the method uses to infer the underlying signed social relations.
Instructions on where to find the data are in the respective R files.

## Code Overview

The repository contains two files that showcase how to apply the method to two different datasets:
  - highschool.R
  - karate.R

The directory 'functions' contains auxiliary functions.
The file `requirements.txt` lists all required R packages to run this tutorial.

## Getting Started

1. Clone the repository to your local machine.
2. Install the required R packages listed in `requirements.txt` from CRAN.
3. Install the R package adjHelpR from github (gi0na/adjHelpR).
4. Run the R notebooks.

## Requirements

- R environment (version >= 4.0.0)
- Listed R packages in `requirements.txt`.

## License

This project is licensed under the GNU AFFERO License - see the [LICENSE.md](LICENSE) file for details.

## Citation

If you use the code or data in your research, please cite our paper as follows:

Andres, G., Casiraghi, G., Vaccario, G., and Schweitzer, F. Reconstructing signed relations from interaction data. Sci Rep 13, 20689 (2023). https://doi.org/10.1038/s41598-023-47822-1
