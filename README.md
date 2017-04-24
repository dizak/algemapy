# algemapy

Python pipeline for mapping alternative genome markers onto phylogenetic tree.

### Installation

1. Requirements.
  * conda users: install dependencies with

  ```
  conda env create --file /path/to/algemapy.yaml
  ```

  or if you do not have access (e.g. conda is installed system-wide)

  ```
  conda env create --file /path/to/algemapy.yaml -p /your/path/to/env/
  ```

  This environment includes:

  * mafft
  * iqtree

  * non-conda users: install dependencies listed in algemapy.yaml by any other means.

2. External scripts/programs.
  * [FLASH](https://ccb.jhu.edu/software/FLASH/)
  * [MAFFT](http://mafft.cbrc.jp/alignment/software/)
  * [iqtree](http://www.iqtree.org/)
  * [mgremap](http://bioputer.mimuw.edu.pl/gorecki/mgremap)
  * [headnode_notifier](https://github.com/dizak/headnode_notifier/releases)

3. How to install.
  1. Use python package manager to download and install dependencies.
  2. Add scripts to system path.

### Usage

###### Formatting the database:
  1. Reference tree:
    * Format: phylip.
    * Download source: [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy):
      1. Search for taxonomy group of interest using [Organism] field.
      2. Display results as common tree. This will redirect you to NCBI Taxonomy Browser.
      3. Check group of interest and click Choose.
      4. Save file as phylip.
    * Sanitize the reference tree using:
    ```
    agmdbf.py downloaded_tree.phy --output formatted_tree.phy --sanitize-ref-tree
    ```
  2. Reference genes:
    * Format: tsv.
    * Download source: [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/):
      1. Search for gene of interest within taxonomy group of interest using [Gene] and [Organism] fields.
      2. Send to file in Tabular (text) format.
      3. Download genes by ids and coordinates using
      ```
      agmdbf.py tabular_summary.txt --output reference_genes_sequences.fasta --download-from-tab-summary
      ```
