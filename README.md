# algemapy

Python pipeline for mapping alternative genome markers onto phylogenetic tree.

### Installation

1. Requirements.
  * os
  * sys
  * glob
  * re
  * time
  * argparse
  * tqdm
  * pathos
  * Biopython
  * jinja2
  * pandas

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
    * Download source: [NCBI Taxonomy Browser](https://www.ncbi.nlm.nih.gov/taxonomy):
      1. Search for taxonomy group of interest using [Organism] field.
      2. Display results as common tree. This will redirect you to NCBI Taxonomy Browser.
      3. Check group of interest and click Choose.
      4. Save file as phylip.
