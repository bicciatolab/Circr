# Circr

`Circr` is a computational tool for the prediction of circRNA-miRNA associations. It combines publicly available algorithms for de novo prediction of miRNA binding sites on target sequences with experimentally validated miRNA-AGO and miRNA-RNA sites in different organisms. `CircR` can be used with either the provided support files or with custom ones, allowing the analysis of novel and not previously annotated circRNAs in virtually any species of interest. `CircR` provides as output an annotated table that allows the user an easy selection of interesting circRNA-miRNA sites for validation and functional investigations.

Key features of `Circr` include:
1. Multifaceted comparison of interaction based on different high level validation software
2. Self defined internal database for assessment of validated interactions
3. Fast calculation and multi thread compatibility
4. Easy to investigate excel structured output table

 #### Contacts:
 silvio.bicciato.at.unimore.it,  mattia.forcato@unimore.it

 # Table of Contents

 - [System requirements](https://github.com/bicciatolab/CircR#System-requirements)
 - [Installation in Python](https://github.com/bicciatolab/CircR#installation-in-r)
 - [Tutorial on example data](https://raw.githack.com/bicciatolab/CircR/main/docs/CircR_tutorial.html)

 ## System requirements

 * Python version: >= 3.6
 * [miRanda](http://cbio.mskcc.org/miRNA2003/miranda.html)
 * [RNAhybrid](https://bibiserv.cebitec.uni-bielefeld.de/download/tools/rnahybrid.html)
 * [TargetScan](http://www.targetscan.org/vert_80/)
 * [Samtools](http://www.htslib.org/)
 * Dependencies: *pandas*, *numpy*, *argparse*, *collections*, *pybedtools*, *functools*, *multiprocessing*,
 *itertools*, *logging*, *re*, *operator*

 ## Installation in Python

 `Circr` is structured as a standalone python script. Once the required python modules have been installed,
 `Circr` can be run in the suitable environment.
