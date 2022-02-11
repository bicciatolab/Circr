# _Circr_

`Circr` is a computational tool for the prediction of miRNA:circRNA associations. It combines publicly available algorithms for de novo prediction of miRNA binding sites on target sequences with experimentally validated miRNA-AGO and miRNA:RNA sites in different organisms. `Circr` can be used with either the provided support files or with custom ones, allowing the analysis of novel and not previously annotated circRNAs in virtually any species of interest. `Circr` provides as output an annotated table that allows the user an easy selection of interesting miRNA:circRNA sites for validation and functional investigations.

Key features of Circr include:

1. Multifaceted comparison of interaction based on different prediction software
2. Self defined internal database for assessment of validated interactions and Argonaute peak occupancy
3. Fast calculation and multi thread compatibility
4. Easy to investigate excel structured output table

#### Contacts:

silvio.bicciato@unimore.it, mattia.forcato@unimore.it


## Table of contents


* [System requirements](https://github.com/bicciatolab/Circr#System-requirements)
* [Installation](https://github.com/bicciatolab/Circr#Installation)
* [Quick start guide](https://github.com/bicciatolab/Circr#Quick-start-guide)


## System requirements

* Python version: >= 3.6
* [miRanda](http://www.microrna.org/microrna/getDownloads.do)
* [RNAhybrid](https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid?id=rnahybrid_view_download)
* [TargetScan](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi)
* [Samtools](http://www.htslib.org/)
* [BEDtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
* Dependencies: _pandas_, _numpy_, _argparse_, _collections_, _pybedtools_, _functools_, _multiprocessing_, _itertools_, _logging_, _re_, _operator_

## Installation

`Circr` is structured as a standalone python script. Once the required python modules, `miRanda` and `RNAhybrid`  have been installed, `Circr` can be run in the suitable environment.

#### Installing `miRanda` and `RNAhybrid`
For the installation of `miRanda` and `RNAhybrid`, please refer to the instructions provided on their website. In case the sites are unavailable, it is still possible to install both tools through [Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) from command line. For `miRanda`, run either one of the following:

```bash
  conda install -c bioconda miranda
  conda install -c bioconda/label/cf201901 miranda
```
while for `RNAhybrid` use:

```bash
  conda install -c genomedk rnahybrid
```
After the installation has been completed for both tools, the executable file must be copied in either `/usr/bin/` or `/usr/local/bin/` folder. For `miRanda` it can be found in the `/path/to/anaconda3/pkgs/miranda<version_number>/bin/` folder, while for `RNAhybrid` it is located in the `/path/to/anaconda3/pkgs/RNAhybrid<version_number>/src/`

#### Installing Pyhton libraries
All required libraries can be installed through `pip` like ```pip3 install pybedtools```

#### Support files
`Circr` requires a set of annotation files to perform the various step of the analysis, which are listed below:

* Reference genome sequence in FASTA format
* Gene annotation in GTF format
* rRNA coordinates in BED format
* miRNA sequences in FASTA format
* validated miRNA:RNA interactions in BED format
* Argonaute peak file in BED format

For simplicity, we provide all the necessary files in a dedicated [drive folder](https://drive.google.com/drive/folders/1zJVyzEFAMtvZTTueWRocxXs63jUxsl-U?usp=sharing) that can be downloaded (at least 8Gb of space is required) within the `Circr` folder downloaded from `GitHub`. This folder contains the support files for 4 different organisms (_human_, _mouse_, _worm_ and _fruitfly_) in 2 different genome versions for each. If only one `organism` or `genome version` is necessary, it is possible to retrieve only the desired one as long as the expected folder tree structure (see below) is maintained.

```
  Circr
   |_ Circr.py
   |_ tools
   |   |_ targetscan_70.pl
   |_ support_files
       |_miRNA
       |_<organism_of_interest>
          |_<assembly_of_interest>

```

After downloading the support_files folder, it is necessary to navigate to the `assembly_of_interest` folder, unzip the genome `FASTA` file and index it

```bash
  cd /path/to/Circr/support_files/organism_of_interest/assembly_of_interest
  gunzip assembly_of_interest.fa.gz
  samtools faidx assembly_of_interest.fa
```

In case a custom set of files is available, it is possible to provide them directly to `Circr` (please, refer to the [Quick start guide](https://github.com/bicciatolab/CircR#Quick-start-guide) section for further information).

## Quick start guide
Below is provided a list of available options together with the default values for each parameter.

```
python3 Circr.py -h


usage: python3 Circr.py [-h] [-i INPUT] [-c] [-s ORGANISM] [-v GENOME_VERSION] [--gtf GTF] [--genome GENOME]
                       [--rRNA RRNA] [--miRNA MIRNA] [--AGO AGO] [--validated_interactions VALIDATED_INTERACTIONS]
                       [--threads THREADS] [-o OUTPUT]

Circr help section

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        List of Circular RNA or their exon coordinates.
  -c, --coord           Defines if the input file contains exon coordinates rather than Circular RNA coordinates
  -s ORGANISM, --organism ORGANISM
                        Defines the selected organism for analysis. Available organisms are: human, mouse, worm,
                        fruitfly. Default is human.
  -v GENOME_VERSION, --genome_version GENOME_VERSION
                        Defines the genome version of the selected organism. Versions available are: human (hg19,
                        hg38), mouse (mm9, mm10), worm (ce10, ce11), fruitfly (dm3, dm6). Default for human is hg38
  --gtf GTF             Alternative location for GTF file for the analysis. Default uses data for selected organism.
  --genome GENOME       Alternative location for genome file for the analysis. Default uses data for selected
                        organism.
  --rRNA RRNA           Alternative location for ribosomal RNA file for the analysis. Default uses data for selected
                        organism.
  --miRNA MIRNA         Alternative location for miRNA file for the analysis. Default uses data for selected organism.
  --AGO AGO             Alternative AGO peaks file. Default uses data for selected organism.
  --validated_interactions VALIDATED_INTERACTIONS
                        Alternative validated interactions file. Default uses data for selected organism.
  --threads THREADS     Set the number of threads for multiprocess. Default is 8 cores
  -o OUTPUT, --output OUTPUT
                        Defines output file name. Default is Circr_Analysis.txt
```

### Test data
In the `docs` folder we provide a set of circRNAs predicted to be expressed in the developing mouse brain (Dori et al., 2019). Data were obtained from RNAseR treated RNA of the three main cell population of the lateral cortex (proliferating and differentiating progenitors and newborn neurons). From this cohort of sequences (available as `Supplemental Data 1` of Dori et al., 2019), we selected 100 circRNAs of different lengths overlapping genes in the sense and antisense strand, and intergenic ones to cover all possible genomic features of circRNAs. This short dataset can be used as input to test that everything has been setup correctly.


### Running `Circr` with the provided annotation files

The `INPUT` file is a list of circRNA coordinates in BED format, an example provided below:

```
chr1	137428925	137433876	CiCo_mm9_circ_000076	.	-
chr1	154702978	154706632	CiCo_mm9_circ_000105	.	-
chr1	159104273	159109383	CiCo_mm9_circ_000118	.	-
chr1	159216861	159229395	CiCo_mm9_circ_000119	.	-
chr17	39980067	39980222	CiCo_mm9_circ_002584	.	-

```

If users have a FASTA file with the circRNA sequence, they must first align the sequence to the genome using tools such as BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi) or BLAT (https://genome.ucsc.edu/cgi-bin/hgBlat) setting the appropriate organism and genome version. Then the resulting genomic coordinates must be formatted into a tab delimited file compatible with Circr. 

`Circr` can be run with the following command

```bash
  python3 /path/to/Circr.py -i circular_mixed.bed -s mouse -v mm9 -o example_circr.csv
```

`-s` and `-v` parameters instruct `Circr` to use the mm9 genome version for mouse and to save the data in the `example_circ.csv` file. By default, the analysis will run with 8 threads but this number can be changed through the `--threads` option.

The first steps of the analysis will evaluate whether each circRNA overlaps a gene annotated on the same strand. In this case, the circRNA will be "split" into its exon/intron coordinates and only the exon will be kept to retrieve the circRNA FASTA sequence. Intergenic circRNAs and antisense to genes will be treated as a single exon (i.e. no splitting into features). Since circRNAs can undergo alternative splicing, it is possible that not all exons are included or the exact features included are not known. In this case, it is possible to run Circr with the `-c` or `--coord` option with the same file as shown before or supply an `INPUT` that includes the coordinates of all the expressed features (see an example below)

```
chr1	137428925	137429291	CiCo_mm9_circ_000076	.	-
chr1	137433773	137433876	CiCo_mm9_circ_000076	.	-
chr1	154702978	154703132	CiCo_mm9_circ_000105	.	-
chr1	154703578	154703741	CiCo_mm9_circ_000105	.	-
chr1	154706008	154706144	CiCo_mm9_circ_000105	.	-
chr1	154706481	154706632	CiCo_mm9_circ_000105	.	-
chr1	159216861	159216985	CiCo_mm9_circ_000119	.	-
chr1	159229267	159229395	CiCo_mm9_circ_000119	.	-
chr17	39980067	39980222	CiCo_mm9_circ_002584	.	-

```

The `--coord` option will instruct `Circr` to skip the initial exon/intron splitting and will directly use the coordinates present in the `INPUT` BED file to retrieve the `FASTA` sequence. An example of how to run `Circr` in the coordinates mode is

```bash
  python3 /path/to/Circr.py -i circular_mixed.bed --coord -s mouse -v mm9 -o example_coord_circr.csv
```

### Running `Circr` with custom annotation files
There are two possible options to run `Circr` with custom files:

1. Specify organism and genome version and then provide one or more custom files
2. Leave the default organism and genome version, but then supply ALL annotation files

#### Providing one or more files
This case applies, for example, if users are performing the analysis on tissue-specific circRNAs and have a set of tissue-specific annotation file for one of the 4 default organisms (human, mouse, worm, fruitfly). In this case, `Circr` can be run as follows:

```bash
  python3 /path/to/Circr.py -i circular_mixed.bed -s mouse -v mm9 \
    --AGO /path/to/my_custom_AGO_peaks.bed --validated_interactions /path/to/my_custom_interactions.bed \
    -o example_custom_anno_circr.csv
```

#### Providing ALL files
This second case applies, for example, if users are analysing circRNAs from an organism other than the 4 provided. In this case, we suggest to omit the `organism` and `genome version` arguments and then manually supply all the necessary files. `Circr` can be run as:

```bash
  python3 /path/to/Circr.py -i circular_mixed.bed \
    --gtf /path/to/my_gene_anno.gtf --genome /path/to/my_genome_ref.fa \
    --rRNA /path/to/my_rRNA_coord.bed --miRNA /path/to/my_miRNA_seq.fa \
    --AGO /path/to/my_custom_AGO_peaks.bed --validated_interactions /path/to/my_custom_interactions.bed \
    -o example_custom_anno_circr.csv
```

### Output description

Indipendently of how `Circr` is run, it will return a comma separated file consisting of 11 fields that are described below:

<p align="center">
<img src="https://github.com/bicciatolab/Circr/blob/main/docs/Table1.png" width="80%" alt="Circr-output">
</p>


and a few example lines are provided

```
Chrom,Start,End,miRNA Name,Circ Name,Strand,Seed Category,ID,Software Matched,Validated,AGO
chr1,10192694,10192715,mmu-let-7g-5p,CiCo_mm9_circ_000003,-,8mer,INT_M_0,1,No,No
chr1,10202507,10202528,mmu-let-7g-3p,CiCo_mm9_circ_000003,-,7mer-m8,INT_M_1,1,No,No
chr1,10201577,10201598,mmu-let-7i-5p,CiCo_mm9_circ_000003,-,No Seed,INT_M_2,1,No,No
chr1,10192680,10192702,mmu-let-7i-3p,CiCo_mm9_circ_000003,-,7mer-m8,INT_M_3,1,No,No
chr1,10201743,10201765,mmu-miR-1a-1-5p,CiCo_mm9_circ_000003,-,7mer-m8,INT_M_4,1,No,No
```
