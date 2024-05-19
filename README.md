# In-Progress Overview
Suite of tools and (will be) a FAIR principled pipeline for phylogenomic analyses.


## PhyG -- Phylogenomic Pipeline
In-progress basic phylogenomic pipeline to assess sequence homology within individual multi-sequence alignments. PhyG incorporates MAFFT and Guidance2 to generate multi-sequence alignments and assess sequence homology, respectively.
By default, sequences identified as "non-homologous" (see seq-cut-off in [Guidance2 paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4489236/)) are removed and another round of alignment construction and assessment is performed (up to a maximum of 10 rounds, by default).

The resulting outputs will be sequence-filtered (not columns!) multi-sequence alignments with high-confidence homologous sequences. 

**Dependencies**
+ [Guidance2 2.02+](https://github.com/anzaika/guidance)
+ [MAFFT 7.5+](https://mafft.cbrc.jp/alignment/software/)
+ [TrimAl 1.4+](https://github.com/inab/trimal)
+ [BioPython 1.8+](https://biopython.org/wiki/Download)

Note: PhyG relies on these software being accessible within your environment's PATH for automatic detection.

### Inputs
- Directory of un-aligned gene family FASTA-formatted files (one FASTA-formatted file per gene family)
- Project name to store the results and intermediate files


```
phyg.py -m <MSA-directory> -p <project-name>
```

### List of all SUPPORTED options (some are "planned" or to be removed)

|    Command                |  Comment  |
|---------------------------|-----------|
| `-h`, `--help`  | Print the help message |
| `-m`, `--msa-dir <STRING>`  | Directory with un-aligned FASTA formatted files. |
| `-n`, `--porject-name <STRING>`  | Name for your project/outputs. |
| `-t`, `--threads <INTEGER>`  | Number of available threads to use. Default value is `1`. |
| `--max-iter <INTEGER>` | Maximum number of Guidance2 iterations to run. Default value is `10`.|
| `--seqscore`,`-ssc <FLOAT>` | Minimum sequence-score for Guidance2 homology assessment. Default value is `0.6`.|


## PhyG-Contam -- Contamination-Evaluation Pipeline -- IN PREP!
In-progress approach to assess potential contamination for individual taxa (genomes/transcriptomes -- transcriptomes only at the moment). 

This uses a set of 400 broadly conserved orthologous gene families, identifies homologous sequences from transcriptomes, and places them into existing single-gene phylogenies. Afterwards, the query taxon's "phylogenetic sisters" are assessed and presented in a summary table, which can be used to generate "contamination rules" for sequence removal if desired.

Ideally, transcriptomes are not contaminated, but if they are, predicted CDSs can be classified using [TIdeS v1.2+](https://github.com/xxmalcala/TIdeS). Future versions will include means to easily generate training data for TIdeS/

**Dependencies**
+ [DIAMOND 2.0.13+](https://github.com/bbuchfink/diamond)
+ [EPA-ng v3.8+](https://github.com/pierrebarbera/epa-ng)
+ [gappa v0.8+](https://github.com/lczech/gappa)
+ [MAFFT 7.5+](https://mafft.cbrc.jp/alignment/software/)
+ [BioPython 1.8+](https://biopython.org/wiki/Download)
+ [ete3 v3.1+](http://etetoolkit.org/download/)
+ [TIdeS v1.2+](https://github.com/xxmalcala/TIdeS) -- OPTIONAL!

### Inputs
- Transcriptome assembly (FASTA-formatted file; can be predicted CDSs or "raw" assembly)
- Taxon name
+ Can be just the genus name or genus and species ("genus_species"). Please ensure the name is found when searching through [NCBI's taxonomy page](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1)
- Directory of Reference MSAs -- link will be provided
- Directory of Reference Phylogenies -- link will be provided
- DIAMOND formatted reference database -- link will be provided

`phyg_contam.py <query-FASTA> <reference-database> <taxon-name> <Reference-MSA-Directory> <Reference-Phylogeny-Directory>`

