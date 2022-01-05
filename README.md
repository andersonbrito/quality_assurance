# IMPACC Quality Assurance

This Snakemake pipeline performs a series of quality assurance steps. It includes five check points, where viral genomes are tagged 'PASS' or 'FAIL', based on: genome coverage; metadata integrity; mutation issues; pango lineage assignment; and molecular clock.

# Requirements

* python=3.6*
* pandas
* augur
* nextalign
* iqtree
* treetime
* pangolin


# Input files

* gisaid_hcov-19.fasta
* genomes.fasta
* metadata_nextstrain.tsv
* metadata_corelab.tsv
* impacc-virology-clin-sample.csv
* impacc-virology-clin-individ.csv
* batch_layout.csv
* reference_seq.fasta
* sequence.gb
* genemap.gff


# Installing and running the pipeline

All one needs to do to run this pipeline is to clone the GitHub directory, and activate the biogrid environment in the IMPACC server. All the packages listed above and their dependencies will be loaded.

```
git clone https://github.com/andersonbrito/quality_assurance.git
source /programs/biogrids.shrc
snakemake assurance
```


# Pipeline overview

![alt text](https://github.com/andersonbrito/quality_assurance/blob/master/images/workflow.svg "Steps of quality assurance")

__Figure 1. Workflow Overview__ 


## Filtering sequences by genome coverage

(*) `filter_coverage` is a checkpoint step. Here, genomes with more than 30% of sites with ambiguities (non-ATCG bases) will be flagged as low coverage. The quality assurance matrix with show the columns: `seq_coverage` (with the proportion of ATGC sites) and `seq_coverage_status` (with `PASS` or `FAIL`). In this step, long sequence headers are renamed to show only `sample_id`:

	`Sequence_Identifier:IM_15202.aid_14729.SARS-CoV-2|External_Sample_Identifier:0865-0071KD00-001_015-0015_1` >>> `0865-0071KD00-001`


## Generating list of representative genomes from SARS-CoV-2 lineages

`lineages_rep` generates a list of 500 genomes from distinct, randomly selected SARS-CoV-2 lineages, and two ancestral genomes to root the phylogeny ('Wuhan/Hu-1/2019', 'Wuhan/WH01/2019'). These genomes will be used in the `root2tip` analysis step. Genome strain names are sampled from the following repository:

	https://raw.githubusercontent.com/cov-lineages/pangoLEARN/master/pangoLEARN/data/lineages.downsample.csv


## Inspecting metadata files

(*) `inspect_metadata` is a checkpoint step. Here, the following files are inspected to detect samples with missing metadata:

- `metadata_corelab.tsv` (metadata file from core lab)
- `impacc-virology-clin-sample.csv` (sample metadata)
- `impacc-virology-clin-individ.csv` (patient metadata)
- `batch_layout.csv` (batch layout file)

Samples that fail at this step are flagged with 'FAIL' under a new column `metadata_status`, which is added in the quality assurance matrix. A `missing_metadata` column will also list metadata field that caused the failure. Likewise, columns `batch_issues` and `batch_status` are added, indicating inconsistencies in the batch layout file.


## Creating base dataset

`base_dataset` prepares the fasta and metadata files required to proceed with the next steps. For this step, a file `gisaid_hcov-19.fasta` containing all the genomes and a `metadata_nextstrain.tsv` file listing the metadata of those genomes must be provided. These files can be downloaded from gisaid.org, using the list of representative genomes generated under `lineages_rep` as query in the 'Search' webpage (as 'Input for the Augur pipeline').


## Combining metadata files

`combine_metadata` generates a combined metadata file, including metadata from publicly available genomes (mentioned above) and from IMPACC core labs.


## Creating a multifasta file

`multifasta` generates a fasta file containing all sequences listed in the combined metadata file.


## Generating MSA

`align` generates a multiple sequence alignment (MSA) of all sequences listed in the metadata above, using `nextalign`.


## Inspecting mutations

(*) `mutations` is a checkpoint step. Here, any genetic changes leading to 'nonstop', 'nonsense', 'frameshift' mutations will be flagged as 'FAIL', representing genomes that require further inspections due to possible sequencing or assembly issues.


## Assigning pangolin lineages

(*) `pangolin` is a checkpoint step. Here, core lab genomes will be assigned with pangolin lineages. Important: the software pangolin must be updated, so that the lineages recently designated can be assigned (`pangolin --update`). Visit the [pangolin website](https://cov-lineages.org/resources/pangolin/updating.html) for more information.

## Masking alignment sites

`mask` masks problematic sites. An updated list of sites to be masked can be provided. This [repository](https://github.com/W-L/ProblematicSites_SARS-CoV2) lists candidate sites to be masked.

## Inferring a Maximum Likelihood phylogenetic tree

`tree` generates a phylogenetic tree using `IQTree`.

## Performing root to tip analysis

(*) `root2tip` is a checkpoint step. Here, molecular clock outlier are identified, flagged with 'FAIL' under a new column `clock_status`, which is added in the quality assurance matrix.


## Generating final files

`assurance` is the last step of quality assurance. Here, the final quality assurance matrix is generated, and a fasta file with quality genomes produced by core labs is created, including only genomes that passed all the check points mentioned above.


# Execution

To run this pipeline, users can run each rule separately (eg. `snakemake filter_coverage`, `snakemake mutations`, etc), or can simply run all steps in one go, using `snakemake assurance`. As a result, a quality assurance matrix similar to the one below will be generated.

**Quality assurance matrix**

![alt text](https://github.com/andersonbrito/quality_assurance/blob/master/images/qamatrix.png "Quality assurance matrix")

