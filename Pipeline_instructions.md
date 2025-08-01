# Tick-borne pathogen surveillance bioinformatics pipeline                    

This repository contains a bioinformatics pipeline for the analysis of amplicon sequencing datasets of tick-associated microbes.  This was developed by the Stenglein lab at Colorado State University in collaboration with researchers in CDC's Division of Vector-Borne Diseases.

This pipeline is [described in this paper](https://pubmed.ncbi.nlm.nih.gov/37247570/).

## Contents

- [Running the pipeline](#Running-the-pipeline)
    - [Running from github](#Running-from-github)
    - [Running test datasets](#Running-test-datasets)
    - [Running a specific version:](#Running-a-specific-version)
    - [Making sure that you are running the latest version](#Making-sure-that-you-are-running-the-latest-version)
    - [Running in different computing environments](#Running-in-different-computing-environments)
    - [Running by cloning the pipeline's repository](#Running-by-cloning-the-pipeline-repository)
- [Inputs](#Inputs)
    - [Metadata file](#Metadata-file)
    - [Input fastq](#Input-fastq)
- [Outputs](#Outputs)
    - [Output directory](#Output-directory)
    - [Output file name prefixes ](#Output-file-name-prefixes)
- [Surveillance Report](#Surveillance-Report)
    - [Surveillance targets](#Surveillance-targets)
    - [Reference sequences](#Reference-sequences)
    - [Assignment of observed sequences to reference sequences](#Assignment-of-observed-sequences-to-reference-sequences)
    - [Mapping of reference sequence to surveillance targets ](#Mapping-of-reference-sequence-to-surveillance-targets)
    - [Calling positive hits](#Calling-positive-hits)
- [Primers and primer trimming](#Primers-and-primer-trimming)
- [QC reports](#QC-reports)
- [BLASTing of unassigned sequences](#BLASTing-of-unassigned-sequences)
- [Dependencies](#Dependencies)
    - [Nextflow](#Nextflow)
    - [Singularity ](#Singularity)
    - [Conda](#Conda)
    - [R libraries](#R-libraries)
    - [Python dependencies](#Python-dependencies)
- [Additional parameter information](#Additional-parameter-information)




## Running the pipeline

See the [dependencies section](#dependencies) below for information about the main dependencies required for running this pipeline(including nextflow and singularity).

### Running from github

The simplest way to run the pipeline is directly from github, like this:

```
nextflow run stenglein-lab/tick_surveillance -resume --metadata /path/to/metadata_xls --fastq_dir /path/to/fastq/directory -profile singularity
```

A copy of this pipeline is [also maintained in the CDCgov repository](https://github.com/CDCgov/tick_surveillance), so could be run from there, by replacing `stenglein-lab` with `CDCgov` in the above command.

You must specify two required inputs to the pipeline: the path to a metadata excel spreadsheet and the path to a directory containing input fastq.  See [this section](#inputs) for more information on required inputs.

### Running test datasets

The pipeline includes a handful of small datasets (<= 1000 read pairs) that are derived from real known positive (or known negative) datasets.  These are included in the [test directory](./test/) of the repository.  These datasets serve as positive and negative controls and allow you to test that the pipeline is working as expected.  To use these test datasets, run with the test profile, for instance:

```
nextflow run stenglein-lab/tick_surveillance -profile singularity,test
```

Or to run with conda:
```
nextflow run stenglein-lab/tick_surveillance -profile conda,test
```

The results of the test run will be placed in a `test/results` sub-directory.

### Running a specific version

To run a specific version of the pipeline, use the -r option, for example:

```
nextflow run stenglein-lab/tick_surveillance -profile singularity,test -r v1.0.7
```

### Making sure that you are running the latest version

Nextflow will cache the pipeline in `$HOME/.nextflow/assets/` and continue to use the cached version, even if the pipeline has newer versions on github.  To remove the locally cached version, which will force nextflow to download and cache the latest version, run:

```
nextflow drop stenglein-lab/tick_surveillance
nextflow run stenglein-lab/tick_surveillance -profile singularity,test
```

Alternatively, you can just delete the cached pipeline directory:
```
rm -rf ~/.nextflow/assets/stenglein-lab/tick_surveillance/
```

Running nextflow pipelines from github is [described in more detail here](https://www.nextflow.io/docs/latest/sharing.html).  


### Running in different computing environments

You will want to use a profile that matches your computing environment.  So, for instance, if running on an SGE HPC environment, you'd run something like:

```
nextflow run stenglein-lab/tick_surveillance -resume --metadata /path/to/metadata_xls --fastq_dir /path/to/fastq/directory -profile singularity,sge 
```

### Running by cloning the pipeline repository

It is also possible to download the pipeline code to a directory of your choosing.  This can be useful if, for instance, you want to modify or debug the code.  You can do this by cloning the repository (or a fork of the repository):

```
git clone https://github.com/stenglein-lab/tick_surveillance.git
cd tick_surveillance
nextflow run main.nf -resume --metadata /path/to/metadata_xls --fastq_dir /path/to/fastq/directory -profile singularity
```

## Inputs

The pipeline requires two inputs:

1. [A metadata file](#metadata-file).

2. [Sequence datasets in fastq format](#Input-fastq).  

### Metadata file

A metadata file in Microsoft Excel format must be input to the pipeline.  

1. The spreadsheet should contain a single tab
2. The single tab should contain a table that defines metadata for samples being analyzed.
3. The first row of the table should contain column names.
4. Subsequent rows should contain metadata with one row per dataset. 
5. One of the columns must be named Index (case sensitive), and the values in this column must match fastq file names.  In other words, the values in this column should match the sample names specified in the Illumina sample sheet.  

An example of a working metadata file [can be found here](./test/test_metadata.xlsx)

### Input fastq

Input sequence data is assumed to be Illumina paired-end data in separate read1 and read2 files.  Files can be compressed or not but it would be preferred to leave them as compressed files to save disk space.

The expected names of the fastq files are defined by the parameter `fastq_pattern`, whose default value is defined in nextflow.config as `*_R[12]_*.fastq*`.  This pattern can be overridden on the nextflow command line using the `--fastq_pattern` parameter.

The location of the fastq files is specified by the required `fastq_dir` parameter.  

It is expected that sample IDs are not repeated in the Illumina sample sheet.  

It is not advised that datasets from multiple sequencing runs be analyzed together because error-correction in dada2 is based on the assumption that different runs have different run-specific error profiles.  This is [discussed in more detail here](https://github.com/benjjneb/dada2/issues/1177).

## Outputs

The main outputs of the pipeline are:

1. [A surveillance report](#Surveillance-report)
2. [QC reports](#qc-reports)
3. Information about [observed sequences that were not assigned to any reference sequences](#BLASTing-of-unassigned-sequences)

### Output directory

Output files are placed by default in a `results` directory (or `test/results` when running with `-profile test`).  The output directory location can be overridden using the `--outdir` parameter.  For instance:

```
nextflow run stenglein-lab/tick_surveillance ... --outdir some_other_results_directory_name
````

### Output file name prefixes 

The main pipeline output file names will be prefixed by a value that is by default the date the pipeline is run (e.g. `2023_04_06_sequencing_report.xlsx`).  This filename prefix can be changed using the `--output_prefix` parameter.  For instance, running:

```
nextflow run stenglein-lab/tick_surveillance ... --output_prefix Run_XYZ
``` 

Will create a file named `Run_XYZ_sequencing_report.xlsx`

## Surveillance Report

The pipeline outputs a surveillance report in Microsoft Excel format that is named `<output_prefix>_sequencing_report.xlsx`.  

This first tab of this spreadsheet contains the main surveillance report table with positive/negative calls.  Other tabs contrain additional information, such as the number of reads assigned to each surveillance target, a copy of the input metadata, detailed information about specific reference sequence assignments, etc.

### Surveillance targets

Surveillance targets are defined in [the surveillance_columns file](./refseq/surveillance_columns.txt).  Each surveillance target corresponds to a column in the main surveillance report table.  

This is a 2-column tab-delimited file.  The first column contains the names of the columns that will form the surveillance report table.  The second column contains optional default text for this column.

It is possible to add or remove surveillance targets (columns in the surveillance table) by adding or removing them from this file. The path to a custom surveillance columns file can be specified using the `--surveillance_columns` parameter.

### Reference sequences

Reference sequences are sequences that are expected to be observed.  Observed sequences that are sufficiently similar to a reference sequences will be assigned to that reference sequence.

Reference sequence are defined in the [targets.tsv](./refseq/targets.tsv) file.  This tab-delimited file contains the following columns: 

| Column                 | Description |
| -----------            | ----------- |
| ref_sequence_name      | The reference sequence name |
| species                | The species for this reference sequence.  This value will be reported but is not used to map reference sequences to surveillance targets (reporting_columns is used for that).
| primer_name            | The name of the primers expected to amplify this target (e.g., FlaB).  Provided for reference only.
| reporting_columns      | A semicolon-separated list of surveillance targets.  Reads that are assigned to this reference sequence will be assigned to these targets, as described [below](#Mapping-of-reference-sequence-to-surveillance-targets). |
| min_percent_identity   | The minimum percent identity of the alignment between an observed sequence and this reference sequence to be assigned as a positive hit |
| min_percent_aligned    | The minimum percent of the observed sequence that must align to this reference sequence to be assigned as a positive hit |
| max_percent_gaps       | The maximum percent gap characters in alignments of observed sequences and this reference sequence to be assigned as a positive hit|
| internal_control       | True if this corresponds to an internal control target, such as tick actin or a "tick ID" amplicon|
| sequence               | The expected reference sequence, including primers |

Note that the names of the reporting_columns specified in this file must match exactly the names of columns defined in the surveillance columns definition file.

To add a new sequence to the targets.tsv file, you will need to edit this file.  It is a plain-text [tab-delimited file](https://en.wikipedia.org/wiki/Tab-separated_values) that can be edited in google sheets or Microsoft Excel or similar software.  Add a new row for the new referencesequence.  From google sheets, download the file in tab-separated value format and transfer it to the computer where you will be running this pipeline.  

The default location of the targets.tsv file can be overriden by specifying the --targets option on the nextflow command line.  For instance:

```
nextflow run main.nf -profile singularity --targets /path/to/targets.tsv
```

### Assignment of observed sequences to reference sequences

Dada2 outputs observed sequences known as amplicon sequence variants, or ASVs.  The pipeline uses BLASTN to align ASVs to the set of reference sequences defined in [targets.tsv](./refseq/targets.tsv).  ASVs that produce BLASTN alignments to a reference sequence that meet the percent identity and length criteria defined in the targets table will be assigned to that reference_sequence.  If an ASV produces alignments to more than one reference sequence, only the highest scoring alignment will be considered.


### Mapping of reference sequence to surveillance targets 

Surveillance targets (for instance Borrelia_burgdorferi_sensu_stricto) are defined in the [surveillance_columns.txt file](#Surveillance_columns_file).  Multiple reference sequences can map to a single surveillance target.  For instance, both the Bor_burgdorferi_B31 and the Bor_burgdorferi_N40 reference sequences map to the Borrelia_burgdorferi_sensu_stricto target.  The ASV counts for all mapped reference sequences are summed to produce the count for each surveillance target for the purpose of making positive calls.  These summed counts are output in the surveillance_counts tab of the main sequencing_report output speadsheet.

Just as multiple reference sequences can be mapped to one surveillance target, each reference sequence can be mapped to multiple surveillance targets. (This is a [many-to-many](https://en.wikipedia.org/wiki/Many-to-many_(data_model)) mapping).

The reporting_columns column in [targets.tsv](./refseq/targets.tsv) defines mapping of reference sequences to surveillance targets.  The value of this column takes the following format for mapping the reference sequence to one surveillance target:

`surveillance_target_name:[count|name]`

Multiple mappings for one reference sequence can be defined, seperated by semicolons.  For instance, the Bor_burgdorferi_B31 reference sequence is mapped to 2 surveillance targets:

`Borrelia_sp:count;Borrelia_burgdorferi_sensu_stricto:count`

This means that reads assigned to Bor_burgdorferi_B31 would be assigned to two surveillance targets: Borrelia_sp, a catch-all target for many different types of *Borrelia*, and Borrelia_burgdorferi_sensu_stricto, a target for *B. burgdorferi sensu stricto*.

**Count mapping**: As in the example above, count type mapping that counts assigned to those reference sequences will contribute to the summed count for that target.

**Name mapping**: Some surveillance targets report species names instead of read counts.  For instance, the Borrelia_Other_species_name target is meant to report the name of any *Borrelia* species that was called positive in that dataset.  Reference sequences can be mapped to both name and count type targets.   For example, the Bor_andersoni_BC_1_AF264908 reference sequence is mapped to two targets:

`Borrelia_sp:count;Borrelia_Other_species_name:name`

This means that reads counts for this reference sequence will be included in the catch-all Borrelia_sp target, and the species name *Borrelia_andersoni* will be added to the Borrelia_Other_species_name target, provided that reference sequence had enough reads to be called positive individually.  

### Calling positive hits

A key output of the pipeline is to define whether specific samples are positive for particular pathogens.  There are two ways to call a positive hit, depending on whether the surveillance target is a positive control target or not.

**Positive control targets**: these targets are defined as `internal_control` in [targets.tsv](./refseq/targets.tsv).  These are targets that are expected to amplify from any tick DNA sample (e.g tick actin, or a "tick ID" target).   These targets are called positive according to the following procedure:

- The mean number of reads for the target is calculated for each batch of datasets.  Samples can be binned into batches in [the metadata file](#Metadata-file) using the batch column.  If batches are not defined in the metadata, then all datasets will be binnned into a single batch.
- For each sample, if the number of reads assigned to this target is ≥ the mean value less 3 standard deviations, the sample is called positive for this target
- For all calculations, the log10(# reads) are used, because read counts for internal control targets exhibit log-normal distributions.

**All other targets**: These are any surveillance target that is not defined as an internal control in [targets.tsv](./refseq/targets.tsv).  If the summed ASV counts for a surveillance target are ≥50, the target will be called positive.  The value of 50 is a default cutoff that can be overridden using the `--min_reads_for_positive_surveillance_call` parameter on the nextflow command line.


## Primers and primer trimming

Primers are defined in the [primers.tsv](./refseq/primers.tsv) file.  Primer sequences defined in this file are used for two purposes:

1. To identify read pairs that have expected forward and reverse primers at their ends in the expected F/R orientation.  Only read pairs with a pair of expected primers in the expected orientation will be kept for further analysis.
2. Primer sequences will be trimmed off of observed sequences, since primer-derived sequences do not reliably reflect the template sequence (in case of amplification despite primer-template mismatches).

To add a new primer pair to the pipeline, you will need to edit this file.  It is a plain-text [tab-delimited file](https://en.wikipedia.org/wiki/Tab-separated_values) that can be edited in google sheets or similar software.  Add a new row for the new primer pair.  From google sheets, download the file in tab-separated value format and transfer it to the computer where you will be running this pipeline.  

The default location of the primers.tsv file can be overwritten by specifying the --primers option on the nextflow command line.  For instance:

```
nextflow run stenglein-lab/tick_surveillance -profile singularity --primers /path/to/primers.tsv
```

If primer sequences are not entered in the correct orientation, trimming will not work and targets amplified by these sequences will not be detected by the pipeline.  The solution in this case will most likely just be to swap the F/R orientation of the primers in this file.  

**Correct primer orientation:**  The forward primer (primer_f) should appear at the beginning of Illumina read 1 and be in the same orientation as R1.  The reverse primer (primer_r) should appear at the beginning of read 2 and be in the same orientation as read 2.  In other words, relative to the PCR product as a whole, the primers should point towards each other.

## QC reports

Two QC reports in HTML format from initial input sequence data (`...initial_qc_report.html`) and from sequence data after trimming low quality and adapter bases (`...post_trim_qc_report.html`) are output.

## BLASTing of unassigned sequences

It is possible that amplicon sequencing will generate sequences that are off-target or not closely related enough to be assigned to one of the predefined reference sequences.  The pipeline can BLAST these "unassigned" sequences against the NCBI nt database to try to figure out what they are.  

Enabling BLASTing of unassigned sequences is controlled by the `blast_unassigned_sequences` parameter, which is turned off by default.  

There are two ways to run the BLAST:

1. **Remote option:** By using BLASTn with the -remote option.  This is very slow but doesn't require a local copy of the BLAST database.  To run this way, set the `remote_blast_nt` parameter to true.  This sends the sequences to a remote NCBI server for BLASTing.
2. **Local option:** Alternatively, if you have a local copy of the nt BLAST database installed, you can specify its location using the `local_nt_database_dir` parameter, which should be the path to a directory containing a local nt database.  The name of this database is expected to be "nt", but this name can be changed by overriding the `local_nt_database_name` parameter.

These options can be configured on the command line, for example:

```
nextflow run stenglein-lab/tick_surveillance -profile test,singularity --blast_unassigned_sequences --remote_blast_nt true
```

Or in a nextflow config file, for instance:

```
  // profile for local BLAST of unassigned sequences
  local_blast {
    params.blast_unassigned_sequences = true
    params.local_nt_database_dir ="/home/NCBI_databases/nt/"
    params.remote_blast_nt = false
  }
```

The output from BLASTing unassigned sequences is contained in a file named `<run_prefix>_non_reference_sequence_assignments.xlsx`

## Dependencies

This pipeline has two main dependencies: nextflow and singularity.  These programs must be installed on your computer to run this pipeline.

### Nextflow

To run the pipeline you will need to be working on a computer that has nextflow installed. Installation instructions can be [found here](https://www.nextflow.io/docs/latest/getstarted.html#installation).  To test if you have nextflow installed, run:

```
nextflow -version
```

This pipeline requires nextflow version > 22.10

### Singularity 

The pipeline uses singularity containers to run programs like cutadapt, BLAST, and R.  To use these containers you must be running the pipeline on a computer that has [singularity](https://sylabs.io/singularity) [installed](https://sylabs.io/guides/latest/admin-guide/installation.html).  To test if singularity is installed, run:

```
singularity --version
```

There is no specified minimum version of singularity, but older versions of singularity (<~3.9) may not work.  the pipeline has been tested with singularity v3.9.5.  

Singularity containers will be automatically downloaded and stored in a directory named `singularity_cacheDir` in your home directory.  They will only be downloaded once.

### Conda

It is possible to run this pipeline using [conda](https://docs.conda.io/en/latest/) to handle dependencies. But it is strongly recommended to use singularity instead of conda.  

### R libraries

Some of the pipeline code is implemented in [R scripts](./scripts/).  Some of these scripts require R packages like [openxlsx](https://www.rdocumentation.org/packages/openxlsx/versions/4.2.5.2), for writing output in Excel format.  When running with singularity, these packages are installed locally, on top of a [Rocker tidyverse singularity image](https://rocker-project.org/images/).  This occurs in nextflow process `setup_R_dependencies`, which invokes [this script](./scripts/install_R_packages.R).

### Python dependencies

Some of the pipeline code is implemented in [Python scripts](./scripts/).  In particular, the tree-building scripts.  These python scripts require various python modules.  This is handled by creating a python virtual environment (venv), which happens in [nextflow process `setup_python_venv`](./modules/local/setup_python_env/main.nf).  Versions of pythons packages are defined [in this file](./lib/requirements.txt)).  This venv is then activated from a basic python singularity image (for instance [in process `CREATE_FASTA_FOR_TREES`](./subworkflows/generate_trees.nf)).


## Additional parameter information

For a complete description of pipeline parameters, see [this page](./parameters.md)
