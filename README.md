# CGR
A Data Processing Pipeline for Whole-genome Sequencing-Based Grapevine Genetic Testing using Chaos Game Representation and Integrative Deep Learning Approaches

Here we provide a list of scripts and python codes for processing grapevine whole genome sequence (WGS) data to be used grapevine species and cultivar identification using a combination of Chaos Game Representation (CGR) and Integrative Deep Learning Approaches. The untitlity of individual script/code is described below:

##perl scripts for aligning Illumina WGS data to the reference genome and for call variants, and then convert the variant genotype to pseudo-DNA sequences
1. fastq2cram2vcf4graham.sbatch: a slurm script for aligning WGS data in Illumina PE reads to a reference genome and generate cram and bcf files using high performacing computing.
2. VCF2FaSeq_v3.pl: a perl script for generating a pseudo fasta sequence to represent the genotype of variants of WGS samples for use in sequence similarity and phylogeny ananlysis

The remaining scripts are used in preprocessing sequence data into CGR images, generate structured training data, and fine-tune models. Similar code is split into separate files to facilitate job execution using a scheduler. Vitis Vinifera is processed separately due to significantly larger requirements resulting for handling the larger dataset available.

```plaintext
- cgr_by_chrom_v2.py
- cgr_by_chrom_v3.py
- cgr_by_species.py
- Finetune_test.py
- Finetune_vvcultivars.py
- generate_training_data.py
- generate_vvcultivar_training_data.py
- generate_wholespecies_training_data.py
- process_complete_sequences.py
- process_complete_sequences_nonvv.py
- process_folders_species.py
- process_sequences_otherspecies_v2.py
- process_sequences_species_file_v2.py
- SVM.py
- performance_analysis_Cs.txt
```

This organized into three main functionalities:
1. **Preprocessing**: Preprocessing sequence data into CGR files
2. **Training Data Generation**: Creating datasets for model training based on class balance, normalization, resolution, and CGR methodology.
3. **Fine-Tuning**: Fine-tuning models using the generated datasets with additional support for interactive execution.

## Script Details

### Preprocessing Scripts
Preprocessing scripts prepare the raw sequence data into a format suitable for training data generation. Each script targets specific categories to facilitate independent job scheduling:

- **`process_complete_sequences.py`**: Processes complete sequences for general use.
- **`process_complete_sequences_nonvv.py`**: Processes complete sequences excluding certain species.
- **`process_folders_species.py`**: Prepares species folders to facilitate organization
- **`process_sequences_otherspecies_v2.py`**: Handles sequences for non-vitis vinifera species in a specialized manner.
- **`process_sequences_species_file_v2.py`**: Processes sequences from species-specific files.
- **`cgr_by_*` Scripts**: Module for generating CGRs (Chaos Game Representations) based on specified parameters.

### Training Data Generation Scripts
These scripts generate structured datasets for model training:

- **`generate_training_data.py`**: Creates training data sets with customizable parameters like class balance and normalization.
- **`generate_vvcultivar_training_data.py`**: Generates training data specific to VV cultivars.
- **`generate_wholespecies_training_data.py`**: Prepares datasets covering whole species.

### Fine-Tuning Scripts
The fine-tuning scripts train models on the generated datasets:

- **`Finetune_test.py`**: General-purpose fine-tuning script. Can be run interactively or with specific parameters.
- **`Finetune_vvcultivars.py`**: Fine-tunes models specifically for VV cultivar datasets.

### SVM scripts
These scripts are used to perform kmer-based Support Vector Machines (SVMs) sample clustering
- ** `SVM.py` **: python code for kmer-based Support Vector Machines (SVMs) clustering
- ** `performance_analysis_Cs.txt` **: Modified code to perform analysis across varying values for C and write to txt file.

## Usage

Finetune_test.py includes an interactive script to facilitate testing in CLI.
Running each of the remaining python files with `--help` or no parameters to view required and optional parameters. For example:
```bash
python process_complete_sequences.py --help
```
