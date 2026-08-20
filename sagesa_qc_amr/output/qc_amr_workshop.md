# Workshop Pathogen Genomics and AMR: Quality Assessment and Interpretation Report

## Table Of Contents
1. [Illumina Raw Sequencing Data Quality](#illuminaqc)
2. [ONT Raw Sequencing Data Quality](#ontqc)
3. [Illumina Read Cleaning and Pre-processing](#illuminacleaning)
4. [ONT Read Cleaning and Pre-processing](#ontcleaning)
5. [Sequencing Depth and Coverage](#depth)
6. [Genome Assembly Quality](#assemblyqc)
7. [Contamination and Data Integrity](#contamination)
8. [AMR detection with AMRFinderPlus](#amrfinderplus)
9. [Interpretation of AMR reports: case studies](#casestudies)

   
## 1. Illumina Raw Sequencing Data Quality <a name="illuminaqc"></a>

### Step 1: Download fastq files

```
conda activate fastq-dl
fastq-dl --accession SRR9027862
```

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR902/002/SRR9027862/SRR9027862_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR902/002/SRR9027862/SRR9027862_2.fastq.gz
```

### Step 2: Running FastQC on Illumina fastqs

```
fastqc
```

### Step 3:  Opening the FastQC Report (GUI)


### Step 4:  Running MultiQC

```
cd output/multiqc
```

```
multiqc ./
```

## 2. ONT Raw Sequencing Data Quality <a name="ontqc"></a>

## 3. Illumina Read Cleaning and Pre-processing <a name="illuminacleaning"></a>

## 4. ONT Read Cleaning and Pre-processing <a name="ontcleaning"></a>

## 5. Sequencing Depth and Coverage <a name="depth"></a>

## 6. Genome Assembly Quality <a name="assemblyqc"></a>

## 7. Contamination and Data Integrity <a name="contamination"></a>

## 8. AMR detection with AMRFinderPlus <a name="amrfinderplus"></a>

## 9. Interpretation of AMR reports: case studies <a name="casestudies"></a>




