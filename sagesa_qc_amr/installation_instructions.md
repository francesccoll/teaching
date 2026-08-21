
# 1. Install conda if you don't already have it

Check if conda is already installed on your system:

```
conda --version
mamba --version
```

If you don't already have conda installed, follow the next steps.

### [See the official install miniconda instructions](https://docs.conda.io/projects/miniconda/en/latest/) for Linux

OR

```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh
```

Follow prompts to accept license, choose install path, and allow the new bin directory to be added to .bashrc

### [See Install miniconda instructions](https://docs.conda.io/projects/miniconda/en/latest/) for MacOS

OR

```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash ./Miniconda3-latest-MacOSX-x86_64.sh
```

Follow prompts to accept license, choose install path, and allow the new bin directory to be added to .bashrc

We recommend you set up bioconda as the default channel, but it is not required:

```
export PATH=$HOME/miniconda3/bin:$PATH # Change to match installation location, if not default.
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

```
source ~/miniconda3/bin/activate
```

# 2. Install tools

## 2.1 Installing [fastq-dl](https://github.com/rpetit3/fastq-dl)

Create a mamba/conda environment with fastq-dl installed:
```
mamba create -n fastq-dl -c conda-forge -c bioconda fastq-dl
conda activate fastq-dl
fastq-dl --version
conda deactivate
```

## 2.2 Installing [fastqc](https://github.com/s-andrews/FastQC/)

Follow FastQC installation instructions [here](https://github.com/s-andrews/FastQC/blob/master/INSTALL.md).

The installation files can be downloaded from [here](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)

## 2.3 Installing [multiqc](https://github.com/multiqc/multiqc)

Create a mamba/conda environment with multiqc installed:
```
mamba create -n multiqc -c conda-forge -c bioconda multiqc
conda activate multiqc
multiqc --version
conda deactivate
```

## 2.4 Installing [nanoQC](https://github.com/wdecoster/nanoQC)

Create a mamba/conda environment with nanoqc installed:
```
mamba create -n nanoqc -c conda-forge -c bioconda nanoqc
conda activate nanoqc
nanoQC --version
conda deactivate
```

## 2.5 Installing [NanoPlot](https://github.com/wdecoster/nanoplot)

Create a mamba/conda environment with NanoPlot installed:
```
mamba create -n nanoplot -c conda-forge -c bioconda nanoplot
conda activate nanoplot
NanoPlot --version
conda deactivate
```

## 2.6 Installing [fastp](https://github.com/opengene/fastp)

Create a mamba/conda environment with fastp installed:
```
mamba create -n fastp -c bioconda fastp
conda activate fastp
fastp --version
conda deactivate
```

## 2.7 Installing [porechop](https://github.com/rrwick/porechop)

Create a conda environment called ```porechop``` and install the tool inside this environment:
```
conda create -n porechop
conda activate porechop
git clone https://github.com/rrwick/Porechop.git
cd Porechop
python3 setup.py install
porechop -h
```


filtlong

atb 

quast

checkm2 

sylph

amrfinder

lre-finder



