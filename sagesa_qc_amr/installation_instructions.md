
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

Before installing the tools below, navigate to a directory in your system where you would like to download the workshop data using the ```cd``` command. Next, clone this GitHub directory and navigate to the ```data``` subdirectory:

```
git clone https://github.com/francesccoll/teaching/
cd teaching/sagesa_qc_amr/
```

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
python3 setup.py install --prefix=$HOME/.local
cd ..
rm -rf Porechop
```

## 2.8 Installing [filtlong](https://github.com/rrwick/filtlong)

Filtlong builds into a stand-alone executable:
```
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong
make -j
bin/filtlong -h
```

It is recommended to copy filtlong executable to a directory in your PATH:
```
cp bin/filtlong /usr/local/bin
cd ..
rm -rf Filtlong
```

## 2.9 Installing [AllTheBacteria command-line tool](https://allthebacteria.org/docs/cli/)

On Linux or macOS:
```
curl -fsSL https://raw.githubusercontent.com/allthebacteria/atb-cli/main/install.sh | bash
```
Pre-built binaries for Linux, macOS, and Windows are available from: https://github.com/allthebacteria/atb-cli/releases/latest.

Choose a data directory on a disk with enough space. The default metadata fetch includes AMR data and builds local indexes; this can use roughly 35 GB. Add another few GB if you fetch the sketchlib index.
```
atb config set general.data_dir /path/to/large/volume/atb-data
atb fetch
```
 
## 2.10 Installing [quast](https://github.com/ablab/quast)

See System requirements and Installation here: https://github.com/ablab/quast

Or install Quast via conda:
```
mamba create -n quast -c conda-forge -c bioconda quast
conda activate quast
quast --version
conda deactivate
```

## 2.11 Installing [checkm2](https://github.com/chklovski/CheckM2)

Create a mamba/conda environment with checkm2 installed:
```
mamba create -n checkm2 -c bioconda -c conda-forge checkm2
conda activate checkm2
checkm2 --version
conda deactivate
```

You will also need to download and install the external DIAMOND database CheckM2 relies on for rapid annotation. Use checkm2 ```database --download``` to install it into a custom location using ```checkm2 database --download --path /custom/path/to/database/```. Replace ```/custom/path/to/database/``` with your chosen custom path. 

```
checkm2 database --download --path /custom/path/to/database
export CHECKM2DB="/custom/path/to/database"
```

## 2.12 Installing [sylph](https://github.com/bluenote-1577/sylph)

Create a mamba/conda environment with checkm2 installed:
```
mamba create -n sylph -c bioconda sylph
conda activate sylph
sylph --version
conda deactivate
```

Download the pre-built sylph database:
```
wget http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r220-c200-dbv1.syldb
```

## 2.13 Installing [amrfinder](https://github.com/ncbi/amr)

Create a mamba/conda environment with amrfinder installed, and download the amrfinder database:
```
mamba create -y -c conda-forge -c bioconda -n amrfinder --strict-channel-priority ncbi-amrfinderplus
conda activate amrfinder
mkdir ./data/amrfinder_db
amrfinder_update -d ./data/amrfinder_db
conda deactivate
```

## 2.15 Installing [lre-finder](https://bitbucket.org/genomicepidemiology/lre-finder/src/master/)

```
git clone https://bitbucket.org/genomicepidemiology/kma.git
cd kma && make
cp kma /usr/local/bin
cd ..
rm -rf kma
```

```
git clone https://bitbucket.org/genomicepidemiology/lre-finder.git
cd lre-finder && make
tar -xvzf elmDB.tar.gz
kma index -i elmDB/elm.fsa -o elmDB/elm
cp getGene LRE-Finder.py /usr/local/bin
```

## 2.15 Installing [Snippy](https://github.com/tseemann/snippy)

```
mamba install -c conda-forge -c bioconda snippy
snippy --version
snippy --check
conda deactivate
```




