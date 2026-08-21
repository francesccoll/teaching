
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


