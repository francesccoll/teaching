# Workshop Pathogen Genomics and AMR: Quality Assessment and AMR Interpretation Report

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

### Step 1.1: Download fastq files

```fastq-dl``` is a command-line tool for downloading sequencing FASTQ files from the ENA (European Nucleotide Archive) or NCBI SRA. You give it an accession (e.g. a BioProject, sample, experiment, or individual sequencing run), and it looks up the associated runs and downloads their FASTQs. By default, it tries ENA first and can fall back to SRA. The usual output is .fastq.gz files plus metadata such as run-info.tsv. If you have an SRA/ENA accession and just want the corresponding FASTQ files without manually finding FTP URLs or dealing with SRA conversion yourself, fastq-dl is a convenient option.

Activate ```fastq-dl``` conda environment and download the fastq.gz files of the following run accessions.
```
conda activate fastq-dl
fastq-dl --accession ERR4095909 --outdir ./data/
```

NOTE: as an alternative, fastq.gz files can also be downloaded with ```wget``` using corresponding FASTQ FTP URLs:
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR409/009/ERR4095909/ERR4095909_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR409/009/ERR4095909/ERR4095909_2.fastq.gz
```

### Step 1.2: Running FastQC on Illumina fastqs

```
ls -lh ./data/ERR4095909*
```

You should see something like
<pre>
-rw-r-----+ 1 resh000371 csic108_res 109M Aug 23 11:42 ./data/ERR4095909_1.fastq.gz
-rw-r-----+ 1 resh000371 csic108_res 136M Aug 23 11:42 ./data/ERR4095909_2.fastq.gz
</pre>

Run the following command to run the fastqc tool on both the read files:
```
fastqc ./data/ERR4095909_1.fastq.gz ./data/ERR4095909_2.fastq.gz
```

Upon successful completion `fastqc` will create an analysis report in html format, one for each read file named after the name of the file you used as input. We can see the report by opening the html file in the web-browser.

### Step 1.3:  Opening the FastQC Report (GUI)




## 2. ONT Raw Sequencing Data Quality <a name="ontqc"></a>

## Step 2.1: Quality Control with NanoQC and NanoPlot

```
nanoQC ERR8282741.fastq.gz -o qc_output
```
```
NanoPlot --fastq ERR8282741.fastq.gz -o qc_output
```

## 3. Illumina Read Cleaning and Pre-processing <a name="illuminacleaning"></a>

## Step 3.1: Run fastp
```
fastp --in1 SRR9027862_1.fastq.gz --in2 SRR9027862_2.fastq.gz --out1 SRR9027862_1.trimmed.fastq.gz --out2 SRR9027862.trimmed.fastq.gz --length_required 40 --cut_front --cut_tail --cut_mean_quality 20
```

## 4. ONT Read Cleaning and Pre-processing <a name="ontcleaning"></a>

## Step 4.1: Adapter Trimming with Porechop
```
porechop -i ERR8282741.fastq.gz -o ERR8282741_trimmed.fastq.gz
```
```
filtlong --keep_percent 90 ERR8282741_trimmed.fastq.gz > ERR8282741_filtered.fastq
```
```
filtlong --min_length 1000 ERR8282741_filtered.fastq > ERR8282741_filtered_min1k.fastq
```

## 5. Sequencing Depth and Coverage <a name="depth"></a>

## 6. Genome Assembly Quality <a name="assemblyqc"></a>

### Step 6.1 Downloading Illumina genome assemblies from AllTheBacteria
```
atb info SAMN11579777
```
```
wget https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com/SAMN11579777.fa.gz
```
```
gunzip SAMN11579777.fa.gz
mv SAMN11579777.fa E_faecium.kerschner2019.Ef-11.fa
```

Count how many contigs the Illumina assembly has:
```
grep -c ">" E_faecium.kerschner2019.Ef-11.fa
```

### Step 6.2 Assessing assembly quality with QUAST
```
quast -h
```
```
quast \
E_faecium.VREN1631.illumina_assembly.fna \
E_faecium.VREN1631.hybrid_assembly.fna \
-o quast/
```

### Step 6.3 Assessing assembly quality with BUSCO markers


## 7. Contamination and Data Integrity <a name="contamination"></a>

### Step 7.1 Using CheckM2 to predict genome completeness
```
checkm2 predict --allmodels --lowmem --database_path {options.database} --remove_intermediates --threads 8
```

### Step 7.2 Using Sylph to identify taxonomic identity and contamination
```
sylph profile $GTDB -1 $fastq1 -2 $fastq2 -t 8 > $run_accession".profiling.tsv"
```

## 8. AMR detection with AMRFinderPlus <a name="amrfinderplus"></a>

```
amrfinder -n $assembly -O Enterococcus_faecium -o $out_file -d $amrfinder_db --threads 8
```

## 9. Interpretation of AMR reports: case studies <a name="casestudies"></a>

### 9.1 Identifying hetero-resistance and mixed infections from *Mycobacterium tuberculosis* genome data

### 9.2 Resolving AMR phenotype-genotype discrepancies from *Salmonella enterica* genomes

In this exercise we will AMRFinder to detect AMR from whole genome sequences of *Salmonella enterica* with a focus on cases of AMR genotype-phenotype discrepancies described in [this study](https://www.sciencedirect.com/science/article/pii/S0740002020301192).

The table below includes the phenotypic antibiogram of four *Salmonella enterica* strains.

| Strain ID | Ampicillin | Piperacillin | Amoxicillin/clavulanate | Chloramphenicol | Colistin | Fosfomycin | Trimethoprim/sulfamethoxazole | Amikacin | Gentamicin | Tobramycin | Tetracycline |
|-----------|------------|--------------|-------------------------|-----------------|----------|------------|-------------------------------|----------|------------|------------|--------------|
| B27 | >8 (R) | >16 (R) | >8/4 (R) | >8 (R) | >2 (R) | >32 (R) | >4/76 (R) | ≤8 (S) | >4 (R) | >4 (R) | >8 (R) |
| WD17 | >8 (R) | 16 (I) | ≤2/1 (S) | ≤8 (S) | ≤2 (S) | ≤32 (S) | ≤2/38 (S) | ≤8 (S) | ≤2 (S) | ≤2 (S) | ≤4 (S) |
| WT5 | >8 (R) | >16 (R) | >8/4 (R) | >8 (R) | ≤2 (S) | ≤32 (S) | >4/76 (R) | ≤8 (S) | ≤2 (S) | ≤2 (S) | >8 (R) |
| WW28B | >8 (R) | >16 (R) | 8/4 (S) | >8 (R) | ≤2 (S) | >32 (R) | >4/76 (R) | ≤8 (S) | ≤2 (S) | ≤2 (S) | >8 (R) |

And the table below reported cases of AMR genotype-phenotype discrepancies:

| Strain ID | Nature of Error | Genotype-Phenotype Disagreement |
|-----------|-----------------|---------------------------------|
| B27 | False negative | Amoxicillin/clavulanate resistance without any resistance determinant |
| B27 | False negative | Colistin resistance without any resistance determinant |
| B27 | False positive | Amikacin susceptibility despite presence of *aac(6′)-Iaa* gene |
| WD17 | False negative | Ampicillin resistance without any resistance determinant |
| WD17 | False positive | Amikacin and tobramycin susceptibility despite presence of *aac(6′)-Iaa* gene |
| WT5 | False negative | Amoxicillin/clavulanate resistance without any resistance determinant |
| WT5 | False positive | Amikacin and tobramycin susceptibility despite presence of *aac(6′)-Iaa* gene |
| WW28B | False negative | Fosfomycin resistance without any resistance determinant |
| WW28B | False positive | Amikacin and tobramycin susceptibility despite presence of *aac(6′)-Iaa* gene |

Now, run AMRFinder for all four *Salmonella* genomes by executing the commands below:

```
conda activate amrfinder

amrfinder -n ./data/Salmonella_enterica.B27.fna -O Salmonella -o Salmonella_enterica.B27.amrfinder.txt -d ./data/amrfinder_db/latest

amrfinder -n ./data/Salmonella_enterica.WD17.fna -O Salmonella -o Salmonella_enterica.WD17.amrfinder.txt -d ./data/amrfinder_db/latest

amrfinder -n ./data/Salmonella_enterica.WT5.fna -O Salmonella -o Salmonella_enterica.WT5.amrfinder.txt -d ./data/amrfinder_db/latest

amrfinder -n ./data/Salmonella_enterica.WW28B.fna -O Salmonella -o Salmonella_enterica.WW28B.amrfinder.txt -d ./data/amrfinder_db/latest
```

**❓ Question**
All four *Salmonella* strains are resistant to ampicillin according to phenotypic testing. Can you identify which resistance determinant(s) are responsible for phenotypic resistance to ampicillin in each strain? Do all four strains carry at least one resistance determinant that could explain ampicillin resistance?

**❓ Question**
Three *Salmonella* strains are resistant to trimethoprim/sulfamethoxazole. Can you identify the AMR genes responsible for resistance to this antibiotic combination?

**❓ Question**
Two *Salmonella* strains (B27 and WW28B) are resistant to fosfomycin. Can you identify the AMR gene(s) responsible for resistance to this antibiotic?

**❓ Question**
For false negatives, i.e. phenotypically resistant strains without any resistance determinant detected by AMRFinder, could you hypothesise what biological phenomenon may be responsible for such discrepancies?


### 9.3 Identifying linezolid-resistance point mutations from Enterococcus faecium genomes

Linezolid, an oxazolidinone antibiotic that inhibits protein synthesis by binding to the 50S subunit of the 23S rRNA, is considered a last-resort treatment option for vancomycin-resistant Enterococci (VRE). However, linezolid resistance has been reported in *E. faecium*, with the most prevalent resistance mechanisms involving mutations within the V domain of the 23S rRNA. In this exercise, we will analyse the genomes of four linezolid resistant VRE strains that originated from a [hospital outbreak in Austria](https://link.springer.com/article/10.1186/s13756-019-0598-z). We will evaluate the ability of different bioinformatic tools to detect linezolid resistance genetic markers from Illumina sequence data. The table below includes the accessions, isolate Ids, and linezolid MIC of the VRE genomes we will analyse:

| Illumina run accession | Isolate ID | Biosample accession | Linezolid MIC | Linezolid resistance |
| ---------------------- | ---------- | ------------------- | ------------- | -------------------- |
| SRR9027862             | Ef-39      | SAMN11579813        | 4             | R                    |
| SRR9027818             | Ef-07      | SAMN11579756        | 16            | R                    |
| SRR9027816             | Ef-09      | SAMN11579758        | 16            | R                    |
| SRR9027815             | Ef-11      | SAMN11579777        | >256          | R                    |

First, move to the working directory and list the genome assembly files we will analyse:
```
ls ./data/E_faecium.kerschner2019.*
```

You should see the four files listed as shown below:
<pre>
./data/E_faecium.kerschner2019.Ef-07.fa  ./data/E_faecium.kerschner2019.Ef-11.fa
./data/E_faecium.kerschner2019.Ef-09.fa  ./data/E_faecium.kerschner2019.Ef-39.fa
</pre>

These four Illumina assemblies were downloaded using AllTheBacteria command-line tool (```atb info``` command) as explained in section 6.1. Now, run the commands below to execute AmrFinderPlus on each of the four *E. faecium* Illumina assemblies:

Remember that the amrfinder database must be downloaded first:
```
mkdir ./data/amrfinder_db
conda activate amrfinder
amrfinder_update -d ./data/amrfinder_db
```

```
conda activate amrfinder
amrfinder -n ./data/E_faecium.kerschner2019.Ef-07.fa -O Enterococcus_faecium -o E_faecium.kerschner2019.Ef-07.amrfinder.txt -d ./data/amrfinder_db/latest
amrfinder -n ./data/E_faecium.kerschner2019.Ef-09.fa -O Enterococcus_faecium -o E_faecium.kerschner2019.Ef-09.amrfinder.txt -d ./data/amrfinder_db/latest
amrfinder -n ./data/E_faecium.kerschner2019.Ef-11.fa -O Enterococcus_faecium -o E_faecium.kerschner2019.Ef-11.amrfinder.txt -d ./data/amrfinder_db/latest
amrfinder -n ./data/E_faecium.kerschner2019.Ef-39.fa -O Enterococcus_faecium -o E_faecium.kerschner2019.Ef-39.amrfinder.txt -d ./data/amrfinder_db/latest
```

Now, pay particular attention to the output files of AMRFinder, those ending in ```.amrfinder.txt```. See if AMRFinder has found any genetic marker responsible for linezolid resistance in these genomes.

Next, we will used a specialized tool called [LRE-Finder](https://pubmed.ncbi.nlm.nih.gov/30863844/)), a tool developed to detect genetic markers of linezolid resistance in enterococci from whole-genome sequences. As this tools accepts raw Illumina reads as input, we will download the Illumina fastq.gz files using ```fastq-dl```:

```
conda activate fastq-dl

fastq-dl --accession SRR9027862
fastq-dl --accession SRR9027818
fastq-dl --accession SRR9027816
fastq-dl --accession SRR9027815

mv *.fastq.gz ./data/
```

After downloading the fastq files, run the following ```LRE-Finder.py``` commands:

```
LRE-Finder.py -ipe ./data/SRR9027862_1.fastq.gz ./data/SRR9027862_2.fastq.gz -o SRR9027862_lre-finder -t_db ./data/elmDB/elm -ID 80 -1t1 -cge -matrix
LRE-Finder.py -ipe ./data/SRR9027818_1.fastq.gz ./data/SRR9027818_2.fastq.gz -o SRR9027818_lre-finder -t_db ./data/elmDB/elm -ID 80 -1t1 -cge -matrix
LRE-Finder.py -ipe ./data/SRR9027816_1.fastq.gz ./data/SRR9027816_2.fastq.gz -o SRR9027816_lre-finder -t_db ./data/elmDB/elm -ID 80 -1t1 -cge -matrix
LRE-Finder.py -ipe ./data/SRR9027815_1.fastq.gz ./data/SRR9027815_2.fastq.gz -o SRR9027815_lre-finder -t_db ./data/elmDB/elm -ID 80 -1t1 -cge -matrix
```

Now, pay attention to the output messages of LRE-Finder on your screen, particularly the linezolid mutations found, and the "Wild type ratio", "Mutant type ratio", and "Predicted phenotype" of each mutation. Check if LRE-Finder has found any genetic marker responsible for linezolid resistance in these genomes, and compare the output of LRE-Finder with that of AMRFinder for the corresponding genomes.

**❓ Question**
Are there any discrepancies in the linezolid resistance markers found? if so, what could explain these discrepancies?

### 9.4 Effect of duplicated and truncated AMR genes on genotypic AMR determination

Tetracycline resistance in *E. faecium* is primarily mediated by ribosomal protection proteins encoded by *tet(M)* and efflux pumps encoded by *tet(L)*. *tet(M)* is among the most prevalent tetracycline resistance genes in this species and is often carried on mobile genetic elements (MGEs). The authors of a [recent study](https://www.biorxiv.org/content/10.64898/2026.06.08.729070v1.full) identified many false negatives when predicting tetracycline resistance using AMRFinder from Illumina genome assemblies. They noticed that many of these false negatives (i.e. tetracycline resistant strains by phenotypic AST lacking *tet* genes in their genome) actually carried partially detected *tet(M)* genes. In this exercise, we will interrogate Illumina and hybrid assemblies to investigate the source of tetracycline genotype-phenotype discrepancies.

First, we will use the original Illumina assemblies to identify resistance markers with AMRFinder:

```
conda activate amrfinder

amrfinder -n ./data/E_faecium.IVR212.illumina_assembly.fna -O Enterococcus_faecium -o E_faecium.IVR212.illumina_assembly.amrfinder.txt -d ./data/amrfinder_db/latest

amrfinder -n ./data/E_faecium.VREN0576.illumina_assembly.fna -O Enterococcus_faecium -o E_faecium.VREN0576.illumina_assembly.amrfinder.txt -d ./data/amrfinder_db/latest

amrfinder -n ./data/E_faecium.VREN1631.illumina_assembly.fna -O Enterococcus_faecium -o E_faecium.VREN1631.illumina_assembly.amrfinder.txt -d ./data/amrfinder_db/latest
```

**❓ Question**
Has AMRFinder identified any *tet* genes that could explain phenotypic tetracycline resistance of these strains?
```grep "tet" *illumina_assembly.amrfinder.txt```

Next, we will lower the minimum coverage to detect an AMR gene with AMRFinder down to 30% (--coverage_min 0.3) from the default 50% to see if any new AMR gene is identified by AMRFinder:

```
amrfinder -n ./data/E_faecium.IVR212.illumina_assembly.fna -O Enterococcus_faecium -o E_faecium.IVR212.illumina_assembly.amrfinder.0.3.txt -d ./data/amrfinder_db/latest --coverage_min 0.3 

amrfinder -n ./data/E_faecium.VREN0576.illumina_assembly.fna -O Enterococcus_faecium -o E_faecium.VREN0576.illumina_assembly.amrfinder.0.3.txt -d ./data/amrfinder_db/latest --coverage_min 0.3 

amrfinder -n ./data/E_faecium.VREN1631.illumina_assembly.fna -O Enterococcus_faecium -o E_faecium.VREN1631.illumina_assembly.amrfinder.0.3.txt -d ./data/amrfinder_db/latest --coverage_min 0.3 
```

**❓ Question**
Has AMRFinder identified any *tet* genes after lowering the minimum gene coverage threshold? In what coverage?
```grep "tet" *illumina_assembly.amrfinder.0.3.txt```

Finally, we will run AMRFinder on hybrid assemblies (closed genomes) produced by sequencing all three strains with both Illumina and ONT:

```
amrfinder -n ./data/E_faecium.IVR212.hybrid_assembly.fna -O Enterococcus_faecium -o E_faecium.IVR212.hybrid_assembly.amrfinder.txt -d ./data/amrfinder_db/latest

amrfinder -n ./data/E_faecium.VREN0576.hybrid_assembly.fna -O Enterococcus_faecium -o E_faecium.VREN0576.hybrid_assembly.amrfinder.txt -d ./data/amrfinder_db/latest

amrfinder -n ./data/E_faecium.VREN1631.hybrid_assembly.fna -O Enterococcus_faecium -o E_faecium.VREN1631.hybrid_assembly.amrfinder.txt -d ./data/amrfinder_db/latest
```

**❓ Question**
Has AMRFinder identified any *tet* genes in the hybrid assemblies? How many copies and in what coverage?
```grep "tet" *hybrid_assembly.amrfinder.txt```

**❓ Question**
Based on the results of AMRFinder from Illumina and hybrid assemblies, which results would you trust more? what could explain the lack of full *tet(M)* detection in the Illumina assemblies?

