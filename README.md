# LncRAnalyzer:
A NextFlow-based pipeline for lncRNAs and Novel Protein Coding Transcripts (NPCTs) identification using RNA-Seq

# Introduction:
LncRAnalyzer-nf is a NextFlow (DSL2) workflow to identify lncRNAs and Novel Protein Coding Transcripts (NPCT) using RNA-Seq. The pipeline contains several steps, including quality control, read alignment to the reference genome, reference-guided transcript assembly, merge annotations, annotation comparison, class code selection, and retrieval of transcripts in FASTA format. The putative class code selected transcripts will be further evaluated for their coding potentials, features, and protein domain homologies using CPC2, CPAT, RNAsamba, LncFinder, LGC, and PfamScan. The final lncRNAs and NPCTs will be selected based on coding potentials, features, and protein domain homologies. Additionally, if LiftOver files for the organism and related species are provided, this pipeline also performs cross-species lncRNA conservation analysis using Slncky. We also integrated the FEELnc plugin to report the mRNA spliced and intergenic lncRNAs in the given RNA-seq samples. For NPCTs, further functional annotations is needed, which include peptide sequences prediction using TransDecoder followed by homology searches using Pfamscan, BLASTP, and BLASTX. The entire workflow is automated and could be implemented in multiple working environments such as Conda, Docker, Singularity, Apptainer, Podman, and Wave.

<p align="center">
  <img src="https://gitlab.com/nikhilshinde0909/LncRAnalyzer-nf/raw/main/bin/LncRAnalyzer.png" width=50% height=25%>
</p>


# Implementation:

## Local environment
1. To execute the steps in the pipeline, download the latest release of LncRAnalyzer to your local system with the following command 
```
git clone https://gitlab.com/nikhilshinde0909/LncRAnalyzer-nf.git
```

2. Change directory to LncRAnalyzer.
```
cd LncRAnalyzer
```

3. Install the latest Miniforge and required software tools as follows
```
bash install_pipeline_tools.sh
```
This will install all required software and tools in the local environment.

4. Pipeline is ready to execute, prepare your inputs and data.groovy in the working directory

```
Working directory
mkdir data
├── data
│   ├── SRR975551_1.fastq.gz
│   ├── SRR975552_1.fastq.gz
│   └── (and other fastq.gz files)
│   ├── SRR975551_2.fastq.gz
│   ├── SRR975552_2.fastq.gz
│   └── (and other fastq.gz files)
│   └── hg38.rRNA.fasta
|   └── hg38.genome.fasta
|   └── hg38.annotation.gtf
|   └── (and other files)
└── data.groovy 
```  
Copy your RNA-seq reads (\*.fastq.gz), rRNA sequences (\*.fa), reference genomes (\*.fa), related sp. reference genome (\*.fa), annotations (\*.gtf) and liftover files in data directory; create file data.txt in the same by using data_template.txt and add paths for raw fastq.gz, rRNA sequences, reference genome, rel sp. reference genome, annotations and liftover files in the same. \
Note: Please refer to the data.groovy template in LncRAnalyzer-nf dir

5. If you don't have a reference genome, annotations, and rRNA sequence information, you can download the same with the script provided with the pipeline as follows
```
python check_ensembl.py org_name
eg. python find_species_in_ensembl.py Sorghum
> sbicolor
python ensembl.py org_name_in_ensembl
eg. python download_datasets_ensembl.py sbicolor
> Ensembl version 56 <- download the datasets
```

6. Similarly, if you don't have liftover files for conservation analysis, then you can generate it through genome alignments of reference and query species genomes as follows
```
python Liftover.py <threads> <genome> <org_name> <genome_related_species> <rel_sp_name> <params_distance>
eg.
python Liftover.py 16 Sorghum_bicolor.dna.toplevel.fa Sbicolor Zea_mays.dna.toplevel.fa Zmays near
```
We also provide an additional script which will take ensembl gtf and produce bed files to run Slncky as follows
```
python ensembl_gtf2bed.py <ensembl_gtf> <output_prefix>
eg.
python ensembl_gtf2bed.py Sorghum_bicolor.58.gtf Sorghum_bicolor
```
This will produce protein-coding, non-coding, mirRNA, and snoRNA bed files for Slncky.

7. The pipeline is ready for execution \
Run the following command to print help
```
nextflow run ~/Path_to_LncRAnalyzer-nf/main.nf --help
```

8. Search for supporting species and configure it in the data.groovy
```
nextflow run ~/Path_to_LncRAnalyzer-nf/main.nf --supporting_species
```

9. Execute pipeline with RNA-seq datasets as follows
```
nextflow run ~/Path_to_LncRAnalyzer-nf/main.nf -c data.groovy -profile standard --threads 16 --memory 40.GB
```

## Containerized environment
1. Run the following commands to execute LncRAnalyzer-nf with pre-built docker image
```
nextflow run ~/Path_to_LncRAnalyzer-nf/main.nf -c data.groovy -profile docker --threads 16 --memory 40.GB
```
2. Run the following commands to execute LncRAnalyzer-nf with pre-built singularity image
```
nextflow run ~/Path_to_LncRAnalyzer-nf/main.nf -c data.groovy -profile singularity --threads 16 --memory 40.GB
```
3. Run the following commands to execute LncRAnalyzer-nf with podman
```
nextflow run ~/Path_to_LncRAnalyzer-nf/main.nf -c data.groovy -profile podman --threads 16 --memory 40.GB
```
4. Run the following commands to execute LncRAnalyzer-nf with apptainer
```
nextflow run ~/Path_to_LncRAnalyzer-nf/main.nf -c data.groovy -profile apptainer --threads 16 --memory 40.GB
```
5. Run the following commands to execute LncRAnalyzer-nf with wave
```
nextflow run ~/Path_to_LncRAnalyzer-nf/main.nf -c data.groovy -profile wave --threads 16 --memory 40.GB
```

## Thanks for using LncRAnalyzer-nf

# Cite:
Nikhil, S., Shaik Mohideen, H., & Natesan Sella, R. (2025). LncRAnalyzer: a robust workflow for long non-coding RNA discovery using RNA-Seq. The Plant journal : for cell and molecular biology, 124(1), e70509. https://doi.org/10.1111/tpj.70509

# Performance:
The performance of coding potential prediction using CPAT, CPC2, LGC, RNAsamba, and FEELnc was estimated with 50 RNA-Seq accessions of sorghum cultivar PR22 from past studies [https://doi.org/10.1186/s12864-019-5734-x] 

<p align="center">
  <img src="https://gitlab.com/nikhilshinde0909/LncRAnalyzer/raw/main/bin/ROC.png" width=70% height=70%>
</p>
