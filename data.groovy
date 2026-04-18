params {

// input reads option
reads_R1 = "data/SRR*_1.fastq.gz" // Note: if reads args mentioned here then remove it from command line
reads_R2 = "data/SRR*_2.fastq.gz" 

// rRNA sequences
rRNAs = "data/Sorghum.rRNAs.fa"

// Organism name
org_name = "Sorghum_bicolor" 
clade = "plants"

// The genome and annotation
genome = "data/Sorghum_bicolor.dna.toplevel.fa"
annotation = "data/Sorghum_bicolor.protein_coding.gtf"
liftover = "data/SbicolortoZmays.over.chain.gz"
noncoding = "data/Sorghum_bicolor.noncoding.bed"
mir = "data/Sorghum_bicolor.mir.bed"
sno = "data/Sorghum_bicolor.sno.bed"
known_lncRNAs_FA = "data/Sorghum_bicolor.PLncDB.fa"  // Optional, In case the organism is not listed in supporting species

// related species name
rel_sp_name = "Zea_mays"

//The genome and annotation of a related species
genome_related_species = "data/Zea_mays.dna.toplevel.fa"
annotation_related_species = "data/Zea_mays.protein_coding.gtf"
rel_liftover = ""
rel_noncoding = "data/Zea_mays.noncoding.bed"
rel_mir = "data/Zea_mays.mir.bed"
rel_sno = "data/Zea_mays.sno.bed"

// Experimental design file 
design = "" // Optional, Provide tab-separated exp design file with columns Run and Conditions
}
