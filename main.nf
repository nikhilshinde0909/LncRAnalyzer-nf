#!/usr/bin/env nextflow
nextflow.enable.dsl=2
VERSION="1.00"

if (params.help) {
    log.info "----------------------------------------------------------------"
    log.info "                      USAGE                                     "
    log.info "----------------------------------------------------------------"
    log.info ""
    log.info "nextflow main.nf -c <config_file>"
    log.info "Please refer data.groovy for config file"
    log.info ""
    log.info "Arguments:"
    log.info "--reads_R1                Files            Forword reads fastq.gz"
    log.info "--reads_R2                Files            Reverse reads fastq.gz"
    log.info "--threads                 Integer          threads (DEFAULT = 16)"
    log.info "--memory                 Integer         memory (DEFAULT = 20.GB)"
    log.info "Flags:"
    log.info "--supporting_species      List            Show supporting species"
    log.info "--version                                           Print version"
    log.info "--help                                       Display this message"
    log.info ""   
    exit 1
}

params.version = null
if (params.version) {
    log.info "Version: $VERSION"
    exit 1
}

params.supporting_species = null
if (params.supporting_species) {
    log.info "-----------------------------------------------------------------"
    log.info "                   SUPPORTING SPECIES                            "
    log.info "-----------------------------------------------------------------"
    log.info "Plants:"
    log.info "Actinidia_chinensis     Amborella_trichopoda       Ananas_comosus"
    log.info "Aquilegia_coerulea    Arabidopsis_lyrata     Arabidopsis_thaliana"
    log.info "Arachis_ipaensis   Brachypodium_distachyon	Brassica_napus"
    log.info "Brassica_rapa          Capsella_rubella            Capsicum_annum"
    log.info "Carica_papaya     Chlamydomonas_reinhardtii       Cicer_arietinum"
    log.info "Citrus_clementina        Citrus_maxima             Citrus_sinesis"
    log.info "Coffea_arabica           Cucumis_sativus            Daucus_carota"
    log.info "Durio_zibethinus       Elaeis_guineensis      Erythranthe_guttata"
    log.info "Eucalyptus_grandis       Eutrema_salsugineum       Fragaria_vesca"
    log.info "Glycine_max        Gossypium_barbadense       Gossypium_raimondii"
    log.info "Hordeum_vulgare     Jatropha_curcas               Malus_domestica"
    log.info "Musa_acuminata          Nicotina_tabacum             Oryza_sativa"
    log.info "Pisum_sativum         Sorghum_bicolor           Triticum_aestivum"
    log.info "Vitis_vinifera            Zea_mays                               "
    log.info ""
    log.info "Vertebrates:"
    log.info "Bos_taurus          Canis_lupus_familiaris           Capra_hircus"
    log.info "Danio_rerio            Equus_caballus               Gallus_gallus"
    log.info "Gorilla_gorilla      Heterocephalus_glaber           Homo_sapiens"
    log.info "Macaca_mulatta         Mus_musculus                     Naja_naja"
    log.info "Ovis_aries             Pan_troglodytes               Pongo_abelii"
    log.info "Rattus_norvegicus        Salmo_salar                   Sus_scrofa"
    log.info ""
    log.info "Invertebrates:                                                    "
    log.info "Drosophila_melanogaster                                          "
    exit 1
}

// define color scheme
ANSI_RESET = "\u001B[0m";
ANSI_RED = "\u001B[31m";
ANSI_GREEN = "\u001B[32m";
def print_red = {  str -> ANSI_RED + str + ANSI_RESET }
def print_green = {  str -> ANSI_GREEN + str + ANSI_RESET }

if (params.reads_R1 && params.rRNAs && params.genome && params.annotation) {
    log.info print_green("Input reads and configuration file is provided via -c ")
    log.info print_green("Starting analysis with LncRAnalyzer-nf version: $VERSION")
    } else {
    log.info print_red("Input data configuration file is not provided")
    log.info "----------------------------------------------------------------"
    log.info "                      USAGE                                     "
    log.info "----------------------------------------------------------------"
    log.info ""
    log.info "nextflow main.nf -c <config_file>"
    log.info "Please refer data.groovy for config file"
    log.info ""
    log.info "Arguments:"
    log.info "--reads_R1                Files            Forword reads fastq.gz"
    log.info "--reads_R2                Files            Reverse reads fastq.gz"
    log.info "--threads                 Integer          threads (DEFAULT = 16)"
    log.info "--supporting_species      List            Show supporting species"
    log.info "--memory                 Integer         memory (DEFAULT = 20.GB)"
    log.info "Flags:"
    log.info "--version                                           Print version"
    log.info "--help                                       Display this message"
    exit 1
}

cpu_free = params.threads
mem_free = params.memory

if(params.rRNAs && params.genome && params.annotation){
log.info """\
         R U N N I N G   A N A L Y S I S   O N   I N P U T S   B E L O W
         ==================================================================
         Species: ${params.org_name}
         Clade : ${params.clade}
         rRNAs: ${params.rRNAs}
         Genome : ${params.genome}
         Annotations : ${params.annotation}
         Related sp: ${params.rel_sp_name}
         Genome related sp.: ${params.genome_related_species}
         Annotations related sp.: ${params.annotation_related_species}
         CPU resources: ${cpu_free}
         Mem resorces: ${mem_free}
         ==================================================================
         """
         .stripIndent()
    } else {
    log.info "ERROR: Missing required input! You must provide rRNAs, genome, and annotations"
    exit 1
}


if (!params.containsKey('reads_R2')) {
    params.reads_R2 = null
}

isPaired = params.reads_R2 ? true : false

read_pairs_ch = isPaired ? 
    Channel.fromFilePairs([params.reads_R1, params.reads_R2], checkIfExists: true).map { id, files -> [id.replaceAll(/_[12]$/, ""), files] } : 
    Channel.fromPath(params.reads_R1, checkIfExists: true).map { f -> [f.simpleName.replaceAll(/_[12]$/, ""), [f]] }


process run_fastp {
    tag "Running FastP on $pair_id"
    cpus cpu_free
    publishDir "FastP", mode: 'copy', pattern: '*.{fastq.gz,html,json}'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*.json"),     emit: json_ch
    tuple val(pair_id), path("*.html"),     emit: html_ch
    tuple val(pair_id), path("*.fastq.gz"), emit: trimmed_fastq_ch

    script:
    if (reads.size() == 2) {
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${pair_id}_trimmed_R1.fastq.gz -O ${pair_id}_trimmed_R2.fastq.gz -h ${pair_id}.html -j ${pair_id}.json
    """
    } else {
    """
    fastp -i ${reads[0]} -o ${pair_id}_trimmed.fastq.gz -h ${pair_id}.html -j ${pair_id}.json 
    """
    }
}

rRNA_input_ch = Channel.fromPath(params.rRNAs, checkIfExists: true)

process rRNA_index {
    tag "Indexing rRNAs"
    publishDir "rRNA_unmapped", pattern: 'rRNA.*', mode: 'copy', overwrite: true
    cpus (cpu_free/4)
    
    input:
    path(rRNA_fa)

    output:
    path("rRNA*"), emit: rRNA_index_ch

    script: 
    """
    hisat2-build -p ${task.cpus} ${rRNA_fa} rRNA
    """
}

process rRNA_unmapped_reads {
    tag "Getting rRNA unmapped reads $pair_id"
    publishDir "rRNA_unmapped", mode: 'copy', pattern: '*.{fastq.gz,summary}'
    cpus cpu_free

    input:
    tuple val(pair_id), path(trimmed)
    path(rRNA_index_files) 

    output:
    tuple val(pair_id), path("${pair_id}_unmapped*.fastq.gz"), emit: clean_reads_ch
    tuple val(pair_id), path("${pair_id}.rRNA.summary"),     emit: rRNA_summary_ch

    script:
    def index_base = "rRNA"
    if (trimmed.size() == 2) {
    """
    hisat2 -p ${task.cpus} --summary-file ${pair_id}.rRNA.summary -x ${index_base} -1 ${trimmed[0]} -2 ${trimmed[1]} | samtools view -f 12 -Su - | samtools sort -n - -o ${pair_id}.u.bam
    bamToFastq -i ${pair_id}.u.bam -fq ${pair_id}_unmapped_1.fastq -fq2 ${pair_id}_unmapped_2.fastq
    gzip ${pair_id}_unmapped_1.fastq ${pair_id}_unmapped_2.fastq
    """
    } else {
    """
    hisat2 -p ${task.cpus} --summary-file ${pair_id}.rRNA.summary -x ${index_base} -U ${trimmed[0]} | samtools view -f 4 -Su -| samtools sort -n - -o ${pair_id}.u.bam
    bamToFastq -i ${pair_id}.u.bam -fq ${pair_id}_unmapped.fastq
    gzip ${pair_id}_unmapped.fastq
    """
    }
}

annotation_input_ch = Channel.fromPath(params.annotation, checkIfExists: true)

process extract_splicesites {
    tag "Known splicesites"
    publishDir "Align_assembly", pattern: 'known_splicesites.*', mode: 'copy', overwrite: true
    cpus 4
    
    input:
    path(annotation_GTF)

    output:
    path("known_splicesites.txt"), emit: known_splicesites_ch

    script: 
    """
    extract_splice_sites.py ${annotation_GTF} > known_splicesites.txt
    """
}

genome_input_ch = Channel.fromPath(params.genome, checkIfExists: true)

process genome_index {
    tag "Genome index"
    publishDir "Align_assembly", pattern: 'genome.*', mode: 'copy', overwrite: true
    cpus (cpu_free/2)
    
    input:
    path(genome_fa)

    output:
    path("genome*"), emit: genome_index_ch

    script: 
    """
    hisat2-build -p ${task.cpus} ${genome_fa} genome
    """
}

process reference_alignment {
    tag "Reference alignment $pair_id"
    publishDir "Align_assembly", mode: 'copy', pattern: '*.{bam,summary}'
    cpus cpu_free

    input:
    tuple val(pair_id), path(rRNA_cleaned)
    path(genome_index_files) 
    path(known_splicesites_file)

    output:
    tuple val(pair_id), path("${pair_id}.bam"), emit: aligned_bam_ch
    tuple val(pair_id), path("${pair_id}.summary"),     emit: alignment_summary

    script:
    def index_base = "genome"
    if (rRNA_cleaned.size() == 2) {
    """
    hisat2 --dta --mp 2 -p ${task.cpus} --known-splicesite-infile ${known_splicesites_file} --summary-file ${pair_id}.summary -x ${index_base} -1 ${rRNA_cleaned[0]} -2 ${rRNA_cleaned[1]} |samtools view -Su - | samtools sort - -o ${pair_id}.bam
    """
    } else {
    """
    hisat2 --dta --mp 2 -p ${task.cpus} --known-splicesite-infile ${known_splicesites_file} --summary-file ${pair_id}.summary -x ${index_base} -U ${rRNA_cleaned[0]} |samtools view -Su - | samtools sort - -o ${pair_id}.bam
    """
    }
}

process reference_guided_assembly {
    tag "Reference guided assembly $pair_id"
    publishDir "Align_assembly", mode: 'copy', pattern: '*.nc.gtf'
    cpus (cpu_free)

    input:
    tuple val(pair_id), path(reference_alignment)

    output:
    tuple val(pair_id), path("${pair_id}.nc.gtf"), emit: assembly_GTF_ch

    script:
    """
    stringtie -p ${task.cpus} ${reference_alignment} -m 200 -a 10 --conservative -g 50 -u -c 3 -o ${pair_id}.nc.gtf
    """
}

process merged_GTF {
    tag "Merging gtf files"
    publishDir "Align_assembly", mode: 'copy', pattern: '*.{merged.gtf,txt}'
    cpus 4

    input:
    path(nc_GTFs)
    
    output:
    path("genome.merged.gtf"), emit: GTF_merged_ch
    path("GTFs_list.txt"), emit: GTFs_list_ch

    script:
    """
    for file in ${nc_GTFs}; do echo \$file >> GTFs_list.txt; done
    stringtie --merge -m 200 -c 3 -T 1 -o genome.merged.gtf GTFs_list.txt
    """
}

process size_selected_GTF {
    tag "Selecting transcript size"
    publishDir "Annotation_compare", mode: 'copy', pattern: '*.size_selected.gtf'
    cpus 4

    input:
    path(merged_GTF_file)

    output:
    path("genome.size_selected.gtf"), emit: GTF_size_select_ch

    script:
    """
    gffread ${merged_GTF_file} -l 200 -U -T -o genome.size_selected.gtf
    """
}

process compare_annotations {
    tag "Comparing annotations"
    publishDir "Annotation_compare", mode: 'copy', pattern: 'compare.*'
    cpus 4

    input:
    path(size_selected_file)
    path(annotation_GTF)

    output:
    path("compare.annotated.gtf"), emit: annotation_compare_ch

    script:
    def index_base = "compare"
    """
    gffcompare ${size_selected_file} -r ${annotation_GTF} -o ${index_base}
    """
}

process classcode_select {
    tag "Selecting transcript classcodes"
    publishDir "Putative_lncRNA-NPCTs", mode: 'copy', pattern: 'Putative_lncRNA-NPCT.*'
    cpus 4

    input:
    path(compared_GTF)
    path(genome_fa)

    output:
    path("Putative_lncRNA-NPCT.list"), emit: lnc_npcts_list_ch
    path("Putative_lncRNA-NPCT.gtf"), emit: lnc_npcts_GTF_ch
    path("Putative_lncRNA-NPCT.fa"), emit: lnc_npcts_FA_ch

    script:
    """
    python $baseDir/bin/classcodes_select.py ${compared_GTF} Putative_lncRNA-NPCT.list
    python $baseDir/bin/subset_gtf.py ${compared_GTF} Putative_lncRNA-NPCT.list Putative_lncRNA-NPCT.gtf
    gffread Putative_lncRNA-NPCT.gtf -g ${genome_fa} -w Putative_lncRNA-NPCT.fa
    """
}

process perform_cpc2 {
    tag "Runnging CPC2"
    publishDir "CPC2", mode: 'copy', pattern: '*cpc2.*'
    cpus 4
    
    input:
    path(lncRNA_npcts_fa)

    output:
    path("Putative.lnc_NPCTs.cpc2.txt"), emit: cpc2_out_ch
    path("final_lncRNAs.cpc2.list"), emit: cpc2_lncrnas_ch
    path("final_NPCTs.cpc2.list"), emit: cpc2_npcts_ch

    script:
    """
    conda run -n cpc2-cpat-slncky CPC2.py -i ${lncRNA_npcts_fa} -o Putative.lnc_NPCTs.cpc2.txt
    grep -E -w 'noncoding' Putative.lnc_NPCTs.cpc2.txt |cut -f 1 > final_lncRNAs.cpc2.list
    grep -E -w 'coding' Putative.lnc_NPCTs.cpc2.txt |cut -f 1 > final_NPCTs.cpc2.list
    """
}

process perform_lgc {
    tag "Runnging LGC"
    publishDir "LGC", mode: 'copy', pattern: '*lgc.*'
    cpus 4
    
    input:
    path(lncRNA_npcts_fa)

    output:
    path("Putative.lnc_NPCTs.lgc.txt"), emit: lgc_out_ch
    path("final_lncRNAs.lgc.list"), emit: lgc_lncrnas_ch
    path("final_NPCTs.lgc.list"), emit: lgc_npcts_ch

    script:
    """
    $baseDir/utils/lgc-1.0.py ${lncRNA_npcts_fa} Putative.lnc_NPCTs.lgc.txt
    sed 1,11d Putative.lnc_NPCTs.lgc.txt | grep -w 'Non-coding' |cut -f1 > final_lncRNAs.lgc.list
    sed 1,11d Putative.lnc_NPCTs.lgc.txt | grep -w 'Coding' |cut -f1 > final_NPCTs.lgc.list
    """
}

hexamer_table = "${baseDir}/Models/CPAT/${params.org_name}_hexamer.TSV"
logit_model = "${baseDir}/Models/CPAT/${params.org_name}.logit.RData"
cutoff_file = "${baseDir}/Models/CPAT/CPAT_cutoffs.TSV"

models_exist = file(hexamer_table).exists() && file(logit_model).exists()

if (models_exist) {
    log.info "CPAT models are exists for species, Mode: Inference Only."
    ch_gtf_for_cds = Channel.empty()
    ch_hexamer_final = Channel.fromPath(hexamer_table)
    ch_logit_final = Channel.fromPath(logit_model)
    ch_feature_for_cutoff = Channel.fromPath(cutoff_file) 
    known_lncRNA_ch = Channel.empty()
    } else {
    log.info "CPAT models not exists for species, Mode: Training + Inference."
    ch_gtf_for_cds = Channel.fromPath(params.annotation)
    ch_hexamer_final = Channel.empty() 
    ch_logit_final = Channel.empty()
    ch_feature_for_cutoff = Channel.empty()
    known_lncRNA_ch = Channel.fromPath(params.known_lncRNAs_FA, checkIfExists: true)
}

process extract_cds {
    tag "Extracting CDS"
    publishDir "CPAT", pattern: '*.cds.fa', mode: 'copy', overwrite: true
    cpus 4

    input: 
    path(annotation_GTF)
    path(genome_FA)
    
    output:
    path("${params.org_name}.cds.fa"), emit: cds
    
    script: 
    """
    gffread ${annotation_GTF} -g ${genome_FA} -x ${params.org_name}.cds.fa
    """
}


process build_hexamer_table {
    tag "Buiding hexamer table"
    publishDir "CPAT", pattern: '*_hexamer.TSV', mode: 'copy', overwrite: true
    
    input: 
    path(cds)
    path(known_lnc_fa)
    
    output: 
    path("${params.org_name}_hexamer.TSV")
    
    script: 
    """
    conda run -n cpc2-cpat-slncky make_hexamer_tab.py -c ${cds} -n ${known_lnc_fa}| grep -v '^[[:space:]]*\$' > ${params.org_name}_hexamer.TSV
    """
}

process build_logit_model {
    tag "Building logit model"
    publishDir "CPAT", pattern: '*.logit.RData', mode: 'copy', overwrite: true
    
    input:
    path(hexamer)
    path(cds)
    path(known_lnc_fa)
    
    output: 
    path("${params.org_name}.logit.RData"), emit: model
    path("${params.org_name}.feature.xls"), emit: feature
    
    script: 
    """
    conda run -n cpc2-cpat-slncky make_logitModel.py -x $hexamer -c $cds -n $known_lnc_fa -o ${params.org_name}
    """
}

process run_CPAT {
    tag "Running CPAT"
    publishDir "CPAT", pattern: '*.cpat.TSV', mode: 'copy', overwrite: true
    
    input:
    path(hexamer)
    path(logit)
    path(query_fa)
    
    output:
    path("Putative.lnc_NPCTs.cpat.TSV"), emit: cpat_output_ch
    
    script:
    """
    conda run -n cpc2-cpat-slncky cpat.py -x $hexamer -g $query_fa -d $logit -o Putative.lnc_NPCTs.cpat.TSV
    """
}

process get_cutoff {
    tag "Getting CPAT cutoffs"
    publishDir "CPAT", pattern: 'cutoff_value', mode: 'copy', overwrite: true
    publishDir "CPAT", pattern: '*.pdf', mode: 'copy', overwrite: true
    input: 
    path(input_file)
    
    output: 
    path("cutoff_value"), emit: cpat_cutoff_ch
    path("*.pdf"), optional: true
    
    script:
    if (models_exist) {
    """
    grep ${params.org_name} $input_file | cut -f2 > cutoff_value 
    """
    } else {
    """
    Rscript $baseDir/bin/10fold_crossval.R $input_file ${params.org_name}.pdf ${params.org_name} cutoff_value
    """
    }
}

process get_cpat_results {
    tag "Getting CPAT results"
    publishDir "CPAT", pattern: 'final_*.cpat.list', mode: 'copy', overwrite: true

    input:
    path(CPAT_output)
    path(cutoff)

    output:
    path("final_lncRNAs.cpat.list"), emit: cpat_lncrnas_ch
    path("final_NPCTs.cpat.list"), emit: cpat_npcts_ch

    script:
    """
    Rscript $baseDir/bin/get_CPAT_results.R ${CPAT_output} ${cutoff} final_lncRNAs.cpat.list final_NPCTs.cpat.list
    """
}

rnasamba_model = "${baseDir}/Models/rnasamba/${params.org_name}.hdf5"

models_exist = file(rnasamba_model).exists()

if (models_exist) {
    log.info "RNAsamba models are exists for species, Mode: Inference Only."
    ch_gtf_for_mRNAs = Channel.empty()
    ch_model_rnasamba_final = Channel.fromPath(rnasamba_model)
    known_lncRNA_ch = Channel.empty()
    } else {
    log.info "RNAsamba models not exists for species, Mode: Training + Inference."
    ch_gtf_for_mRNAs = Channel.fromPath(params.annotation)
    ch_model_rnasamba_final = Channel.empty() 
    known_lncRNA_ch = Channel.fromPath(params.known_lncRNAs_FA, checkIfExists: true)
}

process extract_mRNAs {
    tag "Extracting spliced mRNA sequences"
    publishDir "RNAsamba", pattern: '*.mRNAs.fa', mode: 'copy', overwrite: true
    cpus 4

    input: 
    path(annotation_GTF)
    path(genome_FA)
    
    output:
    path("${params.org_name}.mRNAs.fa"), emit: mRNAs
    
    script: 
    """
    gffread ${annotation_GTF} -g ${genome_FA} -x ${params.org_name}.mRNAs.fa
    """
}


process rnasamba_train {
    tag "Buiding RNAsamba models"
    publishDir "RNAsamba", pattern: '*.hdf5', mode: 'copy', overwrite: true
    
    input: 
    path(mRNA_sequeces)
    path(known_lnc_fa)
    
    output: 
    path("${params.org_name}.hdf5"), emit: model_rnasamba_ch
    
    script: 
    """
    conda run -n rnasamba rnasamba train -s 15 -e 20 -v 3 ${params.org_name}.hdf5 ${mRNA_sequeces} ${known_lnc_fa}
    """
}


process rnasamba_classify {
    tag "Running RNAsamba classify"
    publishDir "RNAsamba", pattern: '*.rnasamba.TSV', mode: 'copy', overwrite: true
    
    input:
    path(RNAsamba_model)
    path(putative_transcripts_fa)
    
    output:
    path("Putative.lnc_NPCTs.rnasamba.TSV"), emit: rnasamba_output_ch
    
    script:
    """
    conda run -n rnasamba rnasamba classify Putative.lnc_NPCTs.rnasamba.TSV ${putative_transcripts_fa} ${RNAsamba_model} 
    """
}


process get_rnasamba_results {
    tag "Getting RNAsamba results"
    publishDir "RNAsamba", pattern: 'final_*.rnasamba.list', mode: 'copy', overwrite: true
    
    input:
    path(RNAsamba_output)

    output:
    path("final_lncRNAs.rnasamba.list"), emit: rnasamba_lncrnas_ch
    path("final_NPCTs.rnasamba.list"), emit: rnasamba_npcts_ch
    
    script:
    """
    grep -w 'noncoding' ${RNAsamba_output}|cut -f1 > final_lncRNAs.rnasamba.list
    grep -w 'coding' ${RNAsamba_output}|cut -f1 > final_NPCTs.rnasamba.list
    """
}

lncfinder_model = "${baseDir}/Models/LncFinder/${params.org_name}.model.RDS"
lncfinder_frequencies = "${baseDir}/Models/LncFinder/${params.org_name}.freq.RDS"

models_exist = file(lncfinder_model).exists() && file(lncfinder_frequencies).exists()

if (models_exist) {
    log.info "LncFinder models are exists for species, Mode: Inference Only."
    ch_gtf_for_mRNAs = Channel.empty()
    ch_model_lncfinder_final = Channel.fromPath(lncfinder_model)
    ch_frequencies_final = Channel.fromPath(lncfinder_frequencies)
    known_lncRNA_ch = Channel.empty()
    } else {
    log.info "LncFinder models not exists for species, Mode: Training + Inference."
    ch_gtf_for_mRNAs = Channel.fromPath(params.annotation)
    ch_model_lncfinder_final = Channel.empty()
    ch_frequencies_final = Channel.empty()
    known_lncRNA_ch = Channel.fromPath(params.known_lncRNAs_FA, checkIfExists: true)
}

process get_spliced_mRNAs {
    tag "Extracting spliced mRNA sequences"
    publishDir "LncFinder", pattern: '*.mRNAs.fa', mode: 'copy', overwrite: true
    cpus 4

    input: 
    path(annotation_GTF)
    path(genome_FA)
    
    output:
    path("${params.org_name}.mRNAs.fa"), emit: mRNAs
    
    script: 
    """
    gffread ${annotation_GTF} -g ${genome_FA} -x ${params.org_name}.mRNAs.fa
    """
}


process lncfinder_train {
    tag "Buiding LncFinder model and fequencies"
    publishDir "LncFinder", pattern: '*.{model,freq}.RDS', mode: 'copy', overwrite: true
    cpus cpu_free
    
    input: 
    path(mRNA_sequeces)
    path(known_lnc_fa)
    
    output: 
    path("${params.org_name}.model.RDS"), emit: model_lncfinder_ch
    path("${params.org_name}.freq.RDS"), emit: freq_lncfinder_ch
    
    script: 
    """
    Rscript $baseDir/bin/Train_LncFinder.R ${task.cpus} ${mRNA_sequeces} ${known_lnc_fa} ${params.org_name}.model.RDS ${params.org_name}.freq.RDS
    """
}


process run_lncfinder {
    tag "Running LncFinder"
    publishDir "LncFinder", pattern: '*.lncfinder.TSV', mode: 'copy', overwrite: true
    publishDir "LncFinder", pattern: 'final_*.lncfinder.list', mode: 'copy', overwrite: true

    cpus cpu_free
    
    input:
    path(lncfinder_model)
    path(lncfinder_frequencies)
    path(putative_transcripts_fa)
    
    output:
    path("Putative.lnc_NPCTs.lncfinder.TSV"), emit: lncfinder_output_ch
    path("final_lncRNAs.lncfinder.list"), emit: lncfinder_lncrnas_ch
    path("final_NPCTs.lncfinder.list"), emit: lncfinder_npcts_ch
    
    script:
    """
    Rscript $baseDir/bin/Run_LncFinder.R ${task.cpus} ${lncfinder_model} ${lncfinder_frequencies} ${putative_transcripts_fa} Putative.lnc_NPCTs.lncfinder.TSV final_lncRNAs.lncfinder.list final_NPCTs.lncfinder.list
    """
}

process download_pfam {
    tag "Downloading PFam database"
    publishDir "PfamScan", mode: 'copy', pattern: 'Pfam-A.*'
    publishDir "PfamScan", mode: 'copy', pattern: '*lncRNA-NPCT.pep'
    cpus cpu_free
    
    input:
    path(lncRNA_npcts_fa)
    
    output:
    path("Pfam-A.hmm.dat") , emit: pfam_dat_ch
    path("Pfam-A.hmm"), emit: pfam_ch
    path("Putative_lncRNA-NPCT.pep"), emit: pfam_pep_ch
    
    script:
    """
    wget -c ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz -O Pfam-A.hmm.dat.gz
    wget -c ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -O Pfam-A.hmm.gz
    gunzip -f Pfam-A.hmm.dat.gz > Pfam-A.hmm.dat
    gunzip -f Pfam-A.hmm.gz > Pfam-A.hmm 
    transeq -sequence ${lncRNA_npcts_fa} -outseq Putative_lncRNA-NPCT.pep
    """
}

process perform_pfamscan {
    tag "Runnging PfamScan"
    publishDir "PfamScan", mode: 'copy', pattern: '*.pfamscan.txt'
    publishDir "PfamScan", mode: 'copy', pattern: 'final_*.pfamscan.list'
    cpus cpu_free 
    
    input:
    path(pfam)
    path(lncRNA_npcts_pep)
    path(lncRNA_npcts_list)

    output:
    path("Putative.lnc_NPCTs.pfamscan.txt"), emit: pfamscan_out_ch
    path("pfam.log"), optional: true
    path("final_NPCTs.pfamscan.list"), emit: pfamscan_npcts_ch
    path("final_lncRNAs.pfamscan.list"), emit: pfamscan_lncrnas_ch

    script:
    """
    hmmpress -f ${pfam}
    hmmscan --cpu ${task.cpus}  -E 10e-5 --tblout Putative.lnc_NPCTs.pfamscan.txt ${pfam} ${lncRNA_npcts_pep} > pfam.log
    python $baseDir/bin/parse_pfam_results.py Putative.lnc_NPCTs.pfamscan.txt final_NPCTs.pfamscan.list
    Rscript $baseDir/bin/get_pfamscan_lncRNAs.R final_NPCTs.pfamscan.list ${lncRNA_npcts_list} final_lncRNAs.pfamscan.list
    """
}

process FEELnc_shuffle {
    tag "Running FEELnc shuffle mode"
    publishDir "FEELnc", pattern: 'FEELnc_shuffle.*', mode: 'copy', overwrite: true
    cpus cpu_free
    
    input:
    path(merged_GTF)
    path(annotation_GTF)
    path(genome_fa)

    output:
    path("FEELnc_shuffle_filter.log"), optional: true
    path("FEELnc_shuffle_filter.gtf"), emit: shuffle_filter_ch
    path("FEELnc_shuffle_codpot.log"), optional: true
    path("FEELnc_shuffle.codpot.lncRNA.gtf"), emit: shuffle_codpot_ch
    path("FEELnc_shuffle_classifier.log"), optional: true
    

    script:
    def out_prefix="FEELnc_shuffle.codpot"
    """
    conda run -n FEELnc FEELnc_filter.pl --proc=${task.cpus} --infile=${merged_GTF} --mRNAfile=${annotation_GTF} --biotype=transcript_biotype=protein_coding --monoex=-1 --outlog=FEELnc_shuffle_filter.log > FEELnc_shuffle_filter.gtf
    conda run -n FEELnc FEELnc_codpot.pl --infile=FEELnc_shuffle_filter.gtf --mRNAfile=${annotation_GTF} --genome=${genome_fa} --mode=shuffle --outname=${out_prefix} --outdir=. > FEELnc_shuffle_codpot.log
    conda run -n FEELnc FEELnc_classifier.pl -i ${out_prefix}.lncRNA.gtf -a ${annotation_GTF} --log=FEELnc_shuffle_classifier.log > FEELnc_shuffle.classifier.txt
    """
}

process FEELnc_intergenic {
    tag "Running FEELnc intergenic mode"
    publishDir "FEELnc", pattern: 'FEELnc_intergenic.*', mode: 'copy', overwrite: true
    cpus cpu_free
    
    input:
    path(merged_GTF)
    path(annotation_GTF)
    path(genome_fa)

    output:
    path("FEELnc_intergenic_filter.log"), optional: true
    path("FEELnc_intergenic_filter.gtf"), emit: intergenic_filter_ch
    path("FEELnc_intergenic_codpot.log"), optional: true
    path("FEELnc_intergenic.codpot.lncRNA.gtf"), emit: intergenic_codpot_ch
    path("FEELnc_intergenic_classifier.log"), optional: true
    

    script:
    def out_prefix="FEELnc_intergenic.codpot"
    """
    conda run -n FEELnc FEELnc_filter.pl --proc=${task.cpus} --infile=${merged_GTF} --mRNAfile=${annotation_GTF} --monoex=0 --size=200 --outlog=FEELnc_intergenic_filter.log > FEELnc_intergenic_filter.gtf
    conda run -n FEELnc FEELnc_codpot.pl --infile=FEELnc_intergenic_filter.gtf --mRNAfile=${annotation_GTF} --genome=${genome_fa} --mode=intergenic --outname=${out_prefix} --outdir=. > FEELnc_intergenic_codpot.log
    conda run -n FEELnc FEELnc_classifier.pl -i ${out_prefix}.lncRNA.gtf -a ${annotation_GTF} --log=FEELnc_intergenic_classifier.log > FEELnc_intergenic.classifier.txt
    """
}

process FEELnc_results {
    tag "Getting FEELnc results"
    publishDir "FEELnc", pattern: 'final_lncRNAs.feelnc.list', mode: 'copy', overwrite: true
    cpus cpu_free
    
    input:
    path(shuffle_GTF)
    path(intergenic_GTF)
    
    output:
    path("final_lncRNAs.feelnc.list"), emit: feelnc_output_ch
    
    script:
    """
    cat ${shuffle_GTF} ${intergenic_GTF} > temp.gtf
    python $baseDir/bin/parse_FEELnc_outputs.py temp.gtf final_lncRNAs.feelnc.list && rm -rf temp.gtf
    """
}

process lncRNA_outputs {
    tag "Getting lncRNA outputs"
    publishDir "LncRAnalyzer-summary", pattern: 'LncRAnalyzer-Lncs-*', mode: 'copy', overwrite: true
    cpus cpu_free
    
    input:
    path(feelnc_lncRNAs)
    path(cpat_lncRNAs)
    path(cpc2_lncRNAs)
    path(rnasamba_lncrnas)
    path(lgc_lncrnas)
    path(pfamscan_lncrnas)
    path(lncfinder_lncrnas)
    path(merged_GTF)
    
    output:
    path("LncRAnalyzer-Lncs-intersect.txt"), emit: final_lncrnas_ch
    path("LncRAnalyzer-Lncs-intersect.gtf"), emit: final_lncrnas_gtf_ch
    path("LncRAnalyzer-Lncs-predictions.TSV"), emit: lncrnas_pred_tsv_ch
    path("LncRAnalyzer-Lncs-predictions.tiff"), optional: true

    script:
    """
    Rscript $baseDir/bin/Lnc-Intersect.R ${feelnc_lncRNAs} ${cpat_lncRNAs} ${cpc2_lncRNAs} ${rnasamba_lncrnas} ${lgc_lncrnas} ${pfamscan_lncrnas} ${lncfinder_lncrnas} LncRAnalyzer-Lncs-intersect.txt
    python $baseDir/bin/subset_gtf.py ${merged_GTF} LncRAnalyzer-Lncs-intersect.txt LncRAnalyzer-Lncs-intersect.gtf
    Rscript $baseDir/bin/Lnc-Upset.R ${cpat_lncRNAs} ${cpc2_lncRNAs} ${rnasamba_lncrnas} ${feelnc_lncRNAs} ${lgc_lncrnas} ${pfamscan_lncrnas} ${lncfinder_lncrnas} LncRAnalyzer-Lncs-predictions.TSV LncRAnalyzer-Lncs-predictions.tiff
    """
}

process NPCTs_outputs {
    tag "Getting NPCTs outputs"
    publishDir "LncRAnalyzer-summary", pattern: 'LncRAnalyzer-NPCTs-*', mode: 'copy', overwrite: true
    cpus cpu_free
    
    input:
    path(cpat_npcts)
    path(cpc2_npcts)
    path(rnasamba_npcts)
    path(lgc_npcts)
    path(pfamscan_npcts)
    path(lncfinder_npcts)
    path(merged_GTF)
    path(genome_FA)
    
    output:
    path("LncRAnalyzer-NPCTs-intersect.txt"), emit: final_npcts_ch
    path("LncRAnalyzer-NPCTs-intersect.gtf"), emit: final_npcts_gtf_ch
    path("LncRAnalyzer-NPCTs-intersect.fa"), emit: final_npcts_fa_ch
    path("LncRAnalyzer-NPCTs-predictions.TSV"), emit: npcts_pred_tsv_ch
    path("LncRAnalyzer-NPCTs-predictions.tiff"), optional: true
    
    script:
    """
    Rscript $baseDir/bin/NPCTs-Intersect.R ${cpat_npcts} ${cpc2_npcts} ${rnasamba_npcts} ${lgc_npcts} ${pfamscan_npcts} ${lncfinder_npcts} LncRAnalyzer-NPCTs-intersect.txt
    python $baseDir/bin/subset_gtf.py ${merged_GTF} LncRAnalyzer-NPCTs-intersect.txt LncRAnalyzer-NPCTs-intersect.gtf
    gffread LncRAnalyzer-NPCTs-intersect.gtf -g ${genome_FA} -w LncRAnalyzer-NPCTs-intersect.fa
    Rscript $baseDir/bin/NPCTs-Upset.R ${cpat_npcts} ${cpc2_npcts} ${rnasamba_npcts} ${lgc_npcts} ${pfamscan_npcts} ${lncfinder_npcts} LncRAnalyzer-NPCTs-predictions.TSV LncRAnalyzer-NPCTs-predictions.tiff
    """
}

process lncRNA_classes {
    tag "Getting lncRNA classes"
    publishDir "LncRAnalyzer-summary", pattern: 'Summary_classification.*', mode: 'copy', overwrite: true
    publishDir "LncRAnalyzer-summary", pattern: 'LncRNA_classes.*', mode: 'copy', overwrite: true
    cpus cpu_free
    
    input:
    path(lncRNA_GTF)
    path(annotation_GTF)
    
    output:
    path("classification.log"), optional: true
    path("lncRNA_classes.txt")
    path("LncRNA_classes.TSV"), emit: lncrna_classes_ch
    path("Summary_classification.TSV"), emit: summary_classification_tsv_ch
    path("Summary_classification.tiff"), optional: true
    
    script:
    """
    conda run -n FEELnc FEELnc_classifier.pl -i ${lncRNA_GTF} -a ${annotation_GTF} --log=classification.log > lncRNA_classes.txt
    Rscript $baseDir/bin/get_lncRNA_classes.R lncRNA_classes.txt LncRNA_classes.TSV Summary_classification.TSV Summary_classification.tiff
    """
}

process get_lncRNA_bed {
    tag "Converting lncRNA GTF to bed"
    publishDir "TE_derived_lncRNAs", pattern: 'LncRAnalyzer-Lncs-*.bed', mode: 'copy', overwrite: true
    cpus cpu_free

    input:
    path(lncRNA_GTF)

    output:
    path("LncRAnalyzer-Lncs-intersect.bed"), emit: final_lncRNA_bed

    script:
    """
    gffread ${lncRNA_GTF} --bed -o temp.bed
    cut -f1-6 temp.bed > LncRAnalyzer-Lncs-intersect.bed
    rm temp.bed
    """
}

process get_TE_annotations {
    tag "Getting TE annotations"
    publishDir "TE_derived_lncRNAs", pattern: '*.{bed,gff3}', mode: 'copy', overwrite: true
    cpus cpu_free

    input:
    val(org_name)

    output:
    path("${org_name}_*.gff3"), emit: te_gffs_ch
    path("${org_name}_LTR.bed"), emit: ltr_bed_ch
    path("${org_name}_LINE.bed"), emit: line_bed_ch
    path("${org_name}_SINE.bed"), emit: sine_bed_ch
    path("${org_name}_MITE.bed"), emit: mite_bed_ch
    path("${org_name}_TIR.bed"), emit: tir_bed_ch
    path("${org_name}_Helitron.bed"), emit: helitron_bed_ch

    script:
    """
    python3 $baseDir/bin/download_TE.py ${org_name} LTR ${org_name}_LTR.gff3
    python3 $baseDir/bin/download_TE.py ${org_name} LINE ${org_name}_LINE.gff3
    python3 $baseDir/bin/download_TE.py ${org_name} SINE ${org_name}_SINE.gff3
    python3 $baseDir/bin/download_TE.py ${org_name} MITE ${org_name}_MITE.gff3
    python3 $baseDir/bin/download_TE.py ${org_name} TIR ${org_name}_TIR.gff3
    python3 $baseDir/bin/download_TE.py ${org_name} Helitron ${org_name}_Helitron.gff3
    python3 $baseDir/bin/format_TE_gff3tobed.py ${org_name}_LTR.gff3 ${org_name}_LTR.bed
    python3 $baseDir/bin/format_TE_gff3tobed.py ${org_name}_LINE.gff3 ${org_name}_LINE.bed
    python3 $baseDir/bin/format_TE_gff3tobed.py ${org_name}_SINE.gff3 ${org_name}_SINE.bed
    python3 $baseDir/bin/format_TE_gff3tobed.py ${org_name}_MITE.gff3 ${org_name}_MITE.bed
    python3 $baseDir/bin/format_TE_gff3tobed.py ${org_name}_TIR.gff3 ${org_name}_TIR.bed
    python3 $baseDir/bin/format_TE_gff3tobed.py ${org_name}_Helitron.gff3 ${org_name}_Helitron.bed
    """
}


process get_TE_derived_lncRNAs {
    tag "Getting TE-derived lncRNAs"
    publishDir "TE_derived_lncRNAs", pattern: '*.TSV', mode: 'copy', overwrite: true
    cpus cpu_free
    
    input:
    val(org_name)
    path(lnc_bed)
    path(ltr_bed)
    path(line_bed)
    path(sine_bed)
    path(mite_bed)
    path(tir_bed)
    path(helitron_bed)
    
    output:
    path("${org_name}_LTR.TSV"), emit: ltr_out_ch
    path("${org_name}_LINE.TSV"), emit: line_out_ch
    path("${org_name}_SINE.TSV"), emit: sine_out_ch
    path("${org_name}_MITE.TSV"), emit: mite_out_ch
    path("${org_name}_TIR.TSV"), emit: tir_out_ch
    path("${org_name}_Helitron.TSV"), emit: helitron_out_ch
    
    script:
    """
    bedtools intersect -wo -a ${lnc_bed} -b ${ltr_bed} | cut -f 1-12 > ${org_name}_LTR.TSV
    bedtools intersect -wo -a ${lnc_bed} -b ${line_bed} | cut -f 1-12 > ${org_name}_LINE.TSV
    bedtools intersect -wo -a ${lnc_bed} -b ${sine_bed} | cut -f 1-12 > ${org_name}_SINE.TSV
    bedtools intersect -wo -a ${lnc_bed} -b ${mite_bed} | cut -f 1-12 > ${org_name}_MITE.TSV
    bedtools intersect -wo -a ${lnc_bed} -b ${tir_bed} | cut -f 1-12 > ${org_name}_TIR.TSV
    bedtools intersect -wo -a ${lnc_bed} -b ${helitron_bed} | cut -f 1-12 > ${org_name}_Helitron.TSV
    """
}

process summary_TE_derived_lncRNAs {
    tag "Summary TE-derived lncRNAs"
    publishDir "LncRAnalyzer-summary", pattern: 'TE_derived_lncRNAs*', mode: 'copy', overwrite: true
    cpus cpu_free
    
    input:
    path(ltr_out)
    path(line_out)
    path(sine_out)
    path(mite_out)
    path(tir_out)
    path(helitron_out)

    output:
    path("TE_derived_lncRNAs.TSV")
    path("TE_derived_lncRNAs_summary.TSV")
    path("TE_derived_lncRNAs_summary.tiff"), optional: true

    script:
    """
    Rscript $baseDir/bin/get_TE_derived_lncRNAs.R ${ltr_out} ${line_out} ${sine_out} ${mite_out} ${tir_out} ${helitron_out} TE_derived_lncRNAs.TSV TE_derived_lncRNAs_summary.TSV TE_derived_lncRNAs_summary.tiff
    """
}

rel_annotation_ch = Channel.fromPath(params.annotation_related_species, checkIfExists: true)
rel_genome_ch =  Channel.fromPath(params.genome_related_species, checkIfExists: true)

process get_bed_for_slncky {
    tag "Getting bed files for Slncky"
    publishDir "Slncky", mode: 'copy', pattern: '*.bed'
    cpus cpu_free
    
    input:
    val(species)
    val(rel_species)
    path(annotation_GTF)
    path(rel_annotation_GTF)
    path(lnc_npcts_GTF)
    path(lncRNA_GTF)

    output:
    path("${species}.protein_coding.bed"), emit: coding_ch
    path("${rel_species}.protein_coding.bed"), emit: rel_coding_ch
    path("Putative_lncRNA-NPCT.bed"), emit: lnc_npcts_bed_ch
    path("LncRAnalyzer-Lncs-intersect.bed"), emit: final_lncrnas_bed_ch
    
    script:
    """
    gffread ${annotation_GTF} --bed -o temp.gtf
    cut -f 1-12 temp.gtf > ${species}.protein_coding.bed && rm temp.gtf
    gffread ${rel_annotation_GTF} --bed -o temp.gtf
    cut -f 1-12 temp.gtf > ${rel_species}.protein_coding.bed && rm temp.gtf
    gffread ${lnc_npcts_GTF} --bed -o temp.gtf
    cut -f 1-12 temp.gtf > Putative_lncRNA-NPCT.bed && rm temp.gtf
    gffread ${lncRNA_GTF} --bed -o temp.gtf
    cut -f 1-12 temp.gtf > LncRAnalyzer-Lncs-intersect.bed && rm temp.gtf
    """
}

liftover_file = params.liftover ? file(params.liftover).exists() : false

if(liftover_file) {
    liftover_ch = Channel.fromPath(params.liftover, checkIfExists: true)
    noncoding_ch = Channel.fromPath(params.noncoding, checkIfExists: true)
    mir_ch = Channel.fromPath(params.mir, checkIfExists: true)
    sno_ch = Channel.fromPath(params.sno, checkIfExists: true)
    rel_noncoding_ch = Channel.fromPath(params.rel_noncoding, checkIfExists: true)
    rel_mir_ch = Channel.fromPath(params.rel_mir, checkIfExists: true)
    rel_sno_ch = Channel.fromPath(params.rel_sno, checkIfExists: true)
    } else {
    liftover_ch = Channel.empty()
    noncoding_ch = Channel.empty()
    mir_ch = Channel.empty()
    sno_ch = Channel.empty()
    rel_noncoding_ch = Channel.empty()
    rel_mir_ch = Channel.empty()
    rel_sno_ch = Channel.empty()
}

process Slncky_config {
    tag "Getting annotation config for Slncky"
    publishDir "Slncky", mode: 'copy', pattern: '*.config'
    cpus cpu_free
    
    input:
    val(species)
    val(rel_species)
    path(genome_FA)
    path(coding)
    path(liftover)
    path(noncoding)
    path(mir)
    path(sno)
    path(rel_genome_FA)
    path(rel_coding)
    path(rel_noncoding)
    path(rel_mir)
    path(rel_sno)
    
    output:
    path("annotation.config"), emit: annotation_config_ch
    path("*.fai"), emit: genome_index_ch

    script:
    """
    samtools faidx ${genome_FA}
    samtools faidx ${rel_genome_FA}
    echo ">${species}" > annotation.config
    echo "CODING=${coding}" >> annotation.config
    echo "GENOME_FA=${genome_FA}" >> annotation.config
    echo "ORTHOLOG=${rel_species}" >> annotation.config
    if [ "${liftover_file}" == "true" ]; then
        echo "LIFTOVER=${liftover}" >> annotation.config
        echo "NONCODING=${noncoding}" >> annotation.config
        echo "MIRNA=${mir}" >> annotation.config
        echo "SNORNA=${sno}" >> annotation.config
    fi
    echo ">${rel_species}" >> annotation.config
    echo "CODING=${rel_coding}" >> annotation.config
    echo "GENOME_FA=${rel_genome_FA}" >> annotation.config
    echo "ORTHOLOG=${species}" >> annotation.config
    
    if [ "${liftover_file}" == "true" ]; then
        echo "NONCODING=${rel_noncoding}" >> annotation.config
        echo "MIRNA=${rel_mir}" >> annotation.config
        echo "SNORNA=${rel_sno}" >> annotation.config
    fi
    """    
}

process Slncky_run {
    tag "Running Slncky"
    publishDir "Slncky", mode: 'copy', pattern: 'slncky_out.*.txt'
    publishDir "Slncky", mode: 'copy', pattern: 'slncky_out.*.bed'
    cpus cpu_free
    
    input:
    val(species)
    path(annotation_config)
    path(lnc_npcts_bed)
    path(indexes)
    path(genome_FA)
    path(coding)
    path(liftover)
    path(noncoding)
    path(mir)
    path(sno)
    path(rel_genome_FA)
    path(rel_coding)
    path(rel_noncoding)
    path(rel_mir)
    path(rel_sno)

    output:
    path("slncky_out.*")
    
    script:
    def out_prefix = "slncky_out"
    """
    conda run -n cpc2-cpat-slncky slncky.v1.0 -n ${task.cpus} -c ${annotation_config} ${lnc_npcts_bed} ${species} ${out_prefix}
    """

}

process ortholog_search {
    tag "Running ortholog search using Slncky"
    publishDir "Slncky", mode: 'copy', pattern: '${rel_sp_name}*.txt'
    publishDir "Slncky", mode: 'copy', pattern: '${rel_sp_name}.*.bed'
    cpus cpu_free
    
    input:
    val(species)
    val(rel_sp_name)
    path(annotation_config)
    path(final_lncRNA_bed)
    path(indexes)
    path(genome_FA)
    path(coding)
    path(liftover)
    path(noncoding)
    path(mir)
    path(sno)
    path(rel_genome_FA)
    path(rel_coding)
    path(rel_noncoding)
    path(rel_mir)
    path(rel_sno)

    output:
    path("${rel_sp_name}.*")
    
    script:
    """
    conda run -n cpc2-cpat-slncky slncky.v1.0 -n ${task.cpus} -c ${annotation_config} --no_filter --minMatch=0.01 --no_orf --pad=100000 ${final_lncRNA_bed} ${species} ${rel_sp_name}
    """

}

process get_counts {
    tag "Getting lncRNA and PCG counts"
    publishDir "Counts", mode: 'copy', pattern: '*.{counts,TSV}'
    cpus cpu_free

    input:
    val(org_name)
    path(annotation_GTF)
    path(lncRNA_GTF)
    path(bam_list)

    output:
    path("${org_name}_PCG.counts"), emit: pcg_counts_ch
    path("${org_name}_lncRNA.counts"), emit: lnc_counts_ch
    path("${org_name}_PCG.TSV"), emit: pcg_tsv_ch
    path("${org_name}_lncRNA.TSV"), emit: lnc_tsv_ch

    script:
    if (isPaired) {
    """
    featureCounts -T ${task.cpus} -p --primary -t exon -g gene_id -F GTF -a ${annotation_GTF} -o ${org_name}_PCG.counts ${bam_list}
    featureCounts -T ${task.cpus} -p --primary -t exon -g transcript_id -F GTF -a ${lncRNA_GTF} -o ${org_name}_lncRNA.counts ${bam_list}
    sed 1,1d ${org_name}_PCG.counts | sed 's/[^[:space:]]*\\///g; s/\\.bam//g' > ${org_name}_PCG.TSV
    sed 1,1d ${org_name}_lncRNA.counts | sed 's/[^[:space:]]*\\///g; s/\\.bam//g' > ${org_name}_lncRNA.TSV
    """
    } else {
    """
    featureCounts -T ${task.cpus} --primary -t exon -g gene_id -F GTF -a ${annotation_GTF} -o ${org_name}_PCG.counts ${bam_list}
    featureCounts -T ${task.cpus} --primary -t exon -g transcript_id -F GTF -a ${lncRNA_GTF} -o ${org_name}_lncRNA.counts ${bam_list}
    sed 1,1d ${org_name}_PCG.counts | sed 's/[^[:space:]]*\\///g; s/\\.bam//g' > ${org_name}_PCG.TSV
    sed 1,1d ${org_name}_lncRNA.counts | sed 's/[^[:space:]]*\\///g; s/\\.bam//g' > ${org_name}_lncRNA.TSV
    """
    }
}

design_file = params.design ? file(params.design).exists() : false
if(design_file) {
    design_ch = Channel.fromPath(params.design, checkIfExists: true)
    } else {
    design_ch = Channel.empty()
}

process perform_DESeq2 {
    tag "Performing DESeq2"
    publishDir "LncRAnalyzer-summary", mode: 'copy', pattern: '${org_name}_*.DESeq2.TSV'
    cpus cpu_free

    input:
    val(org_name)
    path(experimental_design)
    path(PCG_counts)
    path(lncRNA_counts)

    output:
    path("${org_name}_PCG.DESeq2.TSV"), emit: pcg_deseq2_ch
    path("${org_name}_lncRNA.DESeq2.TSV"), emit: lnc_deseq2_ch
    
    script:
    """
    Rscript $baseDir/bin/DESeq2.R ${experimental_design} ${PCG_counts} ${org_name}_PCG.DESeq2.TSV
    Rscript $baseDir/bin/DESeq2.R ${experimental_design} ${lncRNA_counts} ${org_name}_lncRNA.DESeq2.TSV
    """
}

workflow {
    run_fastp(read_pairs_ch)
    rRNA_index(rRNA_input_ch)
    rRNA_unmapped_reads(run_fastp.out.trimmed_fastq_ch,rRNA_index.out.rRNA_index_ch.collect())
    extract_splicesites(annotation_input_ch)
    genome_index(genome_input_ch)
    reference_alignment(rRNA_unmapped_reads.out.clean_reads_ch,genome_index.out.genome_index_ch.collect(),extract_splicesites.out.known_splicesites_ch.collect())
    reference_guided_assembly(reference_alignment.out.aligned_bam_ch)
    merged_GTF(reference_guided_assembly.out.assembly_GTF_ch.map{it[1]}.collect())
    size_selected_GTF(merged_GTF.out.GTF_merged_ch)
    compare_annotations(size_selected_GTF.out.GTF_size_select_ch,annotation_input_ch)
    classcode_select(compare_annotations.out.annotation_compare_ch,genome_input_ch)
    perform_cpc2(classcode_select.out.lnc_npcts_FA_ch)
    perform_lgc(classcode_select.out.lnc_npcts_FA_ch)
    extract_cds(ch_gtf_for_cds, genome_input_ch)
    build_hexamer_table(extract_cds.out.cds, known_lncRNA_ch)
    build_logit_model(build_hexamer_table.out, extract_cds.out.cds, known_lncRNA_ch)
    final_hexamer = ch_hexamer_final.mix(build_hexamer_table.out)
    final_logit = ch_logit_final.mix(build_logit_model.out.model)
    final_feature = ch_feature_for_cutoff.mix(build_logit_model.out.feature)
    run_CPAT(final_hexamer, final_logit, classcode_select.out.lnc_npcts_FA_ch)
    get_cutoff(final_feature)
    get_cpat_results(run_CPAT.out.cpat_output_ch, get_cutoff.out.cpat_cutoff_ch)
    extract_mRNAs(ch_gtf_for_mRNAs, genome_input_ch)
    rnasamba_train(extract_mRNAs.out.mRNAs, known_lncRNA_ch)
    final_model_rnasamba = ch_model_rnasamba_final.mix(rnasamba_train.out.model_rnasamba_ch)
    rnasamba_classify(final_model_rnasamba, classcode_select.out.lnc_npcts_FA_ch)
    get_rnasamba_results(rnasamba_classify.out.rnasamba_output_ch)
    get_spliced_mRNAs(ch_gtf_for_mRNAs, genome_input_ch)
    lncfinder_train(get_spliced_mRNAs.out.mRNAs, known_lncRNA_ch)
    final_model_lncfinder = ch_model_lncfinder_final.mix(lncfinder_train.out.model_lncfinder_ch)
    final_frequencies = ch_frequencies_final.mix(lncfinder_train.out.freq_lncfinder_ch)
    run_lncfinder(final_model_lncfinder, final_frequencies , classcode_select.out.lnc_npcts_FA_ch)
    download_pfam(classcode_select.out.lnc_npcts_FA_ch)
    perform_pfamscan(download_pfam.out.pfam_ch, download_pfam.out.pfam_pep_ch, classcode_select.out.lnc_npcts_list_ch)
    FEELnc_shuffle(merged_GTF.out.GTF_merged_ch,annotation_input_ch,genome_input_ch)
    FEELnc_intergenic(merged_GTF.out.GTF_merged_ch,annotation_input_ch,genome_input_ch)
    FEELnc_results(FEELnc_shuffle.out.shuffle_codpot_ch, FEELnc_intergenic.out.intergenic_codpot_ch)
    lncRNA_outputs(FEELnc_results.out.feelnc_output_ch, get_cpat_results.out.cpat_lncrnas_ch, perform_cpc2.out.cpc2_lncrnas_ch, get_rnasamba_results.out.rnasamba_lncrnas_ch, perform_lgc.out.lgc_lncrnas_ch, perform_pfamscan.out.pfamscan_lncrnas_ch, run_lncfinder.out.lncfinder_lncrnas_ch, merged_GTF.out.GTF_merged_ch)
    NPCTs_outputs(get_cpat_results.out.cpat_npcts_ch, perform_cpc2.out.cpc2_npcts_ch, get_rnasamba_results.out.rnasamba_npcts_ch, perform_lgc.out.lgc_npcts_ch, perform_pfamscan.out.pfamscan_npcts_ch, run_lncfinder.out.lncfinder_npcts_ch, merged_GTF.out.GTF_merged_ch, genome_input_ch)
    lncRNA_classes(lncRNA_outputs.out.final_lncrnas_gtf_ch, annotation_input_ch)
    if (params.clade == "plants") { 
    get_lncRNA_bed(lncRNA_outputs.out.final_lncrnas_gtf_ch)
    get_TE_annotations(params.org_name)
    get_TE_derived_lncRNAs(params.org_name, get_lncRNA_bed.out.final_lncRNA_bed, get_TE_annotations.out.ltr_bed_ch, get_TE_annotations.out.line_bed_ch, get_TE_annotations.out.sine_bed_ch, get_TE_annotations.out.mite_bed_ch, get_TE_annotations.out.tir_bed_ch, get_TE_annotations.out.helitron_bed_ch)
    summary_TE_derived_lncRNAs(get_TE_derived_lncRNAs.out.ltr_out_ch, get_TE_derived_lncRNAs.out.line_out_ch, get_TE_derived_lncRNAs.out.sine_out_ch, get_TE_derived_lncRNAs.out.mite_out_ch, get_TE_derived_lncRNAs.out.tir_out_ch, get_TE_derived_lncRNAs.out.helitron_out_ch)
    }
    get_bed_for_slncky(params.org_name, params.rel_sp_name, annotation_input_ch, rel_annotation_ch, classcode_select.out.lnc_npcts_GTF_ch, lncRNA_outputs.out.final_lncrnas_gtf_ch)
    Slncky_config(params.org_name, params.rel_sp_name, genome_input_ch, get_bed_for_slncky.out.coding_ch, liftover_ch, noncoding_ch, mir_ch, sno_ch, rel_genome_ch, get_bed_for_slncky.out.rel_coding_ch, rel_noncoding_ch, rel_mir_ch, rel_sno_ch)
    Slncky_run(params.org_name, Slncky_config.out.annotation_config_ch, get_bed_for_slncky.out.lnc_npcts_bed_ch, Slncky_config.out.genome_index_ch.collect(), genome_input_ch, get_bed_for_slncky.out.coding_ch, liftover_ch, noncoding_ch, mir_ch, sno_ch, rel_genome_ch, get_bed_for_slncky.out.rel_coding_ch, rel_noncoding_ch, rel_mir_ch, rel_sno_ch)
    if(liftover_file) {
    ortholog_search(params.org_name, params.rel_sp_name, Slncky_config.out.annotation_config_ch, get_bed_for_slncky.out.final_lncrnas_bed_ch, Slncky_config.out.genome_index_ch.collect(), genome_input_ch, get_bed_for_slncky.out.coding_ch, liftover_ch, noncoding_ch, mir_ch, sno_ch, rel_genome_ch, get_bed_for_slncky.out.rel_coding_ch, rel_noncoding_ch, rel_mir_ch, rel_sno_ch)
    }
    get_counts(params.org_name, annotation_input_ch, lncRNA_outputs.out.final_lncrnas_gtf_ch, reference_alignment.out.aligned_bam_ch.map{it[1]}.collect())
    perform_DESeq2(params.org_name, design_ch, get_counts.out.pcg_tsv_ch, get_counts.out.lnc_tsv_ch)
}
