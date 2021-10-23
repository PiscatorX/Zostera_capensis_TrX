#!/usr/bin/env nextflow
nextflow.enable.dsl = 2



params.blastDB = "${DB_REF}/SwissProt/uniprot_sprot.fasta"
params.blastDB_path = "${DB_REF}/BLAST/SwissProt/"
params.blastDB_moltype = "prot"
params.genome_ref = "${DB_REF}/Genomes/GCA_001185155.1_Zosma_marina.v.2.1_genomic.fna"

include{infoseq; infoseq_stats} from'./raw_reads.nf'
include{assembly_stats; transrate; rnaQuast; analyze_blastTophits} from './assembly.nf'
include{makeBlastDB; blast_tophit; busco_auto_euk} from './db_algos.nf'
//include{rsem_eval;  transrate; rnaQuast} from './assembly.nf'
//include{bowtie2_build; bowtie2_SE; bam_index; samtools_fasta_ref} from './db_algos.nf'
//include{assemply_stats; analyze_blastTophits} from './assembly.nf'




workflow{

     SE_reads = Channel.fromPath("${params.WD}/FixedReadNames/*.fastq")
     denovo_ref = Channel.value("${params.WD}/Trinity/Trinity.fasta")

     blastDB = Channel.fromPath(params.blastDB)
     mol_type = Channel.value(params.blastDB_moltype)
     blastDB_path = Channel.fromPath(params.blastDB_path)
     genome_ref = Channel.value(params.genome_ref)
    
     infoseq(SE_reads)
     infoseq_stats(infoseq.out.collect())

    //assembly_stats(denovo_ref)
    //transrate(SE_reads.collect(), denovo_ref, genome_ref)
    //makeBlastDB(blastDB, mol_type, blastDB_path)
    //blast_tophit(denovo_ref, blastDB, blastDB_path)
    //blast_tophit.out.view()
    rsem_eval(SE_reads.collect(), denovo_ref, genome_ref, infoseq_stats.out.readlen_stats)
    //analyze_blastTophits(blast_tophit.out, denovo_ref, blastDB)
    
    // rnaQuast(SE_reads.collect(), denovo_ref, genome_ref, infoseq_stats.out.readlen_stats) 
    // busco_auto_euk(denovo_ref, species)

}


