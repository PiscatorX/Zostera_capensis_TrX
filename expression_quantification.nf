#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

params.bowtie_path = "${DB_REF}/bowtie/Z.capensis_transcriptome/"

include{bowtie2_build; bowtie2_SE; bam_flagstat} from './nf-lib/db_algos.nf'
include{fix_ReadName; infoseq; infoseq_stats} from './nf-lib/raw_reads.nf'
include{transcript_est_alignfree; abundance_estimates_to_matrix; salmon_index; salmon_quant} from './nf-lib/expression_analysis.nf'



workflow{

     SE_reads = Channel.fromPath("${params.WD}/FixedReadNames/*.fastq")
     denovo_ref = Channel.value("${params.WD}/Trinity/Trinity.fasta")
     sample_file = Channel.fromPath(params.sample_file)
     
     bt2_index_path =  Channel.value(params.bowtie_path)
     bowtie2_build(denovo_ref, bt2_index_path)
     bowtie2_SE(SE_reads, denovo_ref, bt2_index_path)
     bam_flagstat(bowtie2_SE.out.BAM)

     infoseq(SE_reads)
     infoseq_stats(infoseq.out.collect())
     transcript_est_alignfree(SE_reads.collect(), denovo_ref, sample_file, infoseq_stats.out.readlen_stats, "salmon")
     salmon_index(denovo_ref)    
     salmon_quant(bowtie2_SE.out.BAM, salmon_index.out, denovo_ref)
     abundance_estimates_to_matrix(transcript_est_alignfree.out.quant_dir, denovo_ref, "quant.sf", "salmon")          

}
