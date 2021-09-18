#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

params.bowtie_path = "${DB_REF}/bowtie/Z.capensis_transcriptome/"

include{bowtie2_build; bowtie2_SE} from './nf-lib/db_algos.nf'
include{fix_ReadName; infoseq; infoseq_stats} from './nf-lib/raw_reads.nf'
include{transcript_est_alignfree; salmon_index; salmon_quant} from './nf-lib/expression_analysis.nf'



workflow{

     SE_reads = Channel.fromPath("${params.WD}/SortmeRNA/*.fastq")
     denovo_ref = Channel.value("${params.WD}/Trinity/Trinity.fasta")
     sample_file = Channel.fromPath(params.sample_file)
     bt2_index_path =  Channel.value(params.bowtie_path)
     fix_ReadName(SE_reads)
     bowtie2_build(denovo_ref, bt2_index_path)
     bowtie2_SE(fix_ReadName.out, denovo_ref, bt2_index_path)    
     infoseq(SE_reads)
     infoseq_stats(infoseq.out.collect())
     transcript_est_alignfree(fix_ReadName.out.collect(), denovo_ref, sample_file, infoseq_stats.out.readlen_stats, "salmon")
     salmon_index(denovo_ref)    
     salmon_quant(bowtie2_SE.out.BAM, salmon_index.out, denovo_ref)
          

}
