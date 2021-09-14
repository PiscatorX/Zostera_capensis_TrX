#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

params.WD = "/home/drewx/Documents/Zoster_capensis_TrX/Quantx"
params.bowtie_path = "${DB_REF}/bowtie/Z.capensis_transcriptome/"


include{infoseq; infoseq_stats} from './nf-lib/raw_reads.nf'
include{transcript_est_alignfree; salmon_index; salmon_quant} from './nf-lib/expression_analysis.nf'



// params.sample_file = "/opt/DB_REF/Metadata/Z_capensis.Transcriptome.tsv"




// include{bowtie2_build; bowtie2_SE} from './db_algos.nf'
// include{sam2bam} from './align_utils.nf'
// include{sortmerRNA_SE} from './db_algos.nf'

// include{infoseq; infoseq_stats; ubam2fastq; fastqc_SE; multiqc; trimmomatic_SE; fix_ReadName} from './raw_reads.nf'
// include{transcript_est_alignfree; salmon_index; salmon_quant} from './expression_analysis.nf'
// include{Trinity_SE} from './assembly.nf'
// include{merge_tsv} from  './misc_utils.nf'



// workflow get_fastq{

// 	take:
// 	    BAM_files

// 	main:
//         	ubam2fastq(BAM_files)
// 		fastqc_SE(ubam2fastq.out, "RawReads")
// 		multiqc(fastqc_SE.out.collect(), "RawReads")
// 	emit:
// 	    ubam2fastq.out
        	
	
// }




// workflow process_fastq{

//     take:
// 	 FastQ_reads

//      main:
// 	 trimmomatic_SE(FastQ_reads)
// 	 fastqc_SE(trimmomatic_SE.out.TrimmedRead, "Trimmomatic")
// 	 multiqc(fastqc_SE.out.collect(), "Trimmomatic")

//      emit:
//          trimmomatic_SE.out.TrimmedRead
	 
// }




 workflow{

     SE_reads = Channel.fromPath("${PWD}/trX_assembly/trimmomatic/*.fastq")
     denovo_ref = Channel.value("${PWD}/trX_assembly/Trinity/Trinity.fasta")
     sample_file = Channel.fromPath(params.sample_file)
     infoseq(SE_reads)
     infoseq_stats(infoseq.out.collect())
     transcript_est_alignfree(SE_reads.collect(), denovo_ref, sample_file, infoseq_stats.out.readlen_stats, "salmon")
//     sample_file = Channel.fromPath(params.sample_file)
//     bowtie2_SE(SE_reads, denovo_ref, )
//     bowtie2_build(denovo_ref, bt2_index_path)
//     bowtie2_SE(SE_reads, denovo_ref, bt2_index_path)    
//     salmon_index(denovo_ref)    
//     salmon_quant(bowtie2_SE.out.BAM, salmon_index.out, denovo_ref)

}




    
    
