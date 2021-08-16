#!/usr/bin/env nextflow
nextflow.enable.dsl=2



include{ubam2fastq; fastqc_SE; multiqc; trimmomatic_SE; fix_ReadName} from './nf-lib/raw_reads.nf'
include{sortmerRNA_SE} from './nf-lib/db_algos.nf'
include{Trinity_SE} from './nf-lib/assembly.nf'


workflow get_fastq{

	take:
	    BAM_files

	main:
        	ubam2fastq(BAM_files)
		fastqc_SE(ubam2fastq.out, "RawReads")
		multiqc(fastqc_SE.out.collect(), "RawReads")
	emit:
	    ubam2fastq.out
        	
	
}




workflow process_fastq{

    take:
	 FastQ_reads

     main:
	 trimmomatic_SE(FastQ_reads)
	 fastqc_SE(trimmomatic_SE.out.TrimmedRead, "Trimmomatic")
	 multiqc(fastqc_SE.out.collect(), "Trimmomatic")

     emit:
         trimmomatic_SE.out.TrimmedRead
}




workflow {


	BAM_files = Channel.fromPath(params.raw_bams, checkIfExists: true) 
        SortmeRNA_Reflist =  Channel.value(params.SortmeRNA_Reflist)
        sample_file = Channel.fromPath(params.sample_file)
	
	get_fastq(BAM_files)
	process_fastq(get_fastq.out)
	sortmerRNA_SE(process_fastq.out, SortmeRNA_Reflist)
	fix_ReadName(sortmerRNA_SE.out.SE_mRNA_read)
	fastqc_SE(fix_ReadName.out, "SortmeRNA")
	multiqc(fastqc_SE.out.collect(), "SortmeRNA")
	Trinity_SE(fix_ReadName.out.collect(), sample_file) 

}






    
    

