#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.WD = "trX_assembly"
params.m_mem = 8
params.mtp_cores = 4

include{ ubam2fastq; fastqc_SE; multiqc; trimmomatic_SE} from './nf-lib/raw_reads.nf' 


workflow{

	uBAM = Channel.fromPath("/home/drewx/Documents/Zoster_capensis_TrX/*.bam")
	
	ubam2fastq(uBAM)
	
	fastqc_SE(ubam2fastq.out)

        multiqc(fastqc_SE.out)

        trimmomatic_SE(ubam2fastq.out)
	
}
