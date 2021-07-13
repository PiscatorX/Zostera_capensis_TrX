#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include{ubam2fastq} from './nf-lib/raw_reads.nf' 


workflow{

	uBAM = Channel.fromPath("/home/drewx/Documents/Zoster_capensis_TrX/*.bam")
	ubam2fastq(uBAM)

}