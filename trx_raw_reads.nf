#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.WD = "trX_assembly"
params.m_mem = 8
params.mtp_cores = 4

include{ ubam2fastq; fastqc_SE; multiqc; trimmomatic_SE} from './nf-lib/raw_reads.nf' 


process fix_ReadName{

        echo true
	input:
	     path SE_read

	output:
	     path SE_read
shell:
"""
     
     sed -i -e '1~4s/\$/\\\\1/g' read.fastq 
     

"""

}



workflow{

	uBAM = Channel.fromPath("/home/drewx/Documents/Zoster_capensis_TrX/*.bam")
	
	ubam2fastq(uBAM)

	fix_ReadName(ubam2fastq.out)
	
	fastqc_SE(fix_ReadName.out)

        multiqc(fastqc_SE.out)

        trimmomatic_SE(fix_ReadName.out)
	
}

