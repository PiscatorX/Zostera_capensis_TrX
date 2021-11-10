#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


denovo_ref = "${params.WD}/Trinity/Trinity.fasta"

swissprot = "${DB_REF}/SwissProt/uniprot_sprot.fasta"
swissprot_path = "${DB_REF}/BLAST/SwissProt/"
PfamA_hmm =  "${DB_REF}/Pfam/Pfam-A.hmm"
mol_type = "prot"
//e_value = "1e-5"
e_value = "100"
min_len = 10

include{makeBlastDB; blast_tophit; hmmscan} from './nf-lib/db_algos.nf'
include{transdecoder} from './nf-lib/prediction_algos.nf'



workflow{

      makeBlastDB(swissprot, mol_type, swissprot_path)
      blast_tophit(denovo_ref, swissprot, makeBlastDB.out.DB_path, e_value)
      transdecoder(denovo_ref, swissprot, swissprot_path,  PfamA_hmm, e_value,  min_len)
      hmmscan(transdecoder.out.pep, PfamA_hmm)

}


