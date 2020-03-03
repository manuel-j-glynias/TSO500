# TSO500
A project to annotate and determine the reporting status of variants from the TSO500 pipeline
This project has three separate executables:
  1) strands.py (not used in this analysis, but I wanted to have a record of it.)  It used the mygene.info REST API to get strand information for the genes in a data set.
  2) vep37.py which annotates each variant in a datafile, and saves the results to MongoDB-Atlas
  3) decision.py, which uses the mongoDB databases of variants, and executes the decision tree to determine report status, saving these results into a tsv file
  
  
  ##Annotation:
  1) Data are read from the TSO500_UniqueVariants_runs1-6.csv file in the data directory using the read_tso_unique_variants method in the TSO500.py file.
  2) A key is created for the variant, chr_pos_ref/alt.  If the mongoDB annvar collection contains a variant with this key, it is fetched
  3) If not, the variant is written into VEP input format (using get_vep_var_as_string) then sent to the VEP annotation service in batches of 100.
  4) The results contain arrays of transcript_consequences for each variant.  An attempt is made to find the consequence that most closely matches the variant decision already made by the OmniSeq pipeline using the method get_best_t_c
  5) The data from this transcript_consequence are added to the variant object, as is additional information from CKB about the gene and variant, and from ClinVar and HotSpot about the variant.
  6) The data are stored as a json object into mongoDB atlas using the key.
  
  For example, EGFR L858R is:
  ```"EGFR","L858R","c.2573T>G","21","7","55259515","T","G","Substitution - Missense","NM_005228.3","NP_005219.2"```
  
  in the TSO500_UniqueVariants_runs1-6.csv file.  Using the chromosome, pos, ref and alt gives the key:
  7_55259515_T/G
  which after searching in mongoDB omni.acc-var gives:
```
  {
  "HGNC_Symbol": "EGFR",
  "omni_c_dot": "c.2573T>G",
  "omni_p_dot": "L858R",
  "omni_MutationType": "Substitution - Missense",
  "chrom": "7",
  "ref": "T",
  "alt": "G",
  "pos_hg19": "55259515",
  "type": "snv",
  "report_status": "not_reported",
  "gene": "EGFR",
  "gene_description": "EGFR (HER1), epidermal growth factor receptor, is a tyrosine kinase receptor, which activates RAS/RAF/MEK and PI3K/AKT/mTOR pathways, leading to increased cell proliferation and growth (PMID: 24312144). EGFR activating mutations, amplification, and overexpression are found in a variety of tumors, including non-small cell lung cancer (PMID: 26609494, PMID: 30284706) and colorectal cancer (PMID: 30243897), and the EGFRvIII variant is commonly found in glioblastoma (PMID: 30201736).",
  "gene_category": "Oncogene",
  "cdot": "c.2573T>G",
  "pdot": "L858R",
  "protein_start": 858,
  "full_name": "EGFR L858R",
  "mutation_type": "missense_variant",
  "is_protein_altering": true,
  "is_truncating_variants": false,
  "polyphen_prediction": "probably_damaging",
  "sift_prediction": "deleterious",
  "predicted_deleterious": true,
  "ckb_id": "273",
  "variant_description": "EGFR L858R lies within the protein kinase domain of the Egfr protein (UniProt.org). L858R results in increased kinase activity, is transforming in cell culture (PMID: 29533785, PMID: 28979142), demonstrates a growth advantage in IL-3-depleted cells in culture, as compared to wild-type Egfr (PMID: 30952700), and promotes tumor formation in mouse models (PMID: 16187797).",
  "protein_effect": "gain of function",
  "variant_type": "missense",
  "gDot": "chr7:g.55191822T>G",
  "cDot": "c.2573T>G",
  "cdot_pos": 2573,
  "pdot_pos": 858,
  "gene_id": 1956,
  "is_gain_of_function": true,
  "is_loss_of_function": false,
  "in_clinvar": true,
  "clinvar_significance": "drug response",
  "is_clinvar_not_benign": true,
  "is_clinvar_benign": false,
  "is_clinvar_pathogenic": true,
  "clinvar_explain": "drug response(6)/pathogenic(2)/likely pathogenic(1)",
  "variant_id": "16609",
  "hotspots": [
    {
      "begin": 858,
      "end": 858,
      "residue": "L858"
    }
  ],
  "is_near_GOF_LOF_mutation": true,
  "num_upstream_GOF_LOF_mutations": 4,
  "num_downstream_GOF_LOF_mutations": 6
}
```

## Decision
  The main method in decision.py reads each annotated variant from mongodb and determines is reporting status.  
  Using the method is_tso_snv_reportable in reportable_variants.py, each variant flows through the decision tree and has these results added to it.
  Results are written to the file decisions.tsv.  The columns in this file are:
 #### gene
 The HGNC identifier for the gene into which this variant falls, according to vep
#### cdot
The cdot for this variant, according to vep
#### pdot
The pdot for this variant, according to vep
#### gene_category
The whether this gene is an oncogene or a tumor suppressor
#### mutation_type
The mutation type (missense, frameshift, etc.) as determined by vep
#### report_status
Final result of the decision tree:  should this variant be reported or not
#### reasons
List of reasons why it should be reported
#### lack_of_reasons
List of reasons why it might not be reported(i.e. 'NO' steps in the decision tree)
#### is_protein_altering
Does this alteration change the protein sequence
#### in_clinvar
Is the variant in clinvar, using the pDot as the key into the Clinvar mongodb collection
#### clinvar_significance
whether clinvar calls the variant benign or pathogenic
#### is_gain_of_function
Does CKB call this variant GOF
#### is_loss_of_function
Does CKB call this variant LOF
#### hotspots
Does this varint fall into a hot spot region
#### predicted_deleterious
Do Sift and Polyphen agree that this variant is deleterious (from vep)
#### is_truncating_variants
Does this variant truncate the protein
####is_near_GOF_LOF_mutation
Is this variant in a region of annotated GOF or LOF variants
#### omni_gene
What was the gene call from the OmniSeq pipeline
#### omni_cdot
What was the cdot call from the OmniSeq pipeline
#### omni_pdot
What was the pdot call from the OmniSeq pipeline
