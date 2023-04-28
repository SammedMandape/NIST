# Bioinformatics pipeline

The snakemake workflow described here generates sequences (for both haplotypes) given a bed 
file of loci and a VCF file. Optionally, you can also convert the multiline fasta into one line.
The image shows the entire bioinformatics pipeline with applications using targeted forensic markers.
However, the snakemake workflow can be broadly applied to any VCF file and loci of interest. 


```
snakemake -s vcf2seq_v2.smk -c32
```

For dry run, use

`snakemake -nps vcf2seq_v2.smk -c32`

To look at the summary of the snakemake outputs, use

`snakemake -s vc2seq_v2 -c32 --summary`