#configfile:"config.yaml"


"""
The purpose of this script is to convert a given vcf file to sequences.
GRCh38 is used as a reference genome. A bed file of loci of interest is
used to first get reference sequences. This is followed by bcftools
consensus to get variant sequences for both haplotype
"""
# INPUTS
reference_hg38="/eva/edatums/reference_materials/reference_genomes/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# a bed file for regions of interest
bedfile="hg38_liftedover_ucsc_correct.final.sorted.042023.dipBed.overlapSitesOnly.subsetFilt.bed"

# a vcf file to get the variants
vcf="HPRC-cur.20211005-align2-GRCh38.dip.vcf.gz"
#vcf="HPRC-cur.20211005-align2-GRCh38.dip.filtered.bcf"



#print(reference_hg38)

"""
TODOs: 
1. Bedtools intersect to find overlap between our loci and de novo assembled coord. or 
   any other bed file to know the regions in which the variants were called. It is a safe
   assumption that within that region is there was no variant then it is a homozygous loci. 
"""

rule all:
    input:
    #    expand("results/{sample}_vcf2fasta.fa", sample=bedfile)
        #expand("results/{bed}.txt",bed=bedfile),
        "results/GRCh38_samtools_faidx_regions.fa"

# intersect vcf file with bedfile to get variants only in the
# region of interest: use bcftools view -R
rule vcf_in_regions_of_interest:
    input:
        vcf
    output:
        expand("results/{INVCF}.filtered.bcf", INVCF=vcf)
    shell:
        "bcftools view -R {bedfile} -Ob -o {output} {input}"

# take the bedfile convert it to 1-based start. This is to convert it
# to an intervals file that will be used with samtools
rule bedfile_to_intervals:
    input:
        bedfile
    output:
        "results/{bed}.txt"
    shell:
        """
            cat {input} | awk '{{sum=$2+1}}{{print $1":"sum"-"$3}}' > {output}
        """

# samtools faidx to get just the regions of interest from GRCh38
rule:
    input:
        expand("results/{bed}.txt",bed=bedfile)
        #expand("results/{{bed}}.txt",bed=bedfile)
    output:
        "results/GRCh38_samtools_faidx_regions.fa"
    shell:
        """
        samtools faidx {reference_hg38} -r {input} -o {output}
        """