#configfile:"config.yaml"


"""
The purpose of this script is to convert a given vcf file to sequences.
GRCh38 is used as a reference genome. A bed file of loci of interest is
used to first get reference sequences. This is followed by bcftools
consensus to get variant sequences for both haplotype
"""

"""
TODOs: 
1. Bedtools intersect to find overlap between our loci and de novo assembled coord. or 
   any other bed file to know the regions in which the variants were called. It is a safe
   assumption that within that region is there was no variant then it is a homozygous loci. 
"""
#######################################
# INPUTS
#######################################
reference_hg38="/eva/edatums/reference_materials/reference_genomes/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# a bed file for regions of interest
bedfile="hg38_liftedover_ucsc_correct.final.sorted.042023.dipBed.overlapSitesOnly.subsetFilt.bed"

# a vcf file to get the variants
vcf="HPRC-cur.20211005-align2-GRCh38.dip.vcf.gz"
#vcf="HPRC-cur.20211005-align2-GRCh38.dip.filtered.bcf"


#######################################
# Rules
#######################################
rule all:
    input:
    #   expand("results/{sample}_vcf2fasta.fa", sample=bedfile)
        #expand("results/{bed}.txt",bed=bedfile),
        #"results/GRCh38_samtools_faidx_regions.fa"
        expand("results/{INVCF}.bcftools.vcf2fasta.hap1.fa", INVCF=vcf),
        expand("results/{INVCF}.bcftools.vcf2fasta.hap2.fa", INVCF=vcf),
        expand ("results/{INVCF}.bcftools.vcf2fasta.hap1hap2.fa", INVCF=vcf),
        expand("results/{INVCF}.bcftools.vcf2fasta.hap1hap2.oneline.fa", INVCF=vcf),
        expand("results/{bed}.txt", bed=bedfile),
        expand("results/{INVCF}.filtered.bcf{ext}", INVCF=vcf, ext=["", ".csi"])


# intersect vcf file with bedfile to get variants only in the
# region of interest: use bcftools view -R
rule vcf_in_regions_of_interest:
    input:
        vcf
    output:
        VCFOUT = "results/{INVCF}.filtered.bcf",
        CSIOUT = "results/{INVCF}.filtered.bcf.csi"
    shell:
        "bcftools view -R {bedfile} -Ob -o {output.VCFOUT} {input} && "
        "bcftools index {output.VCFOUT}"

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
rule samtools_faidx:
    input:
        #expand("results/{bed}.txt",bed=bedfile)
        #expand("results/{{bed}}.txt",bed=bedfile)
        INTERVALS=expand("results/{bed}.txt",bed=bedfile),
        REF=reference_hg38
    output:
        "results/GRCh38_samtools_faidx_regions.fa"
    shell:
        """
        samtools faidx {input.REF} -r {input.INTERVALS} -o {output}
        """

# get both haplotypes from vcf to sequence
rule bcftools_consensus_hap1:
    input:
        fastafile = "results/GRCh38_samtools_faidx_regions.fa",
        inputvcf = "results/{INVCF}.filtered.bcf"
    output:
        "results/{INVCF}.bcftools.vcf2fasta.hap1.fa"
    shell:
        "bcftools consensus --haplotype 1pIu -f {input.fastafile} {input.inputvcf} > {output}"


rule bcftools_consensus_hap2:
    input:
        fastafile="results/GRCh38_samtools_faidx_regions.fa",
        inputvcf="results/{INVCF}.filtered.bcf"
    output:
        "results/{INVCF}.bcftools.vcf2fasta.hap2.fa"
    shell:
        "bcftools consensus --haplotype 2pIu -f {input.fastafile} {input.inputvcf} > {output}"


# combine both files into one haplotype file
rule combine_hap:
    input:
        hap1="results/{INVCF}.bcftools.vcf2fasta.hap1.fa",
        hap2="results/{INVCF}.bcftools.vcf2fasta.hap2.fa"
    output:
        "results/{INVCF}.bcftools.vcf2fasta.hap1hap2.fa"
    shell:
        "cat {input.hap1} {input.hap2} > {output}"


# optionally convert the multiline fasta file to one line fasta file
rule convert_multiline_to_one_line:
    input:
        "results/{INVCF}.bcftools.vcf2fasta.hap1hap2.fa"
    output:
        "results/{INVCF}.bcftools.vcf2fasta.hap1hap2.oneline.fa"
    shell:
        """ 
        perl -pe '$. >1 and /^>/ ? print "\n" : chomp' {input} > {output}
        """