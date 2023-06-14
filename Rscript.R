library(tidyverse)
library(data.table)
#library(spgs)
library(Biostrings)

# function to rev_complement a sequence

reverse_complement <- function (dna) {
  #browser()
  bases = strsplit(dna, '')[[1L]]
  compl = ifelse(bases == 'A', 'T',
                 ifelse(bases == 'C', 'G',
                        ifelse(bases == 'G', 'C',
                               ifelse(bases == 'T', 'A', 'N'))))
  
  paste(rev(compl), collapse = '')
}


#mydata <- fread("hg38_liftedover_ucsc_correct.final.sorted.dipBed.overlapSitesOnly.bed",
#                header = F, sep = "\t")

#mydata %>% 
#  mutate(reg=paste0(V1,":",V2,"-",V3))

# get hg19 data with loci names
mydata <- fread("hg19_data_with_lociNames.txt", sep = "\t", header = F)

# get hg38 lifted over data
mydata_lifted <- fread("hg38_liftedover_ucsc_correct_04202023.bed",
                       sep="\t", header=F)

c("chr","start-0-based","end-1-based","hg19-1-based-coord","junk") -> colnames(mydata_lifted)

# convert cols to integer to make all cols same data type for joins
mydata_lifted_1 <- mydata_lifted %>% separate(`hg19-1-based-coord`, 
                           into = c("chrHg19","hg19-start-1based","hg19-stop-1based")) %>%
                  mutate(`hg19-start-1based`=as.integer(`hg19-start-1based`),
                         `hg19-stop-1based`=as.integer(`hg19-stop-1based`))

# do left join and get loci names with lifted over data
mydata_lifted_1 %>% left_join(mydata, by=c("chrHg19"="V2",
                                           "hg19-start-1based"="V3",
                                           "hg19-stop-1based"="V4" )) %>% 
  mutate("start-1based"=`start-0-based`+1)-> mydata_lifted_2

#####################################
# get fasta files for both haplotypes
#####################################
myfasta_hap1 <- fread("hg38_liftedover_042023_dipBed.overlapSitesOnly_subsetFilt.bcftools_vcf2fasta.oneliner.fa",
                 header=F)
myfasta_hap2 <- fread("hg38_liftedover_042023_dipBed.overlapSitesOnly_subsetFilt.bcftools_vcf2fasta.hap2.oneliner.fa",
                      header=F)
bind_rows(myfasta_hap1, myfasta_hap2) -> myfasta
myfasta %>% 
  separate(V1, into = c(NA,"chr","start-1based","stop-1based","vcf2seq")) %>% 
  mutate(`start-1based`=as.integer(`start-1based`),
         `stop-1based`=as.integer(`stop-1based`))-> myfasta_1

#####################################
# get loci names to the fasta sequence
#####################################
myfasta_1 %>% 
  left_join(mydata_lifted_2, 
            by =c("start-1based"="start-1based",
                  "stop-1based"="end-1-based",
                  "chr"="chr")) %>%
  rename("V1"="Loci") -> myfasta_2

myfasta_2 %>% distinct(.keep_all = T) -> myfasta_3


#####################################
# get the giab straitrazor results
#####################################
myIlmn <- fread("Marshfield_rawILMN.txt", sep = "\t",header = F, 
                col.names = c("Instr","Loci","Length","Seq","FRD","RRD"))
myIlmn %>% separate(Loci, into = c("Loci","CE"), sep = ":") -> myIlmn_1

myPacb <- fread("Marshfield_rawPACB.txt",sep = "\t", header=F,
                col.names = c("Instr","Loci","Length","Seq","FRD","RRD"))
myPacb %>% separate(Loci, into = c("Loci","CE"), sep = ":") %>%
  separate(Loci, into = c("Loci",NA), sep = "_")-> myPacb_1

bind_rows(myIlmn_1,myPacb_1) -> mystrait


myPacb_1 %>% left_join(myfasta_3,by=c("Loci"="Loci")) -> myresults


##########################################################
## the following is to analyze how many of the seq match and how 
## many don't match to the truth data.
## Not everything is substr of vcf2seq, there are 4 cases to be consider:
## 0. There is exact match between vcf2seq and Seq --> 7
## 1. Seq is substring of vcf2seq --> 290
## 2. vcf2seq is substring of Seq --> 51
## 3. Seq is substring of rev_comp of vcf2seq --> 192
## 4. vcf2seq is substring of rev_comp of Seq --> 33
## And there were 5 that didn't match anything or reutrned na in vcf2seq
##########################################################

##########################################################
## 1. Seq is substring of vcf2seq
##########################################################
myresults %>% mutate(weeksWorth=ifelse(Seq==vcf2seq,"Woohoo","Sad"))->myresults_1

myresults %>% mutate(SeqISsubvcf2seq = map2_lgl(vcf2seq, Seq, str_detect)) -> foo

foo %>% group_by(Loci) %>% 
  arrange(desc(FRD), .by_group = TRUE) -> foo1

#count how many are true
foo1 %>% ungroup() %>% filter(SeqISsubvcf2seq=="TRUE") %>% count() #--this # is total loci including hetero
foo1 %>% ungroup() %>% filter(SeqISsubvcf2seq=="TRUE") %>% distinct(Loci) %>% count() #--this # is uniq loci

# choose only top3 that should includes homozyg and heterozyg
foo1 %>% top_n(3, FRD) -> foo2

# look at the ones that didn't match anything or produced false for all loci within groupby loci
foo2 %>% ungroup() %>% 
  group_by(Loci) %>% filter(!(any(SeqISsubvcf2seq=="TRUE"))) -> foo3

# the following is to find what is missing in vcf2seq
foo3 %>% ungroup() %>% distinct(Loci) -> Loci_fasle # 283
foo %>% distinct(Loci) -> Loci_allpacb # 578
foo1 %>% ungroup() %>% filter(SeqISsubvcf2seq=="TRUE") %>% distinct(Loci) -> Loci_true # 290


foo2 %>% filter(is.na(vcf2seq)) %>% count() # 5 -- 

##########################################################
## 2. vcf2seq is substring of Seq
##########################################################

myresults %>% mutate(vcf2seqISsubSeq = map2_lgl(Seq, vcf2seq, str_detect)) -> foovcf2seqISSeq

foovcf2seqISSeq %>% group_by(Loci) %>% 
  arrange(desc(FRD), .by_group = TRUE) %>%  
  filter(vcf2seqISsubSeq=="TRUE") %>% top_n(4, FRD) -> bar  

bar %>% ungroup() %>% distinct(Loci) %>% count() # the count is 51


#########################################################
## 3. Seq is substring of rev_comp of vcf2seq
##########################################################

myresults %>% 
  mutate(vcf2seqRevComp = map_chr(vcf2seq, reverse_complement)) -> foobar

foobar %>% mutate(SeqISsubvcf2seq_revComp = map2_lgl(vcf2seqRevComp, Seq, str_detect)) -> foobar1

foobar1 %>% group_by(Loci) %>%
  arrange(desc(FRD), .by_group = TRUE) %>%
  filter(SeqISsubvcf2seq_revComp == "TRUE") %>% top_n(4,FRD) -> foobar2

foobar2 %>% ungroup() %>% distinct(Loci) %>% count() # the count is 192

#########################################################
## 4. vcf2seq is substring of rev_comp of Seq
##########################################################

foobar %>% mutate(vcf2seqISsubSeq_revComp = map2_lgl(Seq, vcf2seqRevComp, str_detect)) -> foobarbar1

foobarbar1 %>% group_by(Loci) %>%
  arrange(desc(FRD), .by_group = TRUE) %>%
  filter(vcf2seqISsubSeq_revComp == "TRUE") %>% top_n(4,FRD) -> foobarbar2

foobarbar2 %>% ungroup() %>% distinct(Loci) %>% count() # count is 33

#########################################################
## combine all the four above
## exact match will be duplicate of substring match in one or the other
#########################################################

myresults %>% mutate(vcf2seqRevComp = map_chr(vcf2seq, reverse_complement)) %>%
  mutate(SeqISsubvcf2seq = map2_lgl(vcf2seq, Seq, str_detect),
         vcf2seqISsubSeq = map2_lgl(Seq, vcf2seq, str_detect),
         SeqISsubvcf2seq_revComp = map2_lgl(vcf2seqRevComp, Seq, str_detect),
         vcf2seqISsubSeq_revComp = map2_lgl(Seq, vcf2seqRevComp, str_detect)) -> myresults_1
         #exactMatch = ifelse(Seq==vcf2seq,"TRUE","FALSE")

myresults_1 %>% filter(SeqISsubvcf2seq=="FALSE" & vcf2seqISsubSeq=="FALSE"& 
          SeqISsubvcf2seq_revComp=="FALSE" & vcf2seqISsubSeq_revComp=="FALSE")  -> myresults_F

myresults_1 %>% filter(SeqISsubvcf2seq=="TRUE" | vcf2seqISsubSeq=="TRUE" |
              SeqISsubvcf2seq_revComp=="TRUE" | vcf2seqISsubSeq_revComp=="TRUE")  -> myresults_T

myresults_T %>% filter(SeqISsubvcf2seq=="TRUE") %>% distinct(Loci) %>% count()
myresults_T %>% filter(vcf2seqISsubSeq=="TRUE") %>% distinct(Loci) %>% count()
myresults_T %>% filter(SeqISsubvcf2seq_revComp=="TRUE") %>% distinct(Loci) %>% count()
myresults_T %>% filter(vcf2seqISsubSeq_revComp=="TRUE") %>% distinct(Loci) %>% count()
#foo1 %>% filter(exactMatch=="TRUE") %>% distinct(Loci) %>% count()

myresults_1 %>% filter(SeqISsubvcf2seq=="TRUE" & vcf2seqISsubSeq=="TRUE" & exactMatch=="TRUE")
myresults_T %>% distinct(Loci) %>% select(Loci) -> myresults_T_lociList
myresults_F %>% distinct(Loci) %>% select(Loci) -> myresults_F_lociList
