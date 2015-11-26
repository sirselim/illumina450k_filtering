# illumina450k_filtering
A collection of resources to filter 'bad' probes from the Illumina 450k methylation array (http://www.illumina.com/products/methylation_450_beadchip_kits.html).


## BOWTIE2 mapping of 450k probes
All probe sequences were mapped to the human genome (hg19) using BOWTIE2 to identify potential hybridisation issues. 

  - 33,457 probes were identified as aligning greater than once 
  - these are made available in `HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt`

## Additional non-specific probes
Chen *et al.,* identified a series of non-specific probes across the 450k design.

>Chen Y, Lemire M, Choufani S, Butcher DT, Grafodatskaya D, Zanke BW, Gallinger S, Hudson TJ, Weksberg R: *Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray.* **Epigenetics** 2013, 8:203–9.

  - there are a total of 29,233 probes
  - these are available in `48639-non-specific-probes-Illumina450k.csv`

***Note:*** there is overlap between the two probe sets.

## Example filtering strategy (in R)

```R
## generate 'bad' probes filter
# cross-reactive/non-specific
cross.react <- read.csv('48639-non-specific-probes-Illumina450k.csv', head = T, as.is = T)
cross.react.probes <- as.character(cross.react$TargetID)
# BOWTIE2 multi-mapped
multi.map <- read.csv('HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt', head = F, as.is = T)
multi.map.probes <- as.character(multi.map$V1)
# determine unique probes
filter.probes <- unique(c(cross.react.probes, multi.map.probes))
## filter the matrix of beta values (beta_norm)
## CpGs probes (IlmnID) should be rownames
# fitler out 'bad' probes
table(rownames(beta_norm) %in% filter.probes)
filter.bad <- rownames(beta_norm) %in% filter.probes
beta_norm <- beta_norm[!filter.bad,]
```

For a real-world example filtering strategy interested parties can refer to the methods section of our publication: (http://www.genomebiology.com/2015/16/1/8)