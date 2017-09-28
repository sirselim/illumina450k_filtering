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

### remember to include any probes which fail detection

```R
# process failed probes
detP <- detectionP(RGset)
failed <- detP > 0.01
colMeans(failed) # Fraction of failed positions per sample
sum(rowMeans(failed)>0.5) # How many positions failed in >50% of samples?
failed.probes <- rownames(detP[rowMeans(failed)>0.5,])
```

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

## _Update (170928)_ - addition of probes for EPIC/850k processing

Supplementary data from Pidsley *et al*., (2016), suggests cross-reactive and variant containing probes to filter at QC.

>Pidsley, R., Zotenko, E., Peters, T. J., Lawrence, M. G., Risbridger, G. P., Molloy, P., … Clark, S. J. (2016). *Critical evaluation of the Illumina MethylationEPIC BeadChip microarray for whole-genome DNA methylation profiling*. **Genome Biology**, 17(1), 208. https://doi.org/10.1186/s13059-016-1066-1

  - there is overlap between 450k and 850k lists, however this will not cause any issues.

### Extension to the above to filter EPIC data (can apply 450k list as well)

Combine the below with the above 450k process to flter EPIC arrays at QC stage:

```R
# probes from Pidsley 2016 (EPIC)
epic.cross1 <- read.csv('EPIC/13059_2016_1066_MOESM1_ESM.csv', head = T)
# epic.cross2 <- read.csv('EPIC/13059_2016_1066_MOESM2_ESM.csv', head = T)
# epic.cross3 <- read.csv('EPIC/13059_2016_1066_MOESM3_ESM.csv', head = T)
epic.variants1 <- read.csv('EPIC/13059_2016_1066_MOESM4_ESM.csv', head = T)
epic.variants2 <- read.csv('EPIC/13059_2016_1066_MOESM5_ESM.csv', head = T)
epic.variants3 <- read.csv('EPIC/13059_2016_1066_MOESM6_ESM.csv', head = T)
# additional filter probes
epic.add.probes <- c(as.character(epic.cross1$X), as.character(epic.variants1$PROBE), as.character(epic.variants2$PROBE), 
                     as.character(epic.variants3$PROBE))
# final list of unique probes
epic.add.probes <- unique(epic.add.probes)
```

Filtering process follows the same as above, example:

```R
# failed probes (those that fail detection)
beta_norm <- beta_norm[!(rownames(beta_norm) %in% failed.probes),]
colnames(beta_norm) <- pd$ID
# additional epic probes
beta_norm <- beta_norm[!(rownames(beta_norm) %in% epic.add.probes),]
```