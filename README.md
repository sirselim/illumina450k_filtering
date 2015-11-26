# illumina450k_filtering
A collection of resources to filter 'bad' probes from the Illumina 450k methylation array (http://www.illumina.com/products/methylation_450_beadchip_kits.html)


## BOWTIE2 mapping of 450k probes
All probe sequences were mapped to the human genome (hg19) using BOWTIE2 to identify potential hybridisation issues. 

  - 33,457 probes were identified as aligning greater than once 
  - these are made available in `HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt`

## Additional non-specific probes
Chen *et al.,* identified a series of non-specific probes across the 450k design.

Reference: Chen Y, Lemire M, Choufani S, Butcher DT, Grafodatskaya D, Zanke BW, Gallinger S, Hudson TJ, Weksberg R: *Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray.* **Epigenetics** 2013, 8:203–9.

  - there are a total of 29,233 probes
  - these are available in `48639-non-specific-probes-Illumina450k.csv`

***Note:*** there is overlap between the two probe sets.

## Example filtering strategy (in R)

```R
head(data)

```