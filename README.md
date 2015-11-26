# illumina450k_filtering
A collection of resources to filter 'bad' probes from the Illumina 450k methylation array (http://www.illumina.com/products/methylation_450_beadchip_kits.html)


## BOWTIE2 mapping of 450k probes
All probe sequences were mapped to the human genome (hg19) using BOWTIE2 to identify potential hybridisation issues. 

  - 33,457 probes were identified as aligning greater than once 
  - these are made available in `HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt`