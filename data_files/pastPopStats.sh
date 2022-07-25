#!/bin/bash
srun angsd -b bamsNoClones -GL 1 -uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 30 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 27 -doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doGlf 2 -P 1 -out pastPopStats
