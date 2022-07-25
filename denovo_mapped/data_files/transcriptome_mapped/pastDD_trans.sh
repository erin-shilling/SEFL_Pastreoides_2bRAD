#!/bin/bash
angsd -b bamsClones -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 960 -minInd 47 -doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2 -P 1 -out dd
