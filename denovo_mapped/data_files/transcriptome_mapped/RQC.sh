#!/bin/bash
Rscript /mnt/beegfs/home/reckert2017/bin/plotQC.R prefix=dd
gzip -9 dd.counts
