#!/bin/bash

Rscript ../ManhattanGG.R --assoc example.glm.linear --paper-material --y-lim 10 \
			--chr "#CHROM" --pos POS --snp ID --pval P \
			--deco example.deco --deco-color color  --deco-label label --deco-label-color labelcolor \
			--text-angle 90 \
			--suggestiveline 1e-05 --point-size 1.2 \
			--out example

