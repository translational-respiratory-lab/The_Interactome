#!/bin/bash

cd Bray-Curtis/
Rscript bray-curtis.R

cd ..
cd GBLM/
Rscript GBLM.R

cd ..
cd MI/
Rscript MI.R

cd ../Pearsons
Rscript pearsons.R 

cd ../Spearman
Rscript spearman.R
cd ..

echo "All done"
