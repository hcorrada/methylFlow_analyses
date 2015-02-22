#!/usr/bin/bash

### run :  sh cpg.sh par1 par2

### par1 = 0 > new mf
### par1 = 1 > old mf

### par2 = 0 > simple
### par2 = 1 > moderate
### par2 = 2 > Hard


sh cpg.sh 0
./cpgPlot.r 0

sh cpg.sh 1
./cpgPlot.r 1


sh cpg.sh 2
./cpgPlot.r 2




