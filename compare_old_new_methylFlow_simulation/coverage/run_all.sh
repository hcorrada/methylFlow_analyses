#!/usr/bin/bash

### run :  sh coverage.sh par1 par2

### par1 = 0 > new mf
### par1 = 1 > old mf

### par2 = 0 > simple
### par2 = 1 > moderate
### par2 = 2 > Hard


sh coverage.sh 0
./coveragePlot.r 0

sh coverage.sh 1
./coveragePlot.r 1

sh coverage.sh 2
./coveragePlot.r 2



