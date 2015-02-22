#!/usr/bin/bash

### run :  sh readLength.sh par1 par2

### par1 = 0 > new mf
### par1 = 1 > old mf

### par2 = 0 > simple
### par2 = 1 > moderate
### par2 = 2 > Hard


sh readLength.sh 0
sh readLength.sh 1
sh readLength.sh 2



./readLengthPlot.r 0
./readLengthPlot.r 1
./readLengthPlot.r 2
