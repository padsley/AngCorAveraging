#! /bin/bash

#This code runs the AngCor calculations from start to finish for the dipole state

#CHUCK3 calculations


#AngCor inputs generation
cd /home/padsley/codes/AngCor/AngCorNew/AngCorNew/angcor/input
pwd
g++ make_input_PR244_J_1.c -o make_input_PR244_J_1
./make_input_PR244_J_1 .

#AngCor calculations
./DoAngCorPR244J_1.sh

#Creation of the combined AngCor file
cd ../output
./make_final.sh

#Monte Carlo for the dipole state
cd /home/padsley/codes/AngCor/AngCorNew/AngCorNew/Averaging
g++ AverageAngCorResults.cpp -o AverageAngCorResults `root-config --cflags --libs` -O3
./AverageAngCorResults ../chuck3/output/24Mg_alphaInelastic_1-.out ../angcor/output/final.dat DipoleStateAngcorOutput.root