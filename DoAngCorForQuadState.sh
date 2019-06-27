#! /bin/bash

#This code runs the AngCor calculations from start to finish for the quadrupole state

rm /home/padsley/codes/AngCor/AngCorNew/AngCorNew/angcor/input/fort.2

#CHUCK3 calculations
cd /home/padsley/codes/AngCor/AngCorNew/AngCorNew/chuck3
./chuck < input/24Mg_alphaInelastic_2+.com > output/24Mg_alphaInelastic_2+.out
cp -f fort.2 /home/padsley/codes/AngCor/AngCorNew/AngCorNew/angcor/input/fort.2

#AngCor inputs generation
cd /home/padsley/codes/AngCor/AngCorNew/AngCorNew/angcor/input
pwd
g++ make_input_PR244_J_2.c -o make_input_PR244_J_2
./make_input_PR244_J_2 .

#AngCor calculations
./DoAngCorPR244J_2.sh

#Creation of the combined AngCor file
cd ../output
./make_final.sh
mv final.dat finalQuadrupole.dat

#Monte Carlo for the quadrupole state
cd /home/padsley/codes/AngCor/AngCorNew/AngCorNew/Averaging
g++ AverageAngCorResults.cpp -o AverageAngCorResults `root-config --cflags --libs` -O3
./AverageAngCorResults ../chuck3/output/24Mg_alphaInelastic_2+.out ../angcor/output/finalQuadrupole.dat QuadrupoleStateAngcorOutput.root