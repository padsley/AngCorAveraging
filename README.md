# AngCorAveraging
My code to do averaging over the K600's solid angle for the AngCor results

The code is compiled using something approximately like:

g++ AverageAngCorResults.cpp -o AverageAngCorResults `root-config --cflags --libs` -O3

It's run using:

./AverageAngCorResults <CHUCK3CrossSectionFile> <AngCorCombinedResults> <NameOfOutputROOTfile.root>
