/**
adopted from the file OpenMM6.0-Source/platforms/cpu/tests/TestCpuNeighborList.cpp
compiled in there using make command

 **/


#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/ThreadPool.h"
#include "AlignedArray.h"
#include "CpuNeighborList.h"
#include "CpuPlatform.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <set>
#include <utility>
#include <vector>
#include <fcntl.h>


using namespace OpenMM;
using namespace std;

void testNeighborList(int numParticles, AlignedArray<float> positions, float cutoff) {
       
}


double sq(double x)
{ return x*x;}

float sq(float x) 
{return x*x; }

int main(int argc, char * argv[]) {

freopen(NULL, "wb", stdout);
freopen(NULL, "rb", stdin);


int N ;
float dummy;
float cutoff; 
char toprint[sizeof(int)];

cutoff = atof(argv[1]);
N = atoi(argv[2]);
//dummy=fscanf(stdin,"%f",&cutoff);
//dummy=fscanf(stdin,"%d",&N);


double cutoff2 = cutoff * cutoff; 
AlignedArray<float> positions(4*N);






for(int i=0;i<N;i++)
        {
//        fscanf(stdin,"%f %f %f",&positions[4*i],
//                &positions[4*i+1],&positions[4*i+2]);
        for (int j=0; j<3; j++){
            fread(&positions[4*i+j], 1, sizeof(float), stdin);
//            cerr << positions[4*i+j] << endl;


        }
        
        
        }

    try {
        if (!CpuPlatform::isProcessorSupported()) {
            cout << "CPU is not supported.  Exiting." << endl;
            return 1;
        }
 
     float boxSize[3];    
    const int blockSize = 8;
    boxSize[0] = 1000;
    boxSize[1] = 1000;
    boxSize[2] = 1000;     
    vector<set<int> > exclusions(N);
    ThreadPool threads;
    CpuNeighborList neighborList(blockSize);
    neighborList.computeNeighborList(N, positions, exclusions, boxSize, false, cutoff,  threads);

    for (int i = 0; i < (int) neighborList.getSortedAtoms().size(); i++) {
        int blockIndex = i/blockSize;
        int indexInBlock = i-blockIndex*blockSize;
        char mask = 1<<indexInBlock;
        for (int j = 0; j < (int) neighborList.getBlockExclusions(blockIndex).size(); j++) {
            if ((neighborList.getBlockExclusions(blockIndex)[j] & mask) == 0) {
                int atom1 = neighborList.getSortedAtoms()[i];
                int atom2 = neighborList.getBlockNeighbors(blockIndex)[j];
                pair<int, int> entry = make_pair(min(atom1, atom2), max(atom1, atom2));
                int first = entry.first;
                int second = entry.second;
                if (sq(positions[4*first] - positions[4*second]) + sq(positions[4*first+1] - positions[4*second+1]) + sq(positions[4*first+2] - positions[4*second+2]) < cutoff2){
                    
                    std::cout.write(reinterpret_cast<const char*>(&first), sizeof first);
                    std::cout.write(reinterpret_cast<const char*>(&second), sizeof second);
                }                                
            }
        }
    }

 
 
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    return 0; 
    
}
