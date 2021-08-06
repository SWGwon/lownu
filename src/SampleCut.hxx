#ifndef __SAMPLECUT__
#define __SAMPLECUT__

#include <iostream>
#include <iomanip>
#include "TFile.h"
#include "TString.h"
#include "TTree.h"

TFile* file = nullptr;
TTree* tree = nullptr;
std::string outName;

int category;
int numberOfFSNeutron;
int numberOfBranches;
float trackLength;
float recoNu;
float trueNu;
float tof;
float maxAngle;
float maxDistance;
float neighborDistance;
float eDep;

float tofCutValue = 0;
float eDepCutValueCluster = 735;
float eDepCutValueTrack = 3700;
float neighborDCutValue = 20;
float numBranchCutValue = 0;
float trackLCutValue = 20;
float lownuCutValue = 200;
float maxACutValue = 0.7536;
float maxDCutValue = 700;

bool tofCut = true;
bool eDepCut = true;
bool neighborDCut = false;
bool numBranchCut = true;
bool trackLCut = false;
bool lownuCut = false;
bool maxACut = true;
bool maxDCut = true;

double totalNumberOfEvents = 0;
double sig[10] = {};
double bkg[10] = {};

void SampleCut(std::string fileName);

int barSize = 60;

void SetTree();

void PrintProgramInfo();

void PrintSyntax();

void PrintCutInfo();

void PrintInputInfo();

void PrintProgress(int i, int total);

void PrintResult();

void PrintBar(int i);

void GetCommandLineArgs(int argc, char* argv[]);

#endif
