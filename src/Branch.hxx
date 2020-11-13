#ifndef BRANCH_H
#define BRANCH_H

#include <memory>

#include "TChain.h"

extern int t_EvtNum;
extern double t_EvtXSec;
extern double t_EvtDXSec;
extern double t_EvtWght;
extern double t_EvtProb;
extern double t_EvtVtx[4];
extern int t_StdHepN;
extern int t_StdHepPdg[1000];
extern int t_StdHepStatus[1000];
extern int t_StdHepRescat[1000];
extern double t_StdHepX4[1000][4];
extern double t_StdHepP4[1000][4];
extern double t_StdHepPolz[1000][3];
extern int t_StdHepFd[1000];
extern int t_StdHepLd[1000];
extern int t_StdHepFm[1000];
extern int t_StdHepLm[1000];

namespace Branch
{
    void SetBranchAddress(const std::unique_ptr<TChain>&  tree);
}

#endif
