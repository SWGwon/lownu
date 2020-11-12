#include "Branch.hxx"

int t_EvtNum;
double t_EvtXSec;
double t_EvtDXSec;
double t_EvtWght;
double t_EvtProb;
double t_EvtVtx[4];
int t_StdHepN;
int t_StdHepPdg[1000];
int t_StdHepStatus[1000];
int t_StdHepRescat[1000];
double t_StdHepX4[1000][4];
double t_StdHepP4[1000][4];
double t_StdHepPolz[1000][3];
int t_StdHepFd[1000];
int t_StdHepLd[1000];
int t_StdHepFm[1000];
int t_StdHepLm[1000];

void Branch::SetBranchAddress(TChain * tree)
{
    tree->SetBranchAddress("EvtNum",&t_EvtNum);
    tree->SetBranchAddress("EvtXSec",&t_EvtXSec);
    tree->SetBranchAddress("EvtDXSec",&t_EvtDXSec);
    tree->SetBranchAddress("EvtWght",&t_EvtWght);
    tree->SetBranchAddress("EvtProb",&t_EvtProb);
    tree->SetBranchAddress("EvtVtx",&t_EvtVtx);
    tree->SetBranchAddress("StdHepN",&t_StdHepN);
    tree->SetBranchAddress("StdHepPdg",&t_StdHepPdg);
    tree->SetBranchAddress("StdHepStatus",&t_StdHepStatus);
    tree->SetBranchAddress("StdHepRescat",&t_StdHepRescat);
    tree->SetBranchAddress("StdHepX4",&t_StdHepX4);
    tree->SetBranchAddress("StdHepP4",&t_StdHepP4);
    tree->SetBranchAddress("StdHepPolz",&t_StdHepPolz);
    tree->SetBranchAddress("StdHepFd",&t_StdHepFd);
    tree->SetBranchAddress("StdHepLd",&t_StdHepLd);
    tree->SetBranchAddress("StdHepFm",&t_StdHepFm);
    tree->SetBranchAddress("StdHepLm",&t_StdHepLm);
}
