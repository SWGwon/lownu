#include <iostream>

#include "TMath.h"
#include "TH2.h"

#include "Functions.hxx"
#include "Variables.hxx"
#include "Branch.hxx"

int main(int argc, char * argv[])
{
    TH2F * QE_Q2_vs_nu = new TH2F("",";Q^2;nu",100,0,5.5,100,0,5.5);
    TH2F * RES_Q2_vs_nu = new TH2F("",";Q^2;nu",100,0,5.5,100,0,5.5);
    TH2F * DIS_Q2_vs_nu = new TH2F("",";Q^2;nu",100,0,5.5,100,0,5.5);
    TH2F * QE_Q2_vs_nuE = new TH2F("",";Q^2;nuE",100,0,5.5,100,0,5.5);
    TH2F * RES_Q2_vs_nuE = new TH2F("",";Q^2;nuE",100,0,5.5,100,0,5.5);
    TH2F * DIS_Q2_vs_nuE = new TH2F("",";Q^2;nuE",100,0,5.5,100,0,5.5);
    if(!Function::Parsing(argc, argv))
    {
        return 0;
    }

    TChain * tree = new TChain("gRooTracker");

    if(!Function::InputFile(tree))
    {
        return 0;
    }

    Branch::SetBranchAddress(tree);

    for(int i = 0; i < tree->GetEntries(); i++)
    {
        tree->GetEntry(i);
        float Q2 = 0;
        float nu = 0;
        int num_pi = 0;
        int num_pi0 = 0;
        float muonE = 0;
        float norm_muonMomentum[3];
        float muonAngle = 0;
        bool isCC = false;
        if(t_StdHepP4[0][3] == 0)
            continue;
        for(int FS = 0; FS < t_StdHepN; FS++)
        {
            if(t_StdHepPdg[FS] == 13 && t_StdHepStatus[FS] == 1 && isCC == false)
            {
                isCC = true;
                muonE = t_StdHepP4[FS][3];
                norm_muonMomentum[2] = t_StdHepP4[FS][2]/pow(t_StdHepP4[FS][0]*t_StdHepP4[FS][0]+t_StdHepP4[FS][1]*t_StdHepP4[FS][1]+t_StdHepP4[FS][2]*t_StdHepP4[FS][2],0.5);
                muonAngle = TMath::ACos(norm_muonMomentum[2]);
            }
            if(abs(t_StdHepPdg[FS]) == 211 && t_StdHepStatus[FS] == 1)
            {
                num_pi++;
            }
            if(t_StdHepPdg[FS] == 111 && t_StdHepStatus[FS] == 1)
            {
                num_pi0++;
            }
        }
        if(!isCC)
            continue;
        nu = t_StdHepP4[0][3] - muonE;
        Q2 = 4 * t_StdHepP4[0][3] * muonE * pow(TMath::Sin(muonAngle/2),2);
        if(num_pi+num_pi0 == 0)
        {
            QE_Q2_vs_nu->Fill(Q2,nu);
            QE_Q2_vs_nuE->Fill(Q2,t_StdHepP4[0][3]);
        }
        if(num_pi+num_pi0 == 1)
        {
            RES_Q2_vs_nu->Fill(Q2,nu);
            RES_Q2_vs_nuE->Fill(Q2,t_StdHepP4[0][3]);
        }
        if(num_pi+num_pi0 > 1)
        {
            DIS_Q2_vs_nu->Fill(Q2,nu);
            DIS_Q2_vs_nuE->Fill(Q2,t_StdHepP4[0][3]);
        }
    }
    QE_Q2_vs_nu->SaveAs("QE_Q2_vs_nu.C");
    RES_Q2_vs_nu->SaveAs("RES_Q2_vs_nu.C");
    DIS_Q2_vs_nu->SaveAs("DIS_Q2_vs_nu.C");
    QE_Q2_vs_nuE->SaveAs("QE_Q2_vs_nuE.C");
    RES_Q2_vs_nuE->SaveAs("RES_Q2_vs_nuE.C");
    DIS_Q2_vs_nuE->SaveAs("DIS_Q2_vs_nuE.C");

    return 0;
}
