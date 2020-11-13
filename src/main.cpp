#include <iostream>
#include <memory>

#include "TMath.h"
#include "TH2.h"
#include "TCanvas.h"

#include "Functions.hxx"
#include "Variables.hxx"
#include "Branch.hxx"

int main(int argc, char * argv[])
{
    auto QE_Q2_vs_nu = std::make_unique<TH2F> ("","QE;Q^{2};#nu",100,0,5.5,100,0,5.5);
    auto RES_Q2_vs_nu = std::make_unique<TH2F> ("","RES;Q^{2};#nu",100,0,5.5,100,0,5.5);
    auto DIS_Q2_vs_nu = std::make_unique<TH2F> ("","DIS;Q^{2};#nu",100,0,5.5,100,0,5.5);

    auto QE_Q2_vs_nuE = std::make_unique<TH2F> ("","QE;Q^{2};E_{#nu}",100,0,5.5,100,0,5.5);
    auto RES_Q2_vs_nuE = std::make_unique<TH2F> ("","RES;Q^{2};E_{#nu}",100,0,5.5,100,0,5.5);
    auto DIS_Q2_vs_nuE = std::make_unique<TH2F> ("","DIS;Q^{2};E_{#nu}",100,0,5.5,100,0,5.5);

    auto QE_y_vs_nu = std::make_unique<TH2F> ("","QE;y;#nu",100,0,1,100,0,5.5);
    auto RES_y_vs_nu = std::make_unique<TH2F> ("","RES;y;#nu",100,0,1,100,0,5.5);
    auto DIS_y_vs_nu = std::make_unique<TH2F> ("","DIS;y;#nu",100,0,1,100,0,5.5);

    auto QE_y_vs_Q2 = std::make_unique<TH2F> ("","QE;y;Q^{2}",100,0,1,100,0,5.5);
    auto RES_y_vs_Q2 = std::make_unique<TH2F> ("","RES;y;Q^{2}",100,0,1,100,0,5.5);
    auto DIS_y_vs_Q2 = std::make_unique<TH2F> ("","DIS;y;Q^{2}",100,0,1,100,0,5.5);

    auto QE_p2_vs_Q2 = std::make_unique<TH2F> ("","QE;p^{2};Q^{2}",100,0,5.5,100,0,5.5);
    auto RES_p2_vs_Q2 = std::make_unique<TH2F> ("","RES;p^{2};Q^{2}",100,0,5.5,100,0,5.5);
    auto DIS_p2_vs_Q2 = std::make_unique<TH2F> ("","DIS;p^{2};Q^{2}",100,0,5.5,100,0,5.5);

    auto QE_p2_vs_nu = std::make_unique<TH2F> ("","QE;p^{2};#nu",100,0,5.5,100,0,5.5);
    auto RES_p2_vs_nu = std::make_unique<TH2F> ("","RES;p^{2};#nu",100,0,5.5,100,0,5.5);
    auto DIS_p2_vs_nu = std::make_unique<TH2F> ("","DIS;p^{2};#nu",100,0,5.5,100,0,5.5);
    if(!Function::Parsing(argc, argv))
    {
        return 0;
    }

    auto tree = std::make_unique<TChain> ("gRooTracker");

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
        float y = 0;
        float final_state_momentum[3] = {0,0,0};
        float p2 = 0;
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
            if(t_StdHepStatus[FS] == 1 && abs(t_StdHepPdg[FS]) != 13)
            {
                final_state_momentum[0] = final_state_momentum[0]+t_StdHepP4[FS][0];
                final_state_momentum[1] = final_state_momentum[1]+t_StdHepP4[FS][1];
                final_state_momentum[2] = final_state_momentum[2]+t_StdHepP4[FS][2];
            }
        }
        if(!isCC)
            continue;
        nu = t_StdHepP4[0][3] - muonE;
        Q2 = 4 * t_StdHepP4[0][3] * muonE * pow(TMath::Sin(muonAngle/2),2);
        y = nu/t_StdHepP4[0][3];
        p2 = pow(pow(final_state_momentum[0],2)+pow(final_state_momentum[1],2)+pow(final_state_momentum[2],2),0.5);
        std::cout<<p2<<std::endl;
        if(num_pi+num_pi0 == 0)
        {
            QE_Q2_vs_nu->Fill(Q2,nu);
            QE_Q2_vs_nuE->Fill(Q2,t_StdHepP4[0][3]);
            QE_y_vs_nu->Fill(y,nu);
            QE_y_vs_Q2->Fill(y,Q2);
            QE_p2_vs_Q2->Fill(p2,Q2);
            QE_p2_vs_nu->Fill(p2,nu);
        }
        if(num_pi+num_pi0 == 1)
        {
            RES_Q2_vs_nu->Fill(Q2,nu);
            RES_Q2_vs_nuE->Fill(Q2,t_StdHepP4[0][3]);
            RES_y_vs_nu->Fill(y,nu);
            RES_y_vs_Q2->Fill(y,Q2);
            RES_p2_vs_Q2->Fill(p2,Q2);
            RES_p2_vs_nu->Fill(p2,nu);
        }
        if(num_pi+num_pi0 > 1)
        {
            DIS_Q2_vs_nu->Fill(Q2,nu);
            DIS_Q2_vs_nuE->Fill(Q2,t_StdHepP4[0][3]);
            DIS_y_vs_nu->Fill(y,nu);
            DIS_y_vs_Q2->Fill(y,Q2);
            DIS_p2_vs_Q2->Fill(p2,Q2);
            DIS_p2_vs_nu->Fill(p2,nu);
        }
    }
    TCanvas can1;
    can1.Divide(2,2);
    can1.cd(1);
    QE_Q2_vs_nu->Draw("colz");
    can1.cd(2);
    RES_Q2_vs_nu->Draw("colz");
    can1.cd(3);
    DIS_Q2_vs_nu->Draw("colz");
    can1.cd(4);
    can1.SaveAs("Q2_vs_nu.pdf");
    can1.SaveAs("Q2_vs_nu.C");

    TCanvas can2;
    can2.Divide(2,2);
    can2.cd(1);
    QE_Q2_vs_nuE->Draw("colz");
    can2.cd(2);
    RES_Q2_vs_nuE->Draw("colz");
    can2.cd(3);
    DIS_Q2_vs_nuE->Draw("colz");
    can2.cd(4);
    can2.SaveAs("Q2_vs_nuE.pdf");
    can2.SaveAs("Q2_vs_nuE.C");

    TCanvas can3;
    can3.Divide(2,2);
    can3.cd(1);
    QE_y_vs_nu->Draw("colz");
    can3.cd(2);
    RES_y_vs_nu->Draw("colz");
    can3.cd(3);
    DIS_y_vs_nu->Draw("colz");
    can3.cd(4);
    can3.SaveAs("y_vs_nu.pdf");
    can3.SaveAs("y_vs_nu.C");

    TCanvas can4;
    can4.Divide(2,2);
    can4.cd(1);
    QE_p2_vs_Q2->Draw("colz");
    can4.cd(2);
    RES_p2_vs_Q2->Draw("colz");
    can4.cd(3);
    DIS_p2_vs_Q2->Draw("colz");
    can4.cd(4);
    can4.SaveAs("p2_vs_Q2.pdf");
    can4.SaveAs("p2_vs_Q2.C");

    TCanvas can5;
    can5.Divide(2,2);
    can5.cd(1);
    QE_y_vs_Q2->Draw("colz");
    can5.cd(2);
    RES_y_vs_Q2->Draw("colz");
    can5.cd(3);
    DIS_y_vs_Q2->Draw("colz");
    can5.cd(4);
    can5.SaveAs("y_vs_Q2.pdf");
    can5.SaveAs("y_vs_Q2.C");

    TCanvas can6;
    can6.Divide(2,2);
    can6.cd(1);
    QE_p2_vs_nu->Draw("colz");
    can6.cd(2);
    RES_p2_vs_nu->Draw("colz");
    can6.cd(3);
    DIS_p2_vs_nu->Draw("colz");
    can6.cd(4);
    can6.SaveAs("p2_vs_nu.pdf");
    can6.SaveAs("p2_vs_nu.C");

    return 0;
}
