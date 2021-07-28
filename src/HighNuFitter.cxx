#include "HighNuFCN.hxx"
#include "TLegend.h"

int main(int argc, char** argv)
{
    int nBins;
    int binStep;
    if (argv[1] && argv[2]) {
        nBins = std::stoi(argv[1]);
        binStep = std::stoi(argv[2]);
        std::cout << "bin number: " << nBins <<  std::endl;
        std::cout << "bin step: " << binStep <<  std::endl;
    } else {
        std::cout << "invalide args" << std::endl;
        std::cout << "./HighnuFit binNumber binStep" << std::endl;
        return 0;
    }
    if (argv[3]) {
        TOY = argv[3];
        std::cout << "toy == true" << std::endl;
        std::cout << "corr?";
        std::cin >> TOYCORR;
    }

    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("2.3f");

    FCN fcn(nBins, binStep);
    RooMinuit m(fcn);
    m.setStrategy(2);
    RooFitResult* result = m.fit("s");

    TMatrixD* toyCorr = fcn.GetToyCorrelationMatrix();
    TMatrixD* toyCov = fcn.GetToyCovarianceMatrix();
    TMatrixD* corr = fcn.GetCorrelationMatrix();
    TMatrixD* cov = fcn.GetCovarianceMatrix();

    TCanvas can;
    can.Divide(2,2);
    can.cd(1);
    toyCorr->Draw("text colz");
    can.cd(2);
    toyCov->Draw("text colz");
    can.cd(3);
    corr->Draw("text colz");
    can.cd(4);
    cov->Draw("text colz");
    can.SaveAs("asd.pdf");

    TCanvas c2;
    fcn.GetHistGenieShift()->DrawNormalized();
    fcn.GetHistGenieNominal()->Sumw2(false);
    //fcn.GetHistGenieNominal()->SetLineColor(2);
    fcn.GetHistGenieNominal()->DrawNormalized("same");
    TLegend l(0.7,0.7, 0.9,0.9);
    l.AddEntry(fcn.GetHistGenieNominal(), "nominal");
    l.AddEntry(fcn.GetHistGenieShift(), "shift");
    l.Draw();
    c2.SaveAs("slide.pdf");

    fcn.SaveHist("modelHist");
}
