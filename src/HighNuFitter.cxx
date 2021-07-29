#include <unistd.h>
#include "HighNuFCN.hxx"
#include "TLegend.h"

int nBins;
int binStep;
bool ParseArgs(int argc, char* argv[]);
void PrintSyntax();
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    if (!ParseArgs(argc, argv)) return 0;

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
//------------------------------------------------------------------------------
bool ParseArgs(int argc, char* argv[]) {
    bool status = false;
    const char* optstring = "b:s:t";
    char option;

    optind = 1;
    while (-1 != (option = getopt(argc, argv, optstring))) {
        switch (option) {
            case 'b' : 
                {
                    nBins = std::stoi(optarg);
                    break;
                }
            case 's' :
                {
                    binStep = std::stoi(optarg);
                    break;
                }
            case 't' :
                {
                    TOY = true;
                    std::cout << "toy == true" << std::endl;
                    std::cout << "corr?";
                    std::cin >> TOYCORR;
                    break;
                }
            case 'h' :
                {
                    PrintSyntax();
                    break;
                }
        }
    }

    if (nBins != 0 && binStep != 0) status = true;
    if (!status) PrintSyntax();
    return status; 
}
//------------------------------------------------------------------------------
void PrintSyntax() {
    std::cout << "./HighNuFitter\n";
    std::cout << "  -b ${number of bin} (REQUIRED)\n";
    std::cout << "  -s ${bin step size} (REQUIRED)\n";
    std::cout << "  -t                  (OPTIONAL)\n";
    std::cout << "    : use toy model\n";
    std::cout << "  -h\n";
    std::cout << "    : show this message\n";
    std::cout << std::endl;
}
