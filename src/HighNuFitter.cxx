#include <getopt.h>
#include "HighNuFCN.hxx"
#include "TLegend.h"
#include "TGraph.h"

int nBins;
int binStep;
int printLevel = 0;
bool ParseArgs(int argc, char* argv[]);
void PrintSyntax();
std::string fitOpt = "";
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    if (!ParseArgs(argc, argv)) return 0;

    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("2.3f");

    HighNuFCN fcn(nBins, binStep);
    RooMinuit m(fcn);
    m.setPrintLevel(printLevel);
    m.setStrategy(2);
    RooFitResult* result = m.fit(fitOpt.c_str());

    TMatrixD* toyCorr = fcn.GetToyCorrelationMatrix();
    TMatrixD* toyCov = fcn.GetToyCovarianceMatrix();

    TH1D* nominal = fcn.GetHistGenieNominal();
    nominal->Sumw2(false);
    TH1D* shift = fcn.GetHistGenieShift();
    shift->Sumw2(false);
    TH1D* sample = fcn.GetHistSampleResult();
    //TMatrixD corr(*fcn.GetToyCorrelationMatrix());
    TMatrixD corr(*fcn.GetCorrelationMatrix());
    TMatrixD cov(*fcn.GetCovarianceMatrix());

    TCanvas can;
    can.Divide(2,2);
    can.cd(1);
    nominal->SetLineColor(2);
    nominal->Draw();
    shift->Draw("same");
    can.cd(1)->BuildLegend(0.7,0.7, 0.9,0.9);
    can.cd(2);
    sample->Draw();
    can.cd(3);
    corr.Draw("text colz");
    can.cd(4);
    TMatrixD invertCov(corr.Invert());
    cov.Invert().Draw("text colz");
    can.SaveAs("asd.pdf");

    TCanvas c2;
    fcn.GetHistGenieShift()->Draw();
    fcn.GetHistGenieNominal()->Sumw2(false);
    //fcn.GetHistGenieNominal()->SetLineColor(2);
    fcn.GetHistGenieNominal()->Draw("same");
    TLegend l(0.7,0.7, 0.9,0.9);
    l.AddEntry(fcn.GetHistGenieNominal(), "nominal");
    l.AddEntry(fcn.GetHistGenieShift(), "shift");
    l.Draw();
    c2.SaveAs("slide.pdf");

    const int size = fcn.mParV.size();
    double x[size];
    double y[size];
    for (int i = 0; i < size; ++i) {
        x[i] = fcn.mParV.at(i);
        y[i] = fcn.mChi2.at(i);
    }
    TGraph gr(size, x, y);
    TCanvas c3;
    gr.Draw("A*");
    c3.SaveAs("chi2_dist.pdf");


    //fcn.SaveHist("modelHist");
}
//------------------------------------------------------------------------------
bool ParseArgs(int argc, char* argv[]) {
    bool status = false;

    int index;
    int iarg = 0;
    const struct option longopts[] =
    {
        {"num-bin", required_argument, 0, 'b'},
        {"bin-step",  required_argument, 0, 's'},
        {"toy",       no_argument,       0, 't'},
        {"fit-opt",   required_argument, 0, 'f'},
        {"print-level", required_argument, 0, 'p'},
        {"help",       no_argument,       0, 'h'},
        {0,0,0,0},
    };

    while (iarg != -1) {
        iarg = getopt_long(argc, argv, "b:s:tf:p:h", longopts, &index);
        switch (iarg) {
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
            case 'f' :
                {
                    fitOpt = optarg;
                    break;
                }
            case 'p' :
                {
                    printLevel = std::stoi(optarg);
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
    std::cout << "  -b, --num-bin ${number of bin}  (REQUIRED)\n";
    std::cout << "  -s, --bin-step ${bin step size} (REQUIRED)\n";
    std::cout << "  -t, --toy                       (OPTIONAL)\n";
    std::cout << "    : use toy model\n";
    std::cout << "  -f, --fit-opt ${option}         (OPTIONAL)\n";
    std::cout << "    : s - Run Hesse first to estimate initial step size\n";
    std::cout << "      m - Run Migrad only\n";
    std::cout << "      h - Run Hesse to estimate errors\n";
    std::cout << "      v - Verbose mode\n";
    std::cout << "      l - Log parameters after each Minuit steps to file\n";
    std::cout << "      t - Activate profile timer\n";
    std::cout << "      r - Save fit result\n";
    std::cout << "      0 - Run Migrad with strategy 0\n";
    std::cout << "      (default) : migrad() + minos()\n";
    std::cout << "  -p, --print-level ${int}\n";
    std::cout << "    : None =-1 , Reduced =0 , Normal =1 , ExtraForProblem =2 , Maximum =3 \n";
    std::cout << "  -h, --help\n";
    std::cout << "    : show this message\n";
    std::cout << std::endl;
}
