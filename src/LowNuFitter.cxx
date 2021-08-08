#include <getopt.h>
#include "TCanvas.h"

#include "LowNuFCN.hxx"

std::string inputFluxSystematic;
std::string inputData;
int numPars = 0;
double inputError = 0.1;
int printLevel = 0;
std::string fitOpt = "s";
bool ParseArgs(int argc, char* argv[]);
void PrintSyntax();
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    if (!ParseArgs(argc, argv)) return 0;

    LowNuFCN fcn(numPars, inputError, inputFluxSystematic, inputData);

    TH1D beforeFit = fcn.GetPrediction();
    RooMinuit m(fcn);
    m.setPrintLevel(printLevel);
    m.setStrategy(2);
    RooFitResult* result = m.fit(fitOpt.c_str());

    TH1D afterFit = fcn.GetPrediction();
    TH1D fitResult = fcn.GetFittingResult();
    fitResult.SetMinimum(0);
    fitResult.SetMaximum(1);
    TH1D data = fcn.GetData();

    TCanvas c;
    c.Divide(2,2);
    c.cd(1);
    beforeFit.Draw();
    c.cd(2);
    afterFit.Draw();
    c.cd(3);
    data.Draw();
    c.cd(4);
    fitResult.Draw();
    c.SaveAs("asdf.pdf");

    return 0;
}
//------------------------------------------------------------------------------
bool ParseArgs(int argc, char* argv[]) {
    bool status = false;

    int index;
    int iarg = 0;

    const struct option longopts[] =
    {
        {"flux-shift", required_argument, 0, 's'},
        {"data-file",  required_argument, 0, 'd'},
        {"fit-opt",    required_argument, 0, 'f'},
        {"print-level", required_argument, 0, 'p'},
        {"help",       no_argument,       0, 'h'},
        {0,0,0,0},
    };

    while (iarg != -1) {
        iarg = getopt_long(argc, argv, "s:d:n:e:f:p:h", longopts, &index);
        switch (iarg) {
            case 's' : 
                {
                    inputFluxSystematic = optarg;
                    std::cout << "inputFluxSystematic: " << inputFluxSystematic << std::endl;
                    break;
                }
            case 'd' :
                {
                    inputData = optarg;
                    std::cout << "inputData: " << inputData << std::endl;
                    break;
                }
            case 'n' :
                {
                    numPars = std::stoi(optarg);
                    std::cout << "number of flux systematic: " << numPars << std::endl;
                    break;
                }
            case 'e' :
                {
                    inputError = std::stod(optarg);
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

    if (numPars != 0 && !inputFluxSystematic.empty() && !inputData.empty()) 
        status = true;

    if (!status) PrintSyntax();
    return status; 
}
//------------------------------------------------------------------------------
void PrintSyntax() {
    std::cout << "./LowNuFitter\n";
    std::cout << "  -s, --flux-shift ${flux systematic file} (REQUIRED)\n";
    std::cout << "  -d, --data-file  ${data sample file}     (REQUIRED)\n";
    std::cout << "  -n ${number of systematic}               (REQUIRED)\n";
    std::cout << "    : number of flux systematic to use\n";
    std::cout << "  -e ${error}                              (OPTIONAL)\n";
    std::cout << "    : error for lownu cross section, if not set use 0\n";
    std::cout << "  -f, --fit-opt ${option}         (OPTIONAL)\n";
    std::cout << "    : s - Run Hesse first to estimate initial step size\n";
    std::cout << "      m - Run Migrad only\n";
    std::cout << "      h - Run Hesse to estimate errors\n";
    std::cout << "      v - Verbose mode\n";
    std::cout << "      l - Log parameters after each Minuit steps to file\n";
    std::cout << "      t - Activate profile timer\n";
    std::cout << "      r - Save fit result\n";
    std::cout << "      0 - Run Migrad with strategy 0\n";
    std::cout << "      default - s\n";
    std::cout << "  -p, --print-level ${int}\n";
    std::cout << "    : None =-1 , Reduced =0 , Normal =1 , ExtraForProblem =2 , Maximum =3 \n";
    std::cout << "  -h, --help\n";
    std::cout << "    : show this message\n";
    std::cout << std::endl;
}
