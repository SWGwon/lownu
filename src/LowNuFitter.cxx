#include <getopt.h>
#include "TCanvas.h"

#include "LowNuFCN.hxx"

std::string inputFluxSystematic;
std::string inputData;
int numPars = 0;
double inputError = 0;
bool ParseArgs(int argc, char* argv[]);
void PrintSyntax();
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    if (!ParseArgs(argc, argv)) return 0;

    LowNuFCN fcn(numPars, inputError, inputFluxSystematic, inputData);

    TH1D beforeFit = fcn.GetPrediction(fcn.GetPull());
    RooMinuit m(fcn);
    m.setStrategy(2);
    RooFitResult* result = m.fit("s");

    TH1D afterFit = fcn.GetPrediction(fcn.GetPull());
    TH1D data = fcn.GetData();

    TCanvas c;
    c.Divide(2,2);
    c.cd(1);
    beforeFit.Draw();
    c.cd(2);
    afterFit.Draw();
    c.cd(3);
    data.Draw();
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
        {"help",       no_argument,       0, 'h'},
        {0,0,0,0},
    };

    while (iarg != -1) {
        iarg = getopt_long(argc, argv, "s:d:n:e:h", longopts, &index);
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
    std::cout << "  -h, --help\n";
    std::cout << "    : show this message\n";
    std::cout << std::endl;
}
