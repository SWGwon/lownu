#include "TCanvas.h"

#include "LowNuFCN.hxx"

int main(int argc, char* argv[])
{
    std::string inputFluxSystematic = argv[1];
    std::string inputData = argv[2];

    LowNuFCN fcn(5,0, inputFluxSystematic, inputData);

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
