#include "LowNuFCN.hxx"

int main(int argc, char* argv[])
{
    std::string inputFluxSystematic = argv[1];
    std::string inputData = argv[2];

    LowNuFCN fcn(10,0, inputFluxSystematic, inputData);
    RooMinuit m(fcn);
    m.setStrategy(2);
    RooFitResult* result = m.fit("s");

    return 0;
}
