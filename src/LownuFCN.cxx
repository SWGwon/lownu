#include "LowNuFCN.hxx"

LownuFCN::LownuFCN(int numPars, double inError)
    : mNumberOfParameters(numPars), mError(inError) {
    _pulls = new RooListProxy("_pulls", "_pulls", this);

}
//-----------------------------------------------------------------------------
Double_t LownuFCN::evaluate() const {
    double chi2;
    return chi2;
}
//-----------------------------------------------------------------------------
