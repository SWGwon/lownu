#ifndef __LOWNUFCN__
#define __LOWNUFCN__

#include "RooAbsReal.h"
#include "RooListProxy.h"

class LownuFCN : public RooAbsReal {
    public:
        LownuFCN() {};
        LownuFCN(int numPars, double inError);
        virtual TObject* clone(const char* newname=0) const {return new LownuFCN;};
        virtual Double_t evaluate() const;

    private:
        RooListProxy* _pulls;
        int mNumberOfParameters;
        double mError;
};

#endif
