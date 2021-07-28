#ifndef __LOWNUFCN__
#define __LOWNUFCN__

#include <TMatrixT.h>
#include <TVectorT.h>
#include <TH1.h>
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>

#include "RooAbsReal.h"
#include "RooListProxy.h"
#include "RooRealVar.h"
#include "RooMinuit.h"

struct Event {
    float trueNeutrinoE = -1;
    float recoNeutrinoE = -1;
    bool isSignal = false;
};

class LowNuFCN : public RooAbsReal {
    public:
        LowNuFCN() {};
        LowNuFCN(int numPars, double inError, std::string inputFluxSystematic,
                std::string inputDataFile);
        ~LowNuFCN();
        virtual TObject* clone(const char* newname=0) const {return new LowNuFCN;};
        virtual Double_t evaluate() const;

        const int GetNumberOfParameters() const {return this->mNumberOfParameters;};

    private:
        RooArgList _parlist;
        RooListProxy* _pulls;
        std::string mFluxSystematicFileName;
        std::string mDataFileName;
        TH1D* mData;
        const int mBins = 16;
        int mNumberOfParameters;
        double mError;
        const float mNuCut = 200;
        std::vector<TH1D> mFluxSyst;
        std::vector<Event> mEvents;

        TFile* dataFile = nullptr;
        TTree* inputTree = nullptr;
        float recoNeutrinoE;
        float trueNeutrinoE;
        bool isSignal;

        TVectorD* pullCV; 
        TVectorD* pullUnc;

        void LoadFluxSystematics(std::string fluxSystematicFile);
        void LoadData(std::string dataFile);
        void SetInputTree();
        void FillData();
        void SetPullCV();
        void SetPullUnc();
        TMatrixD* PrepareCovMatrix(Int_t nBins, TH1D pred) const;
        TMatrixD* PrepareCovMatrix2(Int_t nBins) const;
        TH1D PreparePrediction(RooListProxy* _pulls) const; 
        Double_t FillEv(RooListProxy* _pulls) const;
        Double_t FillEv2(RooListProxy* _pulls) const;
        Double_t ExtraPull(RooListProxy* _pulls) const;
};

#endif
