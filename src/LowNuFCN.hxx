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
        virtual TObject* clone(const char* newname=0) const {return new LowNuFCN;};
        virtual Double_t evaluate() const;

        const int GetNumberOfParameters() const {return this->mNumberOfParameters;};
        TH1D GetPrediction() const; 
        TH1D GetData() {return this->mData;};
        TMatrixD* GetCovMatrix() const {return this->mCovMat.get();};
        TH1D GetFittingResult();

    private:
        const int mBins = 16;
        int mNumberOfParameters;
        double mError;
        const float mNuCut = 200;
        TH1D mData;
        std::string mFluxSystematicFileName;
        std::string mDataFileName;
        std::vector<TH1D> mFluxSyst;
        std::vector<Event> mEvents;
        std::vector<RooRealVar*> mFluxPars;
        std::vector<RooRealVar*> mParsEnergyBins;

        std::unique_ptr<RooListProxy> mPulls;
        std::unique_ptr<RooListProxy> mPullsbkg;
        std::unique_ptr<TMatrixD> mCovMat;
        std::unique_ptr<TVectorD> pullCV; 
        std::unique_ptr<TVectorD> pullUnc;
        std::unique_ptr<TFile> dataFile;

        TTree* inputTree;
        float recoNeutrinoE;
        float trueNeutrinoE;
        bool isSignal;

        void LoadFluxSystematics(std::string fluxSystematicFile);
        void LoadData(std::string dataFile);
        void SetInputTree();
        void FillData();
        void SetPullCV();
        void SetPullUnc();
        void SetCovMatrix();
        double GetWeight(int inBin) const;
        Double_t FillEv() const;
        Double_t ExtraPull() const;
};

#endif
