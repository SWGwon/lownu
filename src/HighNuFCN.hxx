#ifndef __HIGHNUFCN__
#define __HIGHNUFCN__

#include <memory>

#include "TCanvas.h"
#include "TStyle.h"
#include "TVectorD.h"
#include "TRandom.h"
#include "TMatrix.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "RooAbsReal.h"
#include "RooListProxy.h"
#include "RooRealVar.h"
#include "RooMinuit.h"

extern bool TOY;
extern double TOYCORR;

class HighNuFCN : public RooAbsReal
{
    public:
        HighNuFCN() : mNBins(0), mBinStep(0) {};
        HighNuFCN(const int inNBin, const int inBinStep);

        virtual TObject* clone(const char* newname = 0) const {
            return new HighNuFCN();
        }

        virtual Double_t evaluate() const;

        std::unique_ptr<RooListProxy> mPulls;

        TH1D* GetHistCombinedNominal() {return this->mHistCombinedNominal.get();};

        TH1D* GetHistGenieNominal() {return this->mHistGenieNominal.get();};
        TH1D* GetHistGenieShift() {return this->mHistGenieShift.get();};
        TH1D* GetHistG4Nominal() {return this->mHistG4Nominal.get();};
        TH1D* GetHistG4Shift() {return this->mHistG4Shift.get();};
        TH1D* GetHistSampleResult() {return this->mSampleResult.get();};

        TMatrixD* GetCorrelationMatrix() {return this->mCorrelationMatrix.get();};
        TMatrixD* GetToyCorrelationMatrix() {return this->mToyCorrelationMatrix.get();};

        TMatrixD* GetCovarianceMatrix() {return this->mCovarianceMatrix.get();};
        TMatrixD* GetToyCovarianceMatrix() {return this->mToyCovarianceMatrix.get();};

        void SaveHist(std::string_view name);

    private:
        const int mNBins;
        const int mBinStep;

        std::vector<RooRealVar*> mParVec;
        std::unique_ptr<TH1D> mHistGenieNominal;
        void SetHistGenieNominal(std::string inputFile);

        std::unique_ptr<TH1D> mHistGenieShift;
        void SetHistGenieShift(std::string inputFile);
        void SetHistGenieNominalError();

        std::unique_ptr<TH1D> mHistG4Nominal;
        void SetHistG4Nominal(std::string inputFile);

        std::unique_ptr<TH1D> mHistG4Shift;
        void SetHistG4Shift(std::string inputFile);
        void SetHistG4NominalError();

        std::unique_ptr<TH1D> mHistCombinedNominal;
        void SetHistCombinedError();

        std::unique_ptr<TH1D> mSampleResult;
        void SamplingHistograms(int inSamplingNumber);

        std::vector<TH1D> mSampledHist = {};
        TH1D SamplingEachHistogram();

        std::unique_ptr<TMatrixD> mCorrelationMatrix;
        void SetCorrelationMatrix();

        std::unique_ptr<TMatrixD> mToyCorrelationMatrix;
        void SetToyCorrelationMatrix();

        std::unique_ptr<TMatrixD> mCovarianceMatrix;
        void SetCovarianceMatrix();

        std::unique_ptr<TMatrixD> mToyCovarianceMatrix;
        void SetToyCovarianceMatrix();

        double PredictionAndData() const ;
        double PenaltyForParameters() const;
        double PenaltyForParametersToyModel() const;

        std::unique_ptr<TH1D> mN1Data;
        std::unique_ptr<TH1D> mN2Data;
        std::unique_ptr<TH1D> mN3Data;
};

#endif
