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

class FCN : public RooAbsReal
{
    public:
        FCN() : mNBins(0), mBinStep(0) {};
        FCN(const int inNBin, const int inBinStep) 
            : mNBins(inNBin), mBinStep(inBinStep)
        {
            _pulls = new RooListProxy("_pulls","_pulls",this);
            RooRealVar* Par[this->mNBins];
            for (int i = 0; i < this->mNBins; i++) {
                Par[i] = new RooRealVar(Form("par%d", i), 
                                        Form("par%d", i+1), 
                                        1);
                Par[i]->setConstant(false);
                _parlist.add(*(Par[i]));
            }
            _pulls->add(_parlist);
            this->addServerList(*_pulls);

            this->SetHistGenieNominalError();
            this->SetHistG4NominalError();
            this->SetHistCombinedError();
            this->SamplingHistograms(5000);
            this->SetCorrelationMatrix();
            this->SetCovarianceMatrix();
            this->SetToyCorrelationMatrix();
            this->SetToyCovarianceMatrix();
        }

        virtual TObject* clone(const char* newname = 0) const
        {
            return new FCN();
        }

        virtual Double_t evaluate() const;

        ~FCN()
        {
            delete mHistGenieNominal; 
            delete mHistGenieShift;
            delete mHistG4Nominal;
            delete mHistG4Shift; 
            delete mHistCombinedNominal; 
            delete mSampleResult; 
            delete mCorrelationMatrix; 
            delete mToyCorrelationMatrix; 
            delete mCovarianceMatrix; 
            delete mToyCovarianceMatrix; 
            delete mN1Data;
            delete mN2Data;
            delete mN3Data;
        }

        RooArgList _parlist;
        RooListProxy* _pulls;

        TH1D* GetHistCombinedNominal() {return this->mHistCombinedNominal;};

        TH1D* GetHistGenieNominal() {return this->mHistGenieNominal;};
        TH1D* GetHistGenieShift() {return this->mHistGenieShift;};
        TH1D* GetHistG4Nominal() {return this->mHistG4Nominal;};
        TH1D* GetHistG4Shift() {return this->mHistG4Shift;};
        TH1D* GetHistSampleResult() {return this->mSampleResult;};

        TMatrixD* GetCorrelationMatrix() {return this->mCorrelationMatrix;};
        TMatrixD* GetToyCorrelationMatrix() {return this->mToyCorrelationMatrix;};

        TMatrixD* GetCovarianceMatrix() {return this->mCovarianceMatrix;};
        TMatrixD* GetToyCovarianceMatrix() {return this->mToyCovarianceMatrix;};

        void SaveHist(std::string_view name);

    private:
        const int mNBins;
        const int mBinStep;

        TH1D* mHistGenieNominal = nullptr;
        void SetHistGenieNominal();

        TH1D* mHistGenieShift = nullptr;
        void SetHistGenieShift();
        void SetHistGenieNominalError();

        TH1D* mHistG4Nominal = nullptr;
        void SetHistG4Nominal();

        TH1D* mHistG4Shift = nullptr;
        void SetHistG4Shift();
        void SetHistG4NominalError();

        TH1D* mHistCombinedNominal = nullptr;
        void SetHistCombinedError();

        TH1D* mSampleResult = nullptr;
        void SamplingHistograms(int inSamplingNumber);

        std::vector<TH1D> mSampledHist = {};
        TH1D SamplingEachHistogram();

        TMatrixD* mCorrelationMatrix = nullptr;
        void SetCorrelationMatrix();

        TMatrixD* mToyCorrelationMatrix = nullptr;
        void SetToyCorrelationMatrix();

        TMatrixD* mCovarianceMatrix = nullptr;
        void SetCovarianceMatrix();

        TMatrixD* mToyCovarianceMatrix = nullptr;
        void SetToyCovarianceMatrix();

        double PredictionAndData() const ;
        double PenaltyForParameters() const;
        double PenaltyForParametersToyModel() const;

        TH1D* mN1Data = nullptr;
        TH1D* mN2Data = nullptr;
        TH1D* mN3Data = nullptr;
};

#endif
