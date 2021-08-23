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

struct Event
{
    float recoNu = {};
    float reWeight[100] = {};
    float trueNu = {};      
    int numberOfFSNeutron = {};
};

class HighNuFCN : public RooAbsReal
{

    public:
        HighNuFCN() : mNBins(0), mBinStep(0) {};
        HighNuFCN(const int inNBin, const int inBinStep, bool N1HighNu, bool N1NonNeutron, bool N2, bool N3);

        virtual TObject* clone(const char* newname = 0) const {
            return new HighNuFCN();
        }

        virtual Double_t evaluate() const;

        std::unique_ptr<RooListProxy> mPulls;
        std::unique_ptr<RooListProxy> mPullsN1HighNu;
        std::unique_ptr<RooListProxy> mPullsN1NonNeutron;
        std::unique_ptr<RooListProxy> mPullsN2;
        std::unique_ptr<RooListProxy> mPullsN3;

        TH1D* GetHistCombinedNominal() const {return this->mHistCombinedNominal.get();};

        TH1D* GetHistGenieNominal() const {return this->mHistNominalAll.get();};
        TH1D* GetHistGenieShift() const {return this->mHistShiftAll.get();};
        TH1D* GetHistG4Nominal() const {return this->mHistG4Nominal.get();};
        TH1D* GetHistG4Shift() const {return this->mHistG4Shift.get();};
        TH1D* GetHistSampleResult() const {return this->mSampleResult.get();};

        TMatrixD* GetCorrelationMatrix() const {return this->mCorrelationMatrixAll.get();};
        TMatrixD* GetCorrelationMatrixN1HighNu() const {return this->mCorrelationMatrixN1HighNu.get();};
        TMatrixD* GetCorrelationMatrixN2() const {return this->mCorrelationMatrixN2.get();};
        TMatrixD* GetCorrelationMatrixN3() const {return this->mCorrelationMatrixN3.get();};
        TMatrixD* GetToyCorrelationMatrix() const {return this->mToyCorrelationMatrix.get();};

        TMatrixD* GetCovarianceMatrix() const {return this->mCovarianceMatrixAll.get();};
        TMatrixD* GetCovarianceMatrixN1HighNu() const {return this->mCovarianceMatrixN1HighNu.get();};
        TMatrixD* GetCovarianceMatrixN2() const {return this->mCovarianceMatrixN2.get();};
        TMatrixD* GetCovarianceMatrixN3() const {return this->mCovarianceMatrixN3.get();};
        TMatrixD* GetToyCovarianceMatrix() const {return this->mToyCovarianceMatrix.get();};

        void SaveHist(std::string_view name) const;

        std::vector<RooRealVar*> GetParVec() const {return this->mParVec;};

        mutable std::vector<double> mParV;
        mutable std::vector<double> mChi2;

        void SetParVec(std::vector<double> inVec);

        TMatrixD MyInvert(const TMatrixD& inMatrix) const;

        bool N1HighNu;
        bool N1NonNeutron;
        bool N2;
        bool N3;

        std::vector<Event> mEvents;
        void SetEvents(std::string inputFile);
    private:
        const int mNBins;
        const int mBinStep;

        void InitializeHistograms();

        std::unique_ptr<TH1D> mHistNominalAll;
        std::unique_ptr<TH1D> mHistShiftAll;

        std::unique_ptr<TH1D> mHistNominalN1NonNeutron;
        std::unique_ptr<TH1D> mHistNominalN1HighNu;
        std::unique_ptr<TH1D> mHistNominalN2;
        std::unique_ptr<TH1D> mHistNominalN3;

        std::unique_ptr<TH1D> mHistShiftN1NonNeutron;
        std::unique_ptr<TH1D> mHistShiftN1HighNu;
        std::unique_ptr<TH1D> mHistShiftN2;
        std::unique_ptr<TH1D> mHistShiftN3;

        std::unique_ptr<TH1D> FillHist(std::string inputFile);

        std::vector<RooRealVar*> mParVec;
        std::vector<RooRealVar*> mParVecN1HighNu;
        std::vector<RooRealVar*> mParVecN1NonNeutron;
        std::vector<RooRealVar*> mParVecN2;
        std::vector<RooRealVar*> mParVecN3;
        void SetHistGenieNominal(std::string inputFile);

        void SetHistGenieShift(std::string inputFile);
        void SetHistograms();

        std::unique_ptr<TH1D> mHistG4Nominal;
        void SetHistG4Nominal(std::string inputFile);

        std::unique_ptr<TH1D> mHistG4Shift;
        void SetHistG4Shift(std::string inputFile);
        void SetHistG4NominalError();

        std::unique_ptr<TH1D> mHistCombinedNominal;
        void SetHistCombinedError();

        std::unique_ptr<TH1D> mSampleResult;
        std::vector<TH1D> SamplingHistograms(int inSamplingNumber, const TH1D& inNominal, const TH1D& inShift);

        std::vector<TH1D> mSampledHist = {};
        TH1D SamplingEachHistogram(const TH1D& inNominal, const TH1D& inShift);

        std::unique_ptr<TMatrixD> mCorrelationMatrixAll;
        std::unique_ptr<TMatrixD> mCorrelationMatrixN1HighNu;
        std::unique_ptr<TMatrixD> mCorrelationMatrixN1NonNeutron;
        std::unique_ptr<TMatrixD> mCorrelationMatrixN2;
        std::unique_ptr<TMatrixD> mCorrelationMatrixN3;
        std::unique_ptr<TMatrixD> SetCorrelationMatrix(const TH1D& inNominal, const TH1D& inShift);

        std::unique_ptr<TMatrixD> mToyCorrelationMatrix;
        void SetToyCorrelationMatrix();

        std::unique_ptr<TMatrixD> mCovarianceMatrixAll;
        std::unique_ptr<TMatrixD> mCovarianceMatrixN1HighNu;
        std::unique_ptr<TMatrixD> mCovarianceMatrixN1NonNeutron;
        std::unique_ptr<TMatrixD> mCovarianceMatrixN2;
        std::unique_ptr<TMatrixD> mCovarianceMatrixN3;
        std::unique_ptr<TMatrixD> SetCovarianceMatrix(const TMatrixD& inCorrMat);

        std::unique_ptr<TMatrixD> mToyCovarianceMatrix;
        void SetToyCovarianceMatrix();

        std::unique_ptr<TMatrixD> mTestCov;
        void SetTestCov();
        double TestChi2() const;

        double PredictionMinusData(const TH1D& inNominal) const ;
        double Correlation() const;
        double ExtraPenaltyForParameters() const;
        double CorrelationToyModel() const;

        double DoMatrixCal(const TVectorD& inVec, TMatrixD inMat) const;
};

#endif
