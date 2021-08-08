#include "LowNuFCN.hxx"

#include "TCanvas.h"

double bkgUnc = 0.32;
double bkgWeight = 1;
float kScaling = 2.4E+6 * 0.12; //1year
//-----------------------------------------------------------------------------
LowNuFCN::LowNuFCN(
        int numPars, double inError, std::string inputFluxSystematic,
        std::string inputDataFile)
    : mNumberOfParameters(numPars)
    , mError(inError)
    , mFluxSystematicFileName(inputFluxSystematic) 
    , mDataFileName(inputDataFile) {
    mPulls = std::make_unique<RooListProxy>("mPulls", "mPulls", this);

    for (int i = 0; i < this->GetNumberOfParameters(); ++i) {
        RooRealVar* tempPar = new RooRealVar(Form("flux systematic %d", i+1), 
                Form("par%d", i+1), 0);
        tempPar->setConstant(false);
        tempPar->setError(1);
        //auto tempPar = std::make_unique<RooRealVar>(Form("flux systematic %d", i), 
                //Form("par%d", i+1), 0);
        mFluxPars.push_back(tempPar);
        mPulls->add(*tempPar);
    }

    mPullsbkg = std::make_unique<RooListProxy>("mPullsbkg", "mPullsbkg", this);
    RooRealVar* parBkg = new RooRealVar("background", "background", 0);
    parBkg->setError(bkgUnc);
    parBkg->setConstant(false);
    mPullsbkg->add(*parBkg);

    this->addServerList(*mPulls);
    this->addServerList(*mPullsbkg);
    
    this->LoadFluxSystematics(this->mFluxSystematicFileName);
    this->LoadData(this->mDataFileName);
    this->SetCovMatrix();
    this->SetPullCV();
    this->SetPullUnc();
}
//-----------------------------------------------------------------------------
void LowNuFCN::LoadFluxSystematics(std::string fluxSystematicFile) {
    auto fileFluxShifts = std::make_unique<TFile> (fluxSystematicFile.c_str());
    if (!fileFluxShifts->IsOpen())
        throw std::runtime_error("invalid inputFluxSystematic");

    for (int i = 0; i < this->GetNumberOfParameters(); ++i) {
        TH1D* temp = (TH1D*)fileFluxShifts.get()->Get(
                Form("syst%d/ND_numubar_RHC",i));
        this->mFluxSyst.push_back(*temp);
    }
}
//-----------------------------------------------------------------------------
void LowNuFCN::LoadData(std::string dataFileName) {
    dataFile = std::make_unique<TFile>(dataFileName.c_str());
    if (!dataFile->IsOpen()) 
        throw std::runtime_error("invalid input data");
    std::cout << dataFileName << std::endl;
    this->inputTree = (TTree*)dataFile->Get("tree");
    this->SetInputTree();
    this->FillData();
}
//-----------------------------------------------------------------------------
void LowNuFCN::SetInputTree() {
    this->inputTree->SetBranchAddress("recoNeutrinoE"    , &recoNeutrinoE);
    this->inputTree->SetBranchAddress("trueNeutrinoE"    , &trueNeutrinoE);
    this->inputTree->SetBranchAddress("isSignal"         , &isSignal);
}
//-----------------------------------------------------------------------------
void LowNuFCN::FillData() {
    int numBins = this->mFluxSyst[0].GetNbinsX();
    double minimum = this->mFluxSyst[0].GetBinLowEdge(1);
    double maximum = this->mFluxSyst[0].GetBinLowEdge(numBins) 
                     + this->mFluxSyst[0].GetBinWidth(numBins);
    TH1D tempData("data", "data", numBins, minimum, maximum);
    for (int event = 0; event < this->inputTree->GetEntries(); event++) {
        this->inputTree->GetEntry(event);
        if (isSignal) tempData.Fill(recoNeutrinoE/1000.);
        else tempData.Fill(recoNeutrinoE/1000., bkgWeight);
        Event tempEvent;
        tempEvent.recoNeutrinoE = recoNeutrinoE;
        tempEvent.trueNeutrinoE = trueNeutrinoE;
        tempEvent.isSignal = isSignal;
        this->mEvents.push_back(tempEvent);
    }
    this->mData = tempData;
}
//-----------------------------------------------------------------------------
void LowNuFCN::SetPullCV() {
    this->pullCV = std::make_unique<TVectorD>(this->GetNumberOfParameters());
    for(int i = 0; i < this->GetNumberOfParameters(); ++i) {
        (*this->pullCV)[i] = 0;
    }
}
//-----------------------------------------------------------------------------
void LowNuFCN::SetPullUnc() {
    this->pullUnc = std::make_unique<TVectorD>(this->GetNumberOfParameters());
    for(int i = 0; i < this->GetNumberOfParameters(); ++i) {
        (*this->pullUnc)[i] = 1;
    }
}
//-----------------------------------------------------------------------------
TH1D LowNuFCN::GetPrediction() const {
    //using the same number of bins, flux systematic
    int numBins = this->mFluxSyst[0].GetNbinsX();
    double minimum = this->mFluxSyst[0].GetBinLowEdge(1);
    double maximum = this->mFluxSyst[0].GetBinLowEdge(numBins) 
                     + this->mFluxSyst[0].GetBinWidth(numBins);
    TH1D predNuE("prediction", "prediction", numBins, minimum, maximum);

    for (const auto& event : this->mEvents) {
        int trueNeutrinoEBin = event.trueNeutrinoE/500 + 1;

        double weight = this->GetWeight(trueNeutrinoEBin);
        if (event.isSignal) {
            predNuE.Fill(event.recoNeutrinoE/1000., weight);
        } else if (!event.isSignal) {
            predNuE.Fill(event.recoNeutrinoE/1000., bkgWeight * (1 + ((RooRealVar*)mPullsbkg->at(0))->getVal()));
        }
    }

    return predNuE;
}
//-----------------------------------------------------------------------------
double LowNuFCN::GetWeight(int inBin) const {
    double weight = 1;
    for (int tempPar = 0; tempPar < this->GetNumberOfParameters(); 
            ++tempPar) {
        double tempParV = mFluxPars.at(tempPar)->getVal();
        //std::cout << "par[" << tempPar << "]: " << tempParV << std::endl;
        double tempBinV = this->mFluxSyst[tempPar].GetBinContent(inBin);
        weight *= 1 + tempParV * tempBinV;
    }

    return weight;
}
//-----------------------------------------------------------------------------
//data - prediction
//kSacling : anti neutrino CC 0pi event for 3DST per year, CDR
Double_t LowNuFCN::FillEv() const {
    TH1D tempPred = this->GetPrediction();
    TVectorD difference(this->mBins);

    for (int i = 0; i < this->mBins; ++i) { 
        difference[i] = tempPred.GetBinContent(i+1) 
                           - this->mData.GetBinContent(i+1);
        //std::cout << "P - D[" << i << "]: " <<(*difference)[i] << std::endl;
    }

    TMatrixD covMat(*this->GetCovMatrix());
    covMat.Invert();

    TVectorD mulVec(difference);
    mulVec *= covMat;

    Double_t currentResult = mulVec * difference;

    return (Double_t) currentResult; 
}
//-----------------------------------------------------------------------------
void LowNuFCN::SetCovMatrix() {
    //output covariant matrix
    this->mCovMat = std::make_unique<TMatrixD>(this->mBins , this->mBins);

    //only fill diagonal element
    //(i, i) = statistic
    for(int i = 0; i < this->mBins ; ++i) {
        (*mCovMat)(i, i) = this->mData.GetBinContent(i+1);
        if((*mCovMat)(i, i) == 0) 
            (*mCovMat)(i, i) += 0.0000000001;
    }
}
//-----------------------------------------------------------------------------
Double_t LowNuFCN::ExtraPull() const {
    Double_t pullAdd = 0;

    //signal
    for(int i = 0; i < this->GetNumberOfParameters(); ++i) {
        double parVal = mFluxPars.at(i)->getVal();
        pullAdd += std::pow(parVal - 0, 2)/std::pow(1., 2) ;
    }
    //bkg
    pullAdd += std::pow(((RooRealVar*)mPullsbkg->at(0))->getVal() - 0, 2)/std::pow(bkgUnc, 2);

    return pullAdd;
}
//-----------------------------------------------------------------------------
Double_t LowNuFCN::evaluate() const {
    Double_t matPart = this->FillEv();//original FillEv is matPart
    Double_t extraPull = this->ExtraPull();//same variable extraPull
    std::cout << "p - d side: " << matPart << std::endl;
    std::cout << "extraPull: " << extraPull << std::endl;
    std::cout << "------------------------" << std::endl;
    Double_t chi2 = matPart + extraPull; //If needed, add pull terms here.

    return chi2;
}
//-----------------------------------------------------------------------------
TH1D LowNuFCN::GetFittingResult() {
    int numBins = this->mFluxSyst[0].GetNbinsX();
    double minimum = this->mFluxSyst[0].GetBinLowEdge(1);
    double maximum = this->mFluxSyst[0].GetBinLowEdge(numBins) 
                     + this->mFluxSyst[0].GetBinWidth(numBins);
    TH1D tempResult("result", "result", numBins, minimum, maximum);
    for(int i = 0; i < this->mBins ; ++i) {
        double tempBinValue = 0;
        double tempBinDinominator = 0;
        for (int ii = 0; ii < this->GetNumberOfParameters(); ++ii) {
            tempBinValue += std::pow((this->mFluxSyst[ii].GetBinContent(i+1) * ((RooRealVar*)mPulls->at(ii))->getError()), 2);
            tempBinDinominator += std::pow(this->mFluxSyst[ii].GetBinContent(i+1), 2);
        }
        tempResult.SetBinContent(i+1, std::pow(tempBinValue, 0.5)/std::pow(tempBinDinominator, 0.5));
    }

    return tempResult;
}
