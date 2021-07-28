#include "LowNuFCN.hxx"

#include "TCanvas.h"

float kScaling = 2.4E+6 * 0.25; //1year
const float kCorelation = 0.0;

LowNuFCN::LowNuFCN(
        int numPars, double inError, std::string inputFluxSystematic,
        std::string inputDataFile)
    : mNumberOfParameters(numPars)
    , mError(inError)
    , mFluxSystematicFileName(inputFluxSystematic) 
    , mDataFileName(inputDataFile) {
    _pulls = new RooListProxy("_pulls", "_pulls", this);
    
    //Par[this->GetNumberOfParameters()+1] is for bkg
    //E' = a + b*E_reco
    //Par[this->GetNumberOfParameters()+2] : a
    //Par[this->GetNumberOfParameters()+3] : b
    //Par[this->GetNumberOfParameters()+ 4~20] : each energy bin
    
    RooRealVar* par[this->GetNumberOfParameters()];
    for (int i = 0; i < this->GetNumberOfParameters(); i++) {
        par[i] = new RooRealVar(Form("flux systematic %d", i), Form("par%d", i+1), 0);
        par[i]->setConstant(false);
        _parlist.add(*(par[i]));
    }

    RooRealVar* parBkg = new RooRealVar("background", Form("par%d", this->GetNumberOfParameters()+1), 0);
    parBkg->setConstant(false);
    _parlist.add(*(parBkg));

    _pulls->add(_parlist);
    this->addServerList(*_pulls);
    
    this->LoadFluxSystematics(this->mFluxSystematicFileName);
    this->LoadData(this->mDataFileName);
    this->SetPullCV();
    this->SetPullUnc();

    //this->PreparePrediction(_pulls);
}
//-----------------------------------------------------------------------------
void LowNuFCN::LoadFluxSystematics(std::string fluxSystematicFile) {
    auto fileFluxShifts = std::make_unique<TFile> (fluxSystematicFile.c_str());
    if (!fileFluxShifts->IsOpen())
        throw std::runtime_error("invalid inputFluxSystematic");

    for (int i = 0; i < this->GetNumberOfParameters(); i++) {
        TH1D* temp = (TH1D*)fileFluxShifts.get()->Get(Form("syst%d/ND_numubar_RHC",i));
        this->mFluxSyst.push_back(*temp);
    }

    //TCanvas c;
    //c.Divide(5,5);
    //for (int i = 0; i < this->GetNumberOfParameters(); i++) {
    //    std::cout << i << std::endl;
    //    c.cd(i+1);
    //    mFluxSyst.at(i).Draw();
    //}
    //c.SaveAs("asd.pdf");
}
//-----------------------------------------------------------------------------
void LowNuFCN::LoadData(std::string dataFileName) {
    dataFile = new TFile(dataFileName.c_str());
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
    double maximum = this->mFluxSyst[0].GetBinLowEdge(numBins) + this->mFluxSyst[0].GetBinWidth(numBins);
    this->mData = new TH1D("data", "data", numBins, minimum, maximum);
    for (int event = 0; event < this->inputTree->GetEntries(); event++) {
        this->inputTree->GetEntry(event);
        this->mData->Fill(recoNeutrinoE);
        Event tempEvent;
        tempEvent.recoNeutrinoE = recoNeutrinoE;
        tempEvent.trueNeutrinoE = trueNeutrinoE;
        tempEvent.isSignal = isSignal;
        this->mEvents.push_back(tempEvent);
    }
}
//-----------------------------------------------------------------------------
void LowNuFCN::SetPullCV() {
    this->pullCV = new TVectorD(this->GetNumberOfParameters());
    for(int i = 0; i < this->GetNumberOfParameters(); ++i) {
        (*this->pullCV)[i] = 0;
    }
}
//-----------------------------------------------------------------------------
void LowNuFCN::SetPullUnc() {
    this->pullUnc = new TVectorD(this->GetNumberOfParameters());
    for(int i = 0; i < this->GetNumberOfParameters(); ++i) {
        (*this->pullUnc)[i] = 1;
    }
}
//-----------------------------------------------------------------------------
TH1D LowNuFCN::PreparePrediction(RooListProxy* _pulls) const {
    //using the same number of bins, flux systematic
    int numBins = this->mFluxSyst[0].GetNbinsX();
    double minimum = this->mFluxSyst[0].GetBinLowEdge(1);
    double maximum = this->mFluxSyst[0].GetBinLowEdge(numBins) + this->mFluxSyst[0].GetBinWidth(numBins);
    TH1D predNuE("", "", numBins, minimum, maximum);

    //for (int event = 0; event < this->inputTree->GetEntries(); event++) {
    //    this->inputTree->GetEntry(event);
    //    std::cout << event << std::endl;

    //    int trueNeutrinoEBin = trueNeutrinoE/500 + 1;

    //    double weight = 1;
    //    if (isSignal) {
    //        for (int tempPar = 0; tempPar < this->GetNumberOfParameters(); tempPar++) {
    //            weight *= 1 + ((RooAbsReal*)_pulls->at(tempPar))->getVal() * this->mFluxSyst[tempPar].GetBinContent(trueNeutrinoEBin);
    //        }
    //        predNuE.Fill(recoNeutrinoE/1000., weight);
    //    }

    //    //bkg
    //    if (!isSignal) {
    //        //use 100% (1) error for background
    //        predNuE.Fill(recoNeutrinoE/1000., (1 + ((RooAbsReal*)_pulls->at(this->GetNumberOfParameters()))->getVal()));
    //    }
    //}
    for (const auto& event : this->mEvents) {
        int trueNeutrinoEBin = event.trueNeutrinoE/500 + 1;

        double weight = 1;
        if (event.isSignal) {
            for (int tempPar = 0; tempPar < this->GetNumberOfParameters(); tempPar++) {
                weight *= 1 + ((RooAbsReal*)_pulls->at(tempPar))->getVal() * this->mFluxSyst[tempPar].GetBinContent(trueNeutrinoEBin);
            }
            predNuE.Fill(event.recoNeutrinoE/1000., weight);
        }

        //bkg
        if (!event.isSignal) {
            //use 100% (1) error for background
            predNuE.Fill(event.recoNeutrinoE/1000., (1 + ((RooAbsReal*)_pulls->at(this->GetNumberOfParameters()))->getVal()));
        }
    }
    return predNuE;
}
//-----------------------------------------------------------------------------
TMatrixD* LowNuFCN::PrepareCovMatrix(Int_t nBins, TH1D pred) const {
    //output covariant matrix
    TMatrixD* outMat = new TMatrixD(nBins , nBins);

    //only fill diagonal element
    //(i, i) = statistic
    for(Int_t i = 0; i < nBins ; i++) {
        (*outMat)(i, i) = this->mData->GetEntries();
        if((*outMat)(i, i) == 0) 
            (*outMat)(i, i) += 0.0000000001;
    }
    for(Int_t i = 0; i < nBins ; i++) {
        for(Int_t j = 0; j < nBins ; j++) {
            if (i == j) {
                continue;
            }
        }
    }

    return outMat ;
}
//-----------------------------------------------------------------------------
TMatrixD* LowNuFCN::PrepareCovMatrix2(Int_t nBins) const {
    //output covariant matrix
    TMatrixD* outMat = new TMatrixD(nBins , nBins);

    //only fill diagonal element
    //(i, i) = cross section uncertainty^2
    //{
    for(Int_t i = 0; i < nBins ; i++) {
        (*outMat)(i, i) = std::pow(mError, 2);// + std::pow(detectorSmeartingHighNu, 2);
    }
    //off diagonal : correlation
    for(Int_t i = 0; i < nBins ; i++) {
        for(Int_t j = 0; j < nBins ; j++) {
            if (i == j) {
                continue;
            }
            (*outMat)(i, j) = kCorelation * std::pow((*outMat)(i,i), 0.5) * std::pow((*outMat)(j,j), 0.5);
        }
    }
    //}

    return outMat ;
}
//-----------------------------------------------------------------------------
Double_t LowNuFCN::FillEv(RooListProxy* _pulls) const {
    TH1D tempPred = this->PreparePrediction(_pulls);

    TVectorD* difference = new TVectorD(this->mBins);

    //data - prediciton
    //kSacling : anti neutrino CC 0pi event for 3DST per year, CDR
    for (Int_t i = 0; i < this->mBins; i++) { 
        (*difference)[i] = this->mData->GetBinContent(i+1) - tempPred.GetBinContent(i+1);
        //std::cout << "P - D[" << i << "]: " <<(*difference)[i] << std::endl;
    }

    TMatrixD* covMat = new TMatrixD(this->mBins, this->mBins);
    *covMat = *(this->PrepareCovMatrix(this->mBins, tempPred));
    covMat->Invert();

    TVectorD mulVec(*difference);
    mulVec *= (*covMat);

    Double_t currentResult = TMath::Abs(mulVec*(*difference));

    return (Double_t) currentResult; 
}
//-----------------------------------------------------------------------------
Double_t LowNuFCN::FillEv2(RooListProxy* _pulls) const {
    TVectorD* e_i = new TVectorD(this->mBins);
    TVectorD* centralValue = new TVectorD(this->mBins);
    TVectorD* difference = new TVectorD(this->mBins);

    for (Int_t i = 0; i < this->mBins; i++) {	 
        (*e_i)[i] = ((RooAbsReal*)_pulls->at(this->GetNumberOfParameters()+2+i+1))->getVal();
    }
    for (Int_t i = 0; i < this->mBins; i++) {
        (*centralValue)[i] = 0;
    }

    //e_i - centralValue
    for (Int_t i = 0; i < this->mBins; i++) { 
        (*difference)[i] = (*e_i)[i] - (*centralValue)[i];
    }

    TMatrixD* covMat = new TMatrixD(this->mBins, this->mBins);
    *covMat = *(this->PrepareCovMatrix2(this->mBins));
    covMat->Invert();

    TVectorD mulVec(*difference);
    mulVec *= (*covMat);

    Double_t currentResult = TMath::Abs(mulVec*(*difference));

    return (Double_t) currentResult; 
}
//-----------------------------------------------------------------------------
Double_t LowNuFCN::ExtraPull(RooListProxy* _pulls) const {
    Double_t pullAdd = 0;

    for(Int_t i = 0; i < this->GetNumberOfParameters(); i++) {
        pullAdd += std::pow(((RooAbsReal*)_pulls->at(i))->getVal() - (*pullCV)[i], 2)/TMath::Power((*pullUnc)[i], 2) ;
    }
    //bkg
    double bkgUnc = 1;
    pullAdd += std::pow(((RooAbsReal*)_pulls->at(this->GetNumberOfParameters()))->getVal() - 0, 2)/TMath::Power(bkgUnc, 2);

    return pullAdd;
}
//-----------------------------------------------------------------------------
Double_t LowNuFCN::evaluate() const {
    Double_t matPart = this->FillEv(_pulls);//original FillEv is matPart
    //Double_t energyPart = this->FillEv2(_pulls);//(e_i - CV)^T * ( ) * (e_i - CV)
    Double_t extraPull = this->ExtraPull(_pulls);//same variable extraPull
    std::cout << "p - d side: " << matPart << std::endl;
    //std::cout << "e - CV side: " << energyPart << std::endl;
    std::cout << "extraPull: " << extraPull << std::endl;
    std::cout << "------------------------" << std::endl;
    Double_t chi2 = matPart + extraPull; //If needed, add pull terms here.

    return chi2;
}
//-----------------------------------------------------------------------------
LowNuFCN::~LowNuFCN() {
    delete dataFile;
    delete inputTree;
}
