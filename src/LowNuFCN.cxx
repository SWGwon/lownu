#include "LowNuFCN.hxx"

#include "TCanvas.h"

float kScaling = 2.4E+6 * 0.25; //1year
const float kCorelation = 0.0;
//-----------------------------------------------------------------------------
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
        par[i] = new RooRealVar(Form("flux systematic %d", i), 
                Form("par%d", i+1), 0);
        par[i]->setConstant(false);
        _parlist.add(*(par[i]));
    }

    RooRealVar* parBkg = new RooRealVar("background", 
            Form("par%d", this->GetNumberOfParameters()+1), 0);
    parBkg->setConstant(false);
    _parlist.add(*(parBkg));

    _pulls->add(_parlist);
    this->addServerList(*_pulls);
    
    this->LoadFluxSystematics(this->mFluxSystematicFileName);
    this->LoadData(this->mDataFileName);
    this->SetPullCV();
    this->SetPullUnc();
}
//-----------------------------------------------------------------------------
void LowNuFCN::LoadFluxSystematics(std::string fluxSystematicFile) {
    auto fileFluxShifts = std::make_unique<TFile> (fluxSystematicFile.c_str());
    if (!fileFluxShifts->IsOpen())
        throw std::runtime_error("invalid inputFluxSystematic");

    for (int i = 0; i < this->GetNumberOfParameters(); i++) {
        TH1D* temp = (TH1D*)fileFluxShifts.get()->Get(
                Form("syst%d/ND_numubar_RHC",i));
        this->mFluxSyst.push_back(*temp);
    }
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
    double maximum = this->mFluxSyst[0].GetBinLowEdge(numBins) 
                     + this->mFluxSyst[0].GetBinWidth(numBins);
    TH1D tempData("data", "data", numBins, minimum, maximum);
    for (int event = 0; event < this->inputTree->GetEntries(); event++) {
        this->inputTree->GetEntry(event);
        tempData.Fill(recoNeutrinoE/1000.);
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
TH1D LowNuFCN::GetPrediction(RooListProxy* _pulls) const {
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
            predNuE.Fill(event.recoNeutrinoE/1000., (1 + ((RooAbsReal*)_pulls->at(this->GetNumberOfParameters()))->getVal()));
        }
    }

    return predNuE;
}
//-----------------------------------------------------------------------------
double LowNuFCN::GetWeight(int inBin) const {
    double weight = 1;
    for (int tempPar = 0; tempPar < this->GetNumberOfParameters(); 
            ++tempPar) {
        double tempParV = ((RooAbsReal*)_pulls->at(tempPar))->getVal();
        //std::cout << "par[" << tempPar << "]: " << tempParV << std::endl;
        double tempBinV = this->mFluxSyst[tempPar].GetBinContent(inBin);
        weight *= 1 + tempParV * tempBinV;
    }
    return weight;
}
//-----------------------------------------------------------------------------
Double_t LowNuFCN::FillEv(RooListProxy* _pulls) const {
    TH1D tempPred = this->GetPrediction(_pulls);

    TVectorD* difference = new TVectorD(this->mBins);

    //data - prediciton
    //kSacling : anti neutrino CC 0pi event for 3DST per year, CDR
    for (Int_t i = 0; i < this->mBins; i++) { 
        (*difference)[i] =  tempPred.GetBinContent(i+1) - this->mData.GetBinContent(i+1);
        //std::cout << "P - D[" << i << "]: " <<(*difference)[i] << std::endl;
    }

    TMatrixD* covMat = this->PrepareCovMatrix();
    covMat->Invert();

    TVectorD mulVec(*difference);
    mulVec *= (*covMat);

    Double_t currentResult = mulVec*(*difference);

    return (Double_t) currentResult; 
}
//-----------------------------------------------------------------------------
TMatrixD* LowNuFCN::PrepareCovMatrix() const {
    //output covariant matrix
    TMatrixD* outMat = new TMatrixD(this->mBins , this->mBins);

    //only fill diagonal element
    //(i, i) = statistic
    for(Int_t i = 0; i < this->mBins ; i++) {
        (*outMat)(i, i) = this->mData.GetBinContent(i+1);
        if((*outMat)(i, i) == 0) 
            (*outMat)(i, i) += 0.0000000001;
    }
    for(Int_t i = 0; i < this->mBins ; i++) {
        for(Int_t j = 0; j < this->mBins ; j++) {
            if (i == j) {
                continue;
            }
        }
    }

    return outMat ;
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
    *covMat = *(this->PrepareCovMatrix2());
    covMat->Invert();

    TVectorD mulVec(*difference);
    mulVec *= (*covMat);

    Double_t currentResult = TMath::Abs(mulVec*(*difference));

    return (Double_t) currentResult; 
}
//-----------------------------------------------------------------------------
TMatrixD* LowNuFCN::PrepareCovMatrix2() const {
    //output covariant matrix
    TMatrixD* outMat = new TMatrixD(this->mBins , this->mBins);

    //only fill diagonal element
    //(i, i) = cross section uncertainty^2
    //{
    for(Int_t i = 0; i < this->mBins ; i++) {
        (*outMat)(i, i) = std::pow(mError, 2);// + std::pow(detectorSmeartingHighNu, 2);
    }
    //off diagonal : correlation
    for(Int_t i = 0; i < this->mBins ; i++) {
        for(Int_t j = 0; j < this->mBins ; j++) {
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
Double_t LowNuFCN::ExtraPull(RooListProxy* _pulls) const {
    Double_t pullAdd = 0;

    for(Int_t i = 0; i < this->GetNumberOfParameters(); i++) {
        pullAdd += std::pow(((RooAbsReal*)_pulls->at(i))->getVal() - 0, 2)/std::pow(1., 2) ;
        //pullAdd += std::pow(((RooAbsReal*)_pulls->at(i))->getVal() - (*pullCV)[i], 2)/TMath::Power((*pullUnc)[i], 2) ;
    }
    //bkg
    double bkgUnc = 1.;
    pullAdd += std::pow(((RooAbsReal*)_pulls->at(this->GetNumberOfParameters()))->getVal() - 0, 2)/std::pow(bkgUnc, 2);

    return pullAdd;
}
//-----------------------------------------------------------------------------
Double_t LowNuFCN::evaluate() const {
    Double_t matPart = this->FillEv(_pulls);//original FillEv is matPart
    //Double_t energyPart = this->FillEv2(_pulls);//(e_i - CV)^T * ( ) * (e_i - CV)
    Double_t extraPull = this->ExtraPull(_pulls);//same variable extraPull
    //std::cout << "p - d side: " << matPart << std::endl;
    //std::cout << "e - CV side: " << energyPart << std::endl;
    //std::cout << "extraPull: " << extraPull << std::endl;
    //std::cout << "------------------------" << std::endl;
    Double_t chi2 = matPart + extraPull; //If needed, add pull terms here.

    return chi2;
}
//-----------------------------------------------------------------------------
LowNuFCN::~LowNuFCN() {
    delete dataFile;
    delete _pulls;
}
