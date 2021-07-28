#include "Utils.hxx"
#include "HighNuFCN.hxx"

bool combinHighNu = false;
bool TOY = false;
double TOYCORR = 0;

/////////////////////////////prepare histogram//////////////////////////////////
void FCN::SetHistGenieNominal() {
    TH1D* tempHistNominal = new TH1D("HistGenieNominal", 
            "HistGenieNominal;reco #nu;Normalized fraction", 
            this->mNBins, 0, this->mNBins * this->mBinStep);
    TH1D*  tempHistN1Data = new TH1D("HistN1Data", "HistN1Data", 
            this->mNBins, 0, this->mNBins * this->mBinStep);

    TFile* tempFile = nullptr;
    try {
        tempFile = new TFile ("analysis_output_G1801a00000AfterCut.root");
    } catch (...) {
        std::runtime_error("invalid nominal genie input file");
    }

    TTree* tempTree = (TTree*)tempFile->Get("tree");
    float recoNu;          tempTree->SetBranchAddress("recoNu", &recoNu);
    int numberOfFSNeutron; tempTree->SetBranchAddress("numberOfFSNeutron", 
            &numberOfFSNeutron);

    std::cout << __func__ << std::endl;
    for (int i = 0; i < tempTree->GetEntries(); ++i) {
        PrintProgress(i, tempTree->GetEntries());
        tempTree->GetEntry(i);
        if (numberOfFSNeutron == 1) continue;
        if (combinHighNu) {
            if (300 < recoNu && recoNu < 550) recoNu = 325;
            if (550 < recoNu && recoNu < 800) recoNu = 375;
        }
        if (numberOfFSNeutron > 0)
        tempHistNominal->Fill(recoNu);
        //if (numberOfFSNeutron == 1)
        tempHistN1Data->Fill(recoNu);
    }

    this->mHistGenieNominal = tempHistNominal;
    this->mN1Data = tempHistN1Data;

    TCanvas c;
    tempHistNominal->DrawNormalized();
    c.SaveAs("slides.pdf");

    delete tempFile;
}
//------------------------------------------------------------------------------
void FCN::SetHistGenieShift() {
    std::cout << __func__ << std::endl;
    TH1D* tempHistShift = new TH1D("HistGenieShift", 
            "HistGenieShift;reco #nu;Normalized fraction", 
            this->mNBins, 0, this->mNBins * this->mBinStep);

    TFile* tempFile = nullptr;
    try {
        tempFile = new TFile("analysis_output_G1802a00000AfterCut.root");
    } catch (...) {
        std::runtime_error("invalid nominal genie input file");
    }
    TTree* tempTree = (TTree*)tempFile->Get("tree");
    float recoNu;          tempTree->SetBranchAddress("recoNu", &recoNu);
    int numberOfFSNeutron; tempTree->SetBranchAddress("numberOfFSNeutron", 
                                                      &numberOfFSNeutron);

    for (int i = 0; i < tempTree->GetEntries(); ++i) {
        PrintProgress(i, tempTree->GetEntries());
        tempTree->GetEntry(i);
        if (numberOfFSNeutron == 1) continue;
        if (combinHighNu) {
            if (300 < recoNu && recoNu < 550) recoNu = 325;
            if (550 < recoNu && recoNu < 800) recoNu = 375;
        }
        if (numberOfFSNeutron > 0)
        tempHistShift->Fill(recoNu);
    }
        
    this->mHistGenieShift = tempHistShift;

    delete tempFile;
}
//------------------------------------------------------------------------------
void FCN::SetHistGenieNominalError() {

    SetHistGenieNominal();
    SetHistGenieShift();

    this->GetHistGenieNominal()->Scale(
            1 / this->GetHistGenieNominal()->Integral(), "nosw2");
    this->GetHistGenieShift()->Scale(
            1 / this->GetHistGenieShift()->Integral(), "nosw2");

    for (int i = 0; i < this->mNBins; ++i) {
        double tempDifference = 0;
        tempDifference = this->GetHistGenieNominal()->GetBinContent(i+1) 
                       - this->GetHistGenieShift()->GetBinContent(i+1);

        this->GetHistGenieNominal()->SetBinError(i+1, tempDifference);
    }

}
//------------------------------------------------------------------------------
void FCN::SetHistG4Nominal() {
    TH1D* tempHistNominal = new TH1D("HistG4Nominal", 
            "HistG4Nominal;reco #nu;Normalized fraction", 
            this->mNBins, 0, this->mNBins * this->mBinStep);

    TFile* tempFile = new TFile("G4comparisonDefaultAfterCut.root");
    TTree* tempTree = (TTree*)tempFile->Get("tree");
    float recoNu;          tempTree->SetBranchAddress("recoNu", &recoNu);
    int numberOfFSNeutron; tempTree->SetBranchAddress("numberOfFSNeutron", 
                                                      &numberOfFSNeutron);

    std::cout << __func__ << std::endl;
    for (int i = 0; i < tempTree->GetEntries(); ++i) {
        PrintProgress(i, tempTree->GetEntries());
        tempTree->GetEntry(i);
        if (numberOfFSNeutron == 1) continue;
        if (combinHighNu) {
            if (300 < recoNu && recoNu < 550)
                recoNu = 325;
            if (550 < recoNu && recoNu < 800)
                recoNu = 375;
        }
        if (numberOfFSNeutron > 0)
        tempHistNominal->Fill(recoNu);
        //if (numberOfFSNeutron == 1)
    }

    this->mHistG4Nominal = tempHistNominal;

    delete tempFile;
}
//------------------------------------------------------------------------------
void FCN::SetHistG4Shift() {
    TH1D* tempHistShift = new TH1D("HistG4Shift", 
            "HistG4Shift;reco #nu;Normalized fraction", 
            this->mNBins, 0, this->mNBins * this->mBinStep);

    TFile* tempFile = nullptr;
    try {
        tempFile = new TFile("G4comparisoninclAfterCut.root");
    } catch (...) {
        std::runtime_error("invalid nominal genie input file");
    }
    TTree* tempTree = (TTree*)tempFile->Get("tree");
    float recoNu;          tempTree->SetBranchAddress("recoNu", &recoNu);
    int numberOfFSNeutron; tempTree->SetBranchAddress("numberOfFSNeutron", 
            &numberOfFSNeutron);

    std::cout << __func__ << std::endl;
    for (int i = 0; i < tempTree->GetEntries(); ++i) {
        PrintProgress(i, tempTree->GetEntries());
        tempTree->GetEntry(i);
        if (numberOfFSNeutron == 1) continue;
        if (combinHighNu) {
            if (300 < recoNu && recoNu < 550) recoNu = 325;
            if (550 < recoNu && recoNu < 800) recoNu = 375;
        }
        if (numberOfFSNeutron > 0)
        tempHistShift->Fill(recoNu);
    }
        
    this->mHistG4Shift = tempHistShift;

    delete tempFile;
}
//------------------------------------------------------------------------------
void FCN::SetHistG4NominalError() {

    SetHistG4Nominal();
    SetHistG4Shift();

    this->GetHistG4Nominal()->Scale(
            1 / this->GetHistG4Nominal()->Integral(), "nosw2");
    this->GetHistG4Shift()->Scale(
            1 / this->GetHistG4Shift()->Integral(), "nosw2");

    for (int i = 0; i < this->mNBins; ++i) {
        double tempDifference = 0;
        tempDifference = this->GetHistG4Nominal()->GetBinContent(i+1) 
                       - this->GetHistG4Shift()->GetBinContent(i+1);

        this->GetHistG4Nominal()->SetBinError(i+1, tempDifference);
    }

}
//------------------------------------------------------------------------------
void FCN::SetHistCombinedError() {
    this->mHistCombinedNominal = new TH1D("combinedNominal", 
            "combinedNominal;reco #nu;Normalized fraction", 
            this->mNBins, 0, this->mNBins * this->mBinStep);

    for (int i = 0; i < this->mNBins; ++i) {
        double tempGenieError = this->GetHistGenieNominal()->GetBinError(i+1) 
                              / this->GetHistGenieNominal()->GetBinContent(i+1);
        double tempG4Error = this->GetHistG4Nominal()->GetBinError(i+1) 
                           / this->GetHistG4Nominal()->GetBinContent(i+1);

        this->GetHistCombinedNominal()->SetBinContent(i+1, 
                this->GetHistGenieNominal()->GetBinContent(i+1));
        this->GetHistCombinedNominal()->SetBinError(i+1, 
                std::pow(std::pow(tempGenieError, 2) 
              + std::pow(tempG4Error, 2), 0.5) 
                * this->GetHistCombinedNominal()->GetBinContent(i+1));
    }
}
//------------------------------------------------------------------------------
////////////////////////////////sampling histogram//////////////////////////////
TH1D FCN::SamplingEachHistogram() {
    TH1D tempSample("","tempSample", this->mNBins, 0, 
                    this->mBinStep * this->mNBins);

    TRandom random;
    random.SetSeed(0);
    for (int i = 0; i < this->mNBins; ++i) {
        tempSample.SetBinContent(i+1, std::max((double)0, 
                random.Gaus(this->GetHistCombinedNominal()->GetBinContent(i+1),
                            this->GetHistCombinedNominal()->GetBinError(i+1)))
                );
    }
    tempSample.Scale(
            this->GetHistCombinedNominal()->Integral() / tempSample.Integral(),
            "nosw2");

    return tempSample;
}
//------------------------------------------------------------------------------
void FCN::SamplingHistograms(int inSamplingNumber) {
    std::clock_t START = std::clock();
    std::cout << __func__;
    std::cout << ", inSamplingNumber: " << inSamplingNumber << std::endl;

    this->mSampleResult = new TH1D("","sampling result", 
            this->mNBins, 0, this->mBinStep * this->mNBins);
    
    for (int i = 0; i < inSamplingNumber; ++i) {
        PrintProgress(i, inSamplingNumber);
        this->mSampledHist.push_back(this->SamplingEachHistogram());
    }

    std::vector<TH1D> bins;
    for (int i = 0; i < this->mNBins; ++i) {
        TH1D tempBin("","",1000, -10E+5, 10E+5);
        bins.push_back(tempBin);
    }

    std::cout << "making sample result" << std::endl;
    for (int i = 0; i < this->mNBins; ++i) {
        PrintProgress(i, this->mNBins);
        for (const auto& h : this->mSampledHist) {
            bins.at(i).Fill(h.GetBinContent(i+1));
        }
        this->mSampleResult->SetBinContent(i+1, bins.at(i).GetMean());
        this->mSampleResult->SetBinError(i+1, bins.at(i).GetRMS());
        //this->mSampleResult->SetBinContent(i+1, 
        //        this->GetHistCombinedNominal()->GetBinContent(i+1));
        //this->mSampleResult->SetBinError(i+1, 
        //        this->GetHistCombinedNominal()->GetBinError(i+1));
    }

    std::clock_t END = std::clock();
    std::cout << "took " << (END - START) / CLOCKS_PER_SEC << "s" << std::endl;
    std::cout << "this->mSampledHist.size(): ";
    std::cout << this->mSampledHist.size()  << std::endl;
}
//------------------------------------------------------------------------------
void FCN::SaveHist(std::string_view name) {
    std::cout << __func__ << std::endl;
    TCanvas can;
    can.Divide(2,2);
    can.cd(1);
    this->GetHistGenieNominal()->Draw();
    can.cd(2);
    this->GetHistG4Nominal()->Draw();
    can.cd(3);
    this->GetHistCombinedNominal()->Draw();
    can.cd(4);

    can.SaveAs(name + ".pdf");
}
//------------------------------------------------------------------------------
///////////////////////////////prepare matries//////////////////////////////////
void FCN::SetCorrelationMatrix() {
    std::cout << __func__ << std::endl;
    TMatrixD* tempCorMat = new TMatrixD(this->mNBins, this->mNBins);

    TMatrixD covMat(this->mNBins, this->mNBins);
    TMatrixD XY(this->mNBins, this->mNBins);
    TVectorD X(this->mNBins);

    for (const auto& h : this->mSampledHist) {
        for (int i = 0; i < this->mNBins; ++i) {
            X(i) += h.GetBinContent(i+1);
            for (int j = 0; j < this->mNBins; ++j) {
                XY(i, j) += h.GetBinContent(i+1) * h.GetBinContent(j+1);
            }
        }
    }
    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            double Eij = XY(i, j) / this->mSampledHist.size();
            double Ei = X(i) / this->mSampledHist.size();
            double Ej = X(j) / this->mSampledHist.size();
            double tempNominator = Eij - Ei * Ej;
            double tempDenominator = Ei * Ej;
            covMat(i, j) = tempNominator / tempDenominator;
            //covMat(i, j) = Eij - Ei - Ej;
        }
    }
    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            (*tempCorMat)(i, j) = std::pow(covMat(i, j), 2) 
                                  / (covMat(i, i) * covMat(j, j));
        }
    }

    this->mCorrelationMatrix = tempCorMat;

}
//------------------------------------------------------------------------------
void FCN::SetToyCorrelationMatrix() {
    TMatrixD* tempCorMat = new TMatrixD(this->mNBins, this->mNBins);

    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            if (i == j)
                (*tempCorMat)(i, j) = 1.;
            else if (i == 0 || j == 0)
                (*tempCorMat)(i, j) = -TOYCORR;
            else if (i == 2 || j == 2)
                (*tempCorMat)(i, j) = TOYCORR;
        }
    }

    this->mToyCorrelationMatrix = tempCorMat;
}
//------------------------------------------------------------------------------
void FCN::SetCovarianceMatrix() {
    std::cout << __func__ << std::endl;
    TMatrixD* tempCovMat = new TMatrixD(this->mNBins, this->mNBins);
    TMatrixD parError(this->mNBins, this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        if (this->mBinStep*(i+1) <= 300)
            parError(i,i) = 1.;
        else
            parError(i,i) = 0.1;
    }

    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            (*tempCovMat)(i, j) = parError(i, i) 
                * parError(j, j) * (*this->mCorrelationMatrix)(i,j); 
        }
    }

    this->mCovarianceMatrix = tempCovMat;
}
//------------------------------------------------------------------------------
void FCN::SetToyCovarianceMatrix() {
    TMatrixD* tempCovMat = new TMatrixD(this->mNBins, this->mNBins);
    TMatrixD parError(this->mNBins, this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        if (this->mBinStep*(i+1) <= 300)
            parError(i,i) = 1.;
        else
            parError(i,i) = 0.1;
    }

    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            (*tempCovMat)(i, j) = parError(i, i) 
                * parError(j, j) * (*this->mToyCorrelationMatrix)(i,j); 
        }
    }

    this->mToyCovarianceMatrix = tempCovMat;
}
//------------------------------------------------------------------------------
/////////////////////////////////////chi2///////////////////////////////////////
double FCN::PredictionAndData() const {
    TH1D n1NuPrediction("","1n high nu prediction", 
            this->mNBins, 0, this->mBinStep * this->mNBins);

    for (int i = 0; i < this->mNBins; ++i) {
        n1NuPrediction.SetBinContent(i+1, this->mN1Data->GetBinContent(i+1) 
                * ((RooRealVar*)this->_pulls->at(i))->getValV());
    }

    TVectorD difference(this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        if (this->mBinStep * (i + 1) <= 300)
            difference[i] = 0;
        else
            difference[i] = mN1Data->GetBinContent(i+1) 
                            - n1NuPrediction.GetBinContent(i+1);
        //std::cout << "difference[" << i << "]: ";
        //std::cout << difference[i] << std::endl;
    }
    TMatrixD stat(this->mNBins, this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        stat(i, i) = 0.1 * this->mN1Data->GetBinContent(i+1);
    }
    TVectorD mulVec(difference);
    stat.Invert();
    mulVec *= stat;
    return mulVec * difference;
}
//------------------------------------------------------------------------------
double FCN::PenaltyForParameters() const {
    TVectorD pars(this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        pars[i] = ((RooRealVar*)this->_pulls->at(i))->getValV() - 1;
    }
    TMatrixD invertCov(*(this->mCovarianceMatrix));
    TVectorD mulVec2(pars);
    invertCov.Invert();
    mulVec2 *= invertCov;
    return mulVec2 * pars;
}
//------------------------------------------------------------------------------
double FCN::PenaltyForParametersToyModel() const {
    TVectorD pars(this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        pars[i] = ((RooRealVar*)this->_pulls->at(i))->getValV() - 1;
    }
    TMatrixD invertToyCov(*(this->mToyCovarianceMatrix));
    TVectorD mulVec2(pars);
    invertToyCov.Invert();
    mulVec2 *= invertToyCov;
    return mulVec2 * pars;
}
//------------------------------------------------------------------------------
Double_t FCN::evaluate() const {
    double chi2;
    chi2 += this->PredictionAndData();

    if (!TOY) {
        chi2 += this->PenaltyForParameters();
    } else {
        chi2 += this->PenaltyForParametersToyModel();
    }

    //std::cout << "chi2: " << chi2 << std::endl;
    return chi2;
}
