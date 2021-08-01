#include "Utils.hxx"
#include "HighNuFCN.hxx"

bool combinHighNu = false;
bool TOY = false;
double TOYCORR = 0;

std::vector<double> fracError;

/////////////////////////////prepare histogram//////////////////////////////////
HighNuFCN::HighNuFCN(const int inNBin, const int inBinStep)
    : mNBins(inNBin), mBinStep(inBinStep) {
        mPulls = std::make_unique<RooListProxy>("mPulls","mPulls",this);
        for (int i = 0; i < this->mNBins; i++) {
            RooRealVar* tempPar = new RooRealVar(Form("par%d", i), 
                    Form("par%d", i+1), 1);
            tempPar->setConstant(false);
            mParVec.push_back(tempPar);
            mPulls->add(*tempPar);
        }
        this->addServerList(*mPulls);

        this->SetHistGenieNominalError();
        this->SetHistG4NominalError();
        this->SetHistCombinedError();
        this->SamplingHistograms(1000);
        this->SetCorrelationMatrix();
        this->SetCovarianceMatrix();
        this->SetToyCorrelationMatrix();
        this->SetToyCovarianceMatrix();
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistGenieNominal(std::string inputFile) {
    mHistGenieNominal = std::make_unique<TH1D>("HistGenieNominal", 
            "HistGenieNominal;reco #nu;Normalized fraction", 
            this->mNBins, 0, this->mNBins * this->mBinStep);
    mN1Data = std::make_unique<TH1D>("HistN1Data", "HistN1Data", 
            this->mNBins, 0, this->mNBins * this->mBinStep);

    TFile tempFile(inputFile.c_str());
    if (!tempFile.IsOpen()) {
        std::cout << "in " << __func__ << std::endl;
        std::cout << "invalid input: " << inputFile << std::endl;
    }
    TTree* tempTree = (TTree*)tempFile.Get("tree");
    float recoNu;          tempTree->SetBranchAddress("recoNu", &recoNu);
    int numberOfFSNeutron; tempTree->SetBranchAddress("numberOfFSNeutron", 
            &numberOfFSNeutron);

    std::cout << __func__ << std::endl;
    for (int i = 0; i < tempTree->GetEntries(); ++i) {
        PrintProgress(i, tempTree->GetEntries());
        tempTree->GetEntry(i);
        //if (numberOfFSNeutron == 1) continue;
        if (combinHighNu) {
            if (300 < recoNu && recoNu < 550) recoNu = 325;
            if (550 < recoNu && recoNu < 800) recoNu = 375;
        }
        mHistGenieNominal->Fill(recoNu);
        //if (numberOfFSNeutron == 1)
        mN1Data->Fill(recoNu);
    }
    mN1Data->Scale(1/mN1Data->Integral(), "nosw2");
    mHistGenieNominal->Scale(1/mHistGenieNominal->Integral(), "nosw2");
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistGenieShift(std::string inputFile) {
    mHistGenieShift = std::make_unique<TH1D>("HistGenieShift", 
            "HistGenieShift;reco #nu;Normalized fraction", 
            this->mNBins, 0, this->mNBins * this->mBinStep);
    TFile tempFile(inputFile.c_str());
    if (!tempFile.IsOpen()) {
        std::cout << "in " << __func__ << std::endl;
        std::cout << "invalid input: " << inputFile << std::endl;
    }
    TTree* tempTree = (TTree*)tempFile.Get("tree");
    float recoNu;          tempTree->SetBranchAddress("recoNu", &recoNu);
    int numberOfFSNeutron; tempTree->SetBranchAddress("numberOfFSNeutron", 
            &numberOfFSNeutron);

    std::cout << __func__ << std::endl;
    for (int i = 0; i < tempTree->GetEntries(); ++i) {
        PrintProgress(i, tempTree->GetEntries());
        tempTree->GetEntry(i);
        //if (numberOfFSNeutron == 1) continue;
        if (combinHighNu) {
            if (300 < recoNu && recoNu < 550) recoNu = 325;
            if (550 < recoNu && recoNu < 800) recoNu = 375;
        }
        mHistGenieShift->Fill(recoNu);
    }
    mHistGenieShift->Scale(1/mHistGenieShift->Integral(), "nosw2");
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistGenieNominalError() {

    SetHistGenieNominal("analysis_output_G1801a00000AfterCut.root");
    SetHistGenieShift("analysis_output_G1802a00000AfterCut.root");

    this->GetHistGenieNominal()->Scale(
            1 / this->GetHistGenieNominal()->Integral(), "nosw2");
    this->GetHistGenieShift()->Scale(
            1 / this->GetHistGenieShift()->Integral(), "nosw2");

    for (int i = 0; i < this->mNBins; ++i) {
        double tempDifference = 0;
        tempDifference = this->GetHistGenieShift()->GetBinContent(i+1)
                             - this->GetHistGenieNominal()->GetBinContent(i+1);
        this->GetHistGenieNominal()->SetBinError(i+1, tempDifference);

        double tempFracError = tempDifference/(this->GetHistGenieNominal()->GetBinContent(i+1));
        fracError.push_back(tempFracError);
    }

}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistG4Nominal(std::string inputFile) {
    mHistG4Nominal = std::make_unique<TH1D>("HistG4Nominal", 
            "HistG4Nominal;reco #nu;Normalized fraction", 
            this->mNBins, 0, this->mNBins * this->mBinStep);

    TFile tempFile(inputFile.c_str());
    if (!tempFile.IsOpen()) {
        std::cout << "in " << __func__ << std::endl;
        std::cout << "invalid input: " << inputFile << std::endl;
    }
    TTree* tempTree = (TTree*)tempFile.Get("tree");
    float recoNu;          tempTree->SetBranchAddress("recoNu", &recoNu);
    int numberOfFSNeutron; tempTree->SetBranchAddress("numberOfFSNeutron", 
            &numberOfFSNeutron);

    std::cout << __func__ << std::endl;
    for (int i = 0; i < tempTree->GetEntries(); ++i) {
        PrintProgress(i, tempTree->GetEntries());
        tempTree->GetEntry(i);
        //if (numberOfFSNeutron == 1) continue;
        if (combinHighNu) {
            if (300 < recoNu && recoNu < 550) recoNu = 325;
            if (550 < recoNu && recoNu < 800) recoNu = 375;
        }
        mHistG4Nominal->Fill(recoNu);
    }
    mHistG4Nominal->Scale(1/mHistG4Nominal->Integral(), "nosw2");
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistG4Shift(std::string inputFile) {
    mHistG4Shift = std::make_unique<TH1D>("HistG4Shift", 
            "HistG4Shift;reco #nu;Normalized fraction", 
            this->mNBins, 0, this->mNBins * this->mBinStep);

    TFile tempFile(inputFile.c_str());
    if (!tempFile.IsOpen()) {
        std::cout << "in " << __func__ << std::endl;
        std::cout << "invalid input: " << inputFile << std::endl;
    }
    TTree* tempTree = (TTree*)tempFile.Get("tree");
    float recoNu;          tempTree->SetBranchAddress("recoNu", &recoNu);
    int numberOfFSNeutron; tempTree->SetBranchAddress("numberOfFSNeutron", 
            &numberOfFSNeutron);

    std::cout << __func__ << std::endl;
    for (int i = 0; i < tempTree->GetEntries(); ++i) {
        PrintProgress(i, tempTree->GetEntries());
        tempTree->GetEntry(i);
        //if (numberOfFSNeutron == 1) continue;
        if (combinHighNu) {
            if (300 < recoNu && recoNu < 550) recoNu = 325;
            if (550 < recoNu && recoNu < 800) recoNu = 375;
        }
        mHistG4Shift->Fill(recoNu);
    }
    mHistG4Shift->Scale(1/mHistG4Shift->Integral(), "nosw2");
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistG4NominalError() {

    SetHistG4Nominal("G4comparisonDefaultAfterCut.root");
    SetHistG4Shift("G4comparisoninclAfterCut.root");

    this->GetHistG4Nominal()->Scale(
            1 / this->GetHistG4Nominal()->Integral(), "nosw2");
    this->GetHistG4Shift()->Scale(
            1 / this->GetHistG4Shift()->Integral(), "nosw2");

    for (int i = 0; i < this->mNBins; ++i) {
        double tempDifference = 0;
        tempDifference = this->GetHistG4Shift()->GetBinContent(i+1)
                         - this->GetHistG4Nominal()->GetBinContent(i+1);

        this->GetHistG4Nominal()->SetBinError(i+1, tempDifference);
    }
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistCombinedError() {
    mHistCombinedNominal = std::make_unique<TH1D>("combinedNominal", 
            "combinedNominal;reco #nu;Normalized fraction", 
            this->mNBins, 0, this->mNBins * this->mBinStep);

    for (int i = 0; i < this->mNBins; ++i) {
        //fractional   
        double tempGenieError = this->GetHistGenieNominal()->GetBinError(i+1) 
                              / this->GetHistGenieNominal()->GetBinContent(i+1);
        //fractional   
        double tempG4Error = this->GetHistG4Nominal()->GetBinError(i+1) 
                           / this->GetHistG4Nominal()->GetBinContent(i+1);

        this->GetHistCombinedNominal()->SetBinContent(i+1, 
                this->GetHistGenieNominal()->GetBinContent(i+1));
        this->GetHistCombinedNominal()->SetBinError(i+1, 
                //std::pow(std::pow(tempGenieError, 2) + std::pow(tempG4Error, 2), 0.5) 
                tempGenieError
                * this->GetHistCombinedNominal()->GetBinContent(i+1));
    }
}
//------------------------------------------------------------------------------
////////////////////////////////sampling histogram//////////////////////////////
TH1D HighNuFCN::SamplingEachHistogram() {
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
void HighNuFCN::SamplingHistograms(int inSamplingNumber) {
    mSampleResult = std::make_unique<TH1D>("","sampling result", 
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

    std::cout << this->mSampledHist.size()  << std::endl;
}
//------------------------------------------------------------------------------
void HighNuFCN::SaveHist(std::string_view name) {
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
void HighNuFCN::SetCorrelationMatrix() {
    std::cout << __func__ << std::endl;
    mCorrelationMatrix = std::make_unique<TMatrixD>(this->mNBins, this->mNBins);

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
    TCanvas c;
    covMat.Draw("colz text");
    c.SaveAs("cov.pdf");
    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            (*mCorrelationMatrix)(i, j) = std::pow(covMat(i, j), 2) 
                                  / (covMat(i, i) * covMat(j, j));

        }
    }
}
//------------------------------------------------------------------------------
void HighNuFCN::SetToyCorrelationMatrix() {
    mToyCorrelationMatrix = std::make_unique<TMatrixD>(this->mNBins, this->mNBins);

    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            if (i == j)
                (*mToyCorrelationMatrix)(i, j) = 1.;
            else if (i == 0 || j == 0)
                (*mToyCorrelationMatrix)(i, j) = -TOYCORR;
            else if (i == 2 || j == 2)
                (*mToyCorrelationMatrix)(i, j) = TOYCORR;
        }
    }
}
//------------------------------------------------------------------------------
void HighNuFCN::SetCovarianceMatrix() {
    std::cout << __func__ << std::endl;
    mCovarianceMatrix = std::make_unique<TMatrixD>(this->mNBins, this->mNBins);
    TMatrixD parError(this->mNBins, this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        if (this->mBinStep*(i+1) <= 300)
            parError(i,i) = 1.;
        else
            parError(i,i) = 0.1;
    }

    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            (*mCovarianceMatrix)(i, j) = parError(i, i) 
                * parError(j, j) * (*this->mCorrelationMatrix)(i,j); 
        }
    }
}
//------------------------------------------------------------------------------
void HighNuFCN::SetToyCovarianceMatrix() {
    mToyCovarianceMatrix = std::make_unique<TMatrixD>(this->mNBins, this->mNBins);
    TMatrixD parError(this->mNBins, this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        if (this->mBinStep*(i+1) <= 300)
            parError(i,i) = 1.;
        else
            parError(i,i) = 0.1;
    }

    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            (*mToyCovarianceMatrix)(i, j) = parError(i, i) 
                * parError(j, j) * (*this->mToyCorrelationMatrix)(i,j); 
        }
    }
}
//------------------------------------------------------------------------------
/////////////////////////////////////chi2///////////////////////////////////////
double HighNuFCN::PredictionAndData() const {
    TH1D n1NuPrediction("","1n high nu prediction", 
            this->mNBins, 0, this->mBinStep * this->mNBins);

    for (int i = 0; i < this->mNBins; ++i) {
        n1NuPrediction.SetBinContent(i+1, this->mN1Data->GetBinContent(i+1) 
                * mParVec.at(i)->getValV());
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
double HighNuFCN::PenaltyForParameters() const {
    TVectorD pars(this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        pars[i] = ((RooRealVar*)this->mPulls->at(i))->getValV() - 1;
    }
    TMatrixD invertCov(*(this->mCovarianceMatrix));
    TVectorD mulVec2(pars);
    invertCov.Invert();
    mulVec2 *= invertCov;
    return mulVec2 * pars;
}
//------------------------------------------------------------------------------
double HighNuFCN::PenaltyForParametersToyModel() const {
    TVectorD pars(this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        pars[i] = ((RooRealVar*)this->mPulls->at(i))->getValV() - 1;
    }
    TMatrixD invertToyCov(*(this->mToyCovarianceMatrix));
    TVectorD mulVec2(pars);
    invertToyCov.Invert();
    mulVec2 *= invertToyCov;
    return mulVec2 * pars;
}
//------------------------------------------------------------------------------
Double_t HighNuFCN::evaluate() const {
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
