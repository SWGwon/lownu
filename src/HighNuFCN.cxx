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
                    Form("par%d", i+1), 0);
            tempPar->setConstant(false);
            if (this->mBinStep*(i+1) <= 300)
                tempPar->setError(1);
            else
                tempPar->setError(1);
            mParVec.push_back(tempPar);
            mPulls->add(*tempPar);
        }
        this->addServerList(*mPulls);

        this->SetHistGenieNominalError();
        //this->SetHistG4NominalError();
        this->SetHistCombinedError();
        this->SamplingHistograms(1000);
        this->SetCorrelationMatrix();
        this->SetCovarianceMatrix();
        this->SetToyCorrelationMatrix();
        this->SetToyCovarianceMatrix();

        this->SetTestCov();
}
//------------------------------------------------------------------------------
std::unique_ptr<TH1D> HighNuFCN::FillHist(std::string inputFile) {
    std::unique_ptr<TH1D> tempHist = std::make_unique<TH1D>(inputFile.c_str(), 
            "nominal", this->mNBins, 0, this->mNBins * this->mBinStep);
            //inputFile.c_str(), this->mNBins, 0, this->mNBins * this->mBinStep);

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
        tempHist->Fill(recoNu);
    }
    //tempHist->Scale(1/tempHist->Integral(), "nosw2");

    return std::move(tempHist);
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistGenieNominal(std::string inputFile) {
    std::cout << __func__ << std::endl;
    mHistGenieNominal = FillHist(inputFile);
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistGenieShift(std::string inputFile) {
    std::cout << __func__ << std::endl;
    mHistGenieShift = std::make_unique<TH1D>("", 
            "reweighted", this->mNBins, 0, this->mNBins * this->mBinStep);

    TFile tempFile(inputFile.c_str());
    if (!tempFile.IsOpen()) {
        std::cout << "in " << __func__ << std::endl;
        std::cout << "invalid input: " << inputFile << std::endl;
    }
    TTree* tempTree = (TTree*)tempFile.Get("tree");
    float recoNu;          tempTree->SetBranchAddress("recoNu", &recoNu);
    float reWeight;          tempTree->SetBranchAddress("reWeight", &reWeight);
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
        mHistGenieShift->Fill(recoNu, reWeight);
    }
    mHistGenieShift->Scale(mHistGenieNominal->Integral()/mHistGenieShift->Integral(), "nosw2");
    //mHistGenieShift = FillHist(inputFile);
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistGenieNominalError() {

    //SetHistGenieNominal("analysis_output_G1801a00000AfterCut.root");
    //SetHistGenieShift("analysis_output_G1802a00000AfterCut.root");
    SetHistGenieNominal("reweightAfterCut.root");
    SetHistGenieShift("reweightAfterCut.root");

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
    std::cout << __func__ << std::endl;
    mHistG4Nominal = FillHist(inputFile);
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistG4Shift(std::string inputFile) {
    std::cout << __func__ << std::endl;
    mHistG4Shift = FillHist(inputFile);
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistG4NominalError() {

    SetHistG4Nominal("G4comparisonDefaultAfterCut.root");
    SetHistG4Shift("G4comparisoninclAfterCut.root");

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
        double tempGenieFracError = this->GetHistGenieNominal()->GetBinError(i+1) 
                                    / this->GetHistGenieNominal()->GetBinContent(i+1);
        //fractional   
        //double tempG4FracError = this->GetHistG4Nominal()->GetBinError(i+1) 
        //                  / this->GetHistG4Nominal()->GetBinContent(i+1);

        this->GetHistCombinedNominal()->SetBinContent(i+1, 
                this->GetHistGenieNominal()->GetBinContent(i+1));
        this->GetHistCombinedNominal()->SetBinError(i+1, 
                //std::pow(std::pow(tempGenieFracError, 2) + std::pow(tempG4FracError, 2), 0.5) 
                tempGenieFracError * this->GetHistCombinedNominal()->GetBinContent(i+1));
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
        tempSample.SetBinContent(i+1, //std::max((double)0, 
                random.Gaus(this->GetHistCombinedNominal()->GetBinContent(i+1),
                            this->GetHistCombinedNominal()->GetBinError(i+1))
                //)
                );
    }
    tempSample.Scale(
            this->GetHistCombinedNominal()->Integral() / tempSample.Integral(),
            "nosw2"
            );

    return tempSample;
}
//------------------------------------------------------------------------------
void HighNuFCN::SamplingHistograms(int inSamplingNumber) {
    std::cout << __func__ << std::endl;
    mSampleResult = std::make_unique<TH1D>("","sampling result", 
            this->mNBins, 0, this->mBinStep * this->mNBins);
    
    for (int i = 0; i < inSamplingNumber; ++i) {
        PrintProgress(i, inSamplingNumber);
        this->mSampledHist.push_back(this->SamplingEachHistogram());
    }

    std::vector<TH1D> bins;
    for (int i = 0; i < this->mNBins; ++i) {
        TH1D tempBin("","",20000, 0, 20000);
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
void HighNuFCN::SaveHist(std::string_view name) const {
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
    TMatrixD tempCorr(this->mNBins, this->mNBins);
    mCorrelationMatrix = std::make_unique<TMatrixD>(this->mNBins, this->mNBins);

    TMatrixD covMat(this->mNBins, this->mNBins);
    TMatrixD XX(this->mNBins, this->mNBins);
    TVectorD X(this->mNBins);

    for (const auto& h : this->mSampledHist) {
        for (int i = 0; i < this->mNBins; ++i) {
            X(i) += h.GetBinContent(i+1);
            for (int j = 0; j < this->mNBins; ++j) {
                XX(i, j) += h.GetBinContent(i+1) * h.GetBinContent(j+1);
            }
        }
    }
    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            double Eij = XX(i, j) / this->mSampledHist.size();
            double Ei = X(i) / this->mSampledHist.size();
            double Ej = X(j) / this->mSampledHist.size();
            covMat(i, j) = Eij - Ei * Ej;
        }
    }
    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            tempCorr(i, j) = covMat(i, j) / (fracError.at(i) * fracError.at(j));
        }
    }
    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            //normalize , use absolute sign(if not doesn't work)
            (*mCorrelationMatrix)(i, j) = std::abs(tempCorr(i, j) / std::pow(std::abs(tempCorr(i, i) * tempCorr(j, j)), 0.5));
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
            parError(i,i) = 1;
    }

    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            (*mCovarianceMatrix)(i, j) = (*this->mCorrelationMatrix)(i,j) * (parError(i, i) * parError(j, j)); 
        }
    }
}
//------------------------------------------------------------------------------
/////////////////////////////////////chi2///////////////////////////////////////
double HighNuFCN::PredictionMinusData() const {
    TH1D n1NuPrediction("","1n high nu prediction", 
            this->mNBins, 0, this->mBinStep * this->mNBins);

    for (int i = 0; i < this->mNBins; ++i) {
        n1NuPrediction.SetBinContent(i+1, this->mHistGenieNominal->GetBinContent(i+1) 
                * (1 + mParVec.at(i)->getValV()));
    }

    TVectorD difference(this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        if (this->mBinStep * (i + 1) <= 300)
            difference[i] = 0;
        else
            difference[i] = mHistGenieNominal->GetBinContent(i+1) 
                            - n1NuPrediction.GetBinContent(i+1);
        //std::cout << "difference[" << i << "]: ";
        //std::cout << difference[i] << std::endl;
    }
    TMatrixD stat(this->mNBins, this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            if (i == j)
                stat(i, j) = 1.0 * this->mHistGenieNominal->GetBinContent(i+1);
            else
                stat(i, j) = 0.;
        }
    }
    TVectorD mulVec(difference);
    stat.Invert();
    mulVec *= stat;

    return mulVec * difference;
}
//------------------------------------------------------------------------------
double HighNuFCN::Correlation() const {
    TVectorD pars(this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        pars[i] = mParVec.at(i)->getValV();
    }
    TMatrixD invertCov(*(this->mCovarianceMatrix));
    TVectorD mulVec2(pars);
    invertCov.Invert();
    mulVec2 *= invertCov;
    return mulVec2 * pars;
}
//------------------------------------------------------------------------------
double HighNuFCN::ExtraPenaltyForParameters() const {
    double a = 0;
    for (int i = 0; i < this->mNBins; ++i) {
        if (this->mBinStep*(i+1) <= 300)
            a += std::pow(mParVec.at(i)->getValV() - 1, 2)/1;
        else
            a += std::pow(mParVec.at(i)->getValV() - 1, 2)/0.1;
    }
    return a;
}
//------------------------------------------------------------------------------
Double_t HighNuFCN::evaluate() const {
    double chi2 = 0;
    chi2 += this->PredictionMinusData();

    if (!TOY) {
        chi2 += this->Correlation();
    } else {
        chi2 += this->CorrelationToyModel();
    }

    //for (int i = 0; i < this->mNBins; ++i) {
    //    std::cout << "p[" << i << "]: " << mParVec.at(i)->getValV() << std::endl;
    //}
    //std::cout << "P-D: " << this->PredictionMinusData() << std::endl;
    //std::cout << "P-CV: " << this->Correlation() << std::endl;
    //std::cout << "chi2: " << chi2 << std::endl;
    
    //chi2 += this->TestChi2();
    //chi2 += this->ExtraPenaltyForParameters();

    mParV.push_back(mParVec.at(1)->getValV());
    mChi2.push_back(chi2);

    return chi2;
}
//------------------------------------------------------------------------------
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
//for test
void HighNuFCN::SetTestCov() {
    TMatrixD tempCorr(this->mNBins, this->mNBins);
    mTestCov = std::make_unique<TMatrixD>(this->mNBins, this->mNBins);

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
            (*mTestCov)(i, j) = tempNominator;
        }
    }

    TCanvas c;
    mTestCov->Draw("colz text");
    c.SaveAs("testCov.pdf");

    TMatrixD invertCov(*(this->mTestCov));
    invertCov.Invert();
    TCanvas c1;
    invertCov.Draw("colz text");
    c1.SaveAs("inverttestCov.pdf");
}
//------------------------------------------------------------------------------
double HighNuFCN::TestChi2() const {
    TH1D n1NuPrediction("","1n high nu prediction", 
            this->mNBins, 0, this->mBinStep * this->mNBins);

    for (int i = 0; i < this->mNBins; ++i) {
        n1NuPrediction.SetBinContent(i+1, this->mHistGenieNominal->GetBinContent(i+1) 
                * mParVec.at(i)->getValV());
    }

    TVectorD difference(this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        if (this->mBinStep * (i + 1) <= 300)
            difference[i] = 0;
        else
            difference[i] = mHistGenieNominal->GetBinContent(i+1) 
                            - n1NuPrediction.GetBinContent(i+1);
        //std::cout << "difference[" << i << "]: ";
        //std::cout << difference[i] << std::endl;
    }
    TMatrixD invertCov(*(this->mTestCov));
    //for (int i = 0; i < this->mNBins; ++i) {
    //    invertCov(i, i) = 1.0 * this->mHistGenieNominal->GetBinContent(i+1);
    //}
    //for (int i = 0; i < this->mNBins; ++i) {
    //    for (int j = 0; j < this->mNBins; ++j) {
    //        if (i == j) continue;
    //        invertCov(i, j) = (*mCorrelationMatrix)(i, j) * std::pow(invertCov(i, i), 0.5) * std::pow(invertCov(j, j), 0.5);
    //    }
    //}
    TVectorD mulVec(difference);
    invertCov.Invert();
    mulVec *= invertCov;

    return mulVec * difference;
}
//------------------------------------------------------------------------------
void HighNuFCN::SetToyCorrelationMatrix() {
    mToyCorrelationMatrix = std::make_unique<TMatrixD>(this->mNBins, this->mNBins);

    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            if (i == j)
                (*mToyCorrelationMatrix)(i, j) = 1.;
            else if (i == 0 || j == 0)
                (*mToyCorrelationMatrix)(i, j) = TOYCORR;
            else if (i == 2 || j == 2)
                (*mToyCorrelationMatrix)(i, j) = TOYCORR;
        }
    }
}
//------------------------------------------------------------------------------
double HighNuFCN::CorrelationToyModel() const {
    TVectorD pars(this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        pars[i] = mParVec.at(i)->getValV();
    }
    TMatrixD invertToyCov(*(this->mToyCovarianceMatrix));
    TVectorD mulVec2(pars);
    invertToyCov.Invert();
    mulVec2 *= invertToyCov;
    return mulVec2 * pars;
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
            (*mToyCovarianceMatrix)(i, j) = (*this->mToyCorrelationMatrix)(i,j)
                * parError(i, i) * parError(j, j)  ; 
        }
    }
}
//------------------------------------------------------------------------------
