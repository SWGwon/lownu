#include "Utils.hxx"
#include "HighNuFCN.hxx"

bool combinHighNu = false;
bool TOY = false;
double TOYCORR = 0;


/////////////////////////////prepare histogram//////////////////////////////////
HighNuFCN::HighNuFCN(const int inNBin, const int inBinStep, bool N1HighNu, bool N1NonNeutron, bool N2, bool N3)
    : mNBins(inNBin), mBinStep(inBinStep), N1HighNu(N1HighNu), N1NonNeutron(N1NonNeutron), N2(N2), N3(N3) {
        std::cout << "N1HighNu: " << N1HighNu << std::endl;
        std::cout << "N1NonNeutron: " << N1NonNeutron << std::endl;
        std::cout << "N2: " << N2 << std::endl;
        std::cout << "N3: " << N3 << std::endl;
        SetEvents("output_10_reweightAfterCut.root");
        InitializeHistograms();

        if (!N1HighNu && !N1NonNeutron && !N2 && !N3) {
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
            this->mCorrelationMatrixAll = this->SetCorrelationMatrix(*this->mHistNominalAll, *this->mHistShiftAll);
            this->mCovarianceMatrixAll = this->SetCovarianceMatrix(*this->mCorrelationMatrixAll);
        }

        if (N1HighNu) {
            mPullsN1HighNu = std::make_unique<RooListProxy>("mPullsN1HighNu","mPullsN1HighNu",this);
            for (int i = 0; i < this->mNBins; i++) {
                RooRealVar* tempPar = new RooRealVar(Form("parN1HighNu%d", i), 
                        Form("parN1HighNu%d", i+1), 0);
                tempPar->setConstant(false);
                if (this->mBinStep*(i+1) <= 300) {
                    tempPar->setError(1);
                } else {
                    tempPar->setError(1);
                }
                mParVecN1HighNu.push_back(tempPar);
                mPullsN1HighNu->add(*tempPar);
            }
            this->addServerList(*mPullsN1HighNu);
            this->mCorrelationMatrixN1HighNu = this->SetCorrelationMatrix(*this->mHistNominalN1HighNu, *this->mHistShiftN1HighNu);
            this->mCovarianceMatrixN1HighNu = this->SetCovarianceMatrix(*this->mCorrelationMatrixN1HighNu);
        }

        if (N1NonNeutron) {
            mPullsN1NonNeutron = std::make_unique<RooListProxy>("mPullsN1NonNeutron","mPullsN1NonNeutron",this);
            for (int i = 0; i < this->mNBins; i++) {
                RooRealVar* tempPar = new RooRealVar(Form("parN1NonNeutron%d", i), 
                        Form("parN1NonNeutron%d", i+1), 0);
                tempPar->setConstant(false);
                if (this->mBinStep*(i+1) <= 300)
                    tempPar->setError(1);
                else
                    tempPar->setError(1);
                mParVecN1NonNeutron.push_back(tempPar);
                mPullsN1NonNeutron->add(*tempPar);
            }
            this->addServerList(*mPullsN1NonNeutron);
            this->mCorrelationMatrixN1NonNeutron = this->SetCorrelationMatrix(*this->mHistNominalN1NonNeutron, *this->mHistShiftN1NonNeutron);
            this->mCovarianceMatrixN1NonNeutron = this->SetCovarianceMatrix(*this->mCorrelationMatrixN1NonNeutron);
        }

        if (N2) {
            mPullsN2 = std::make_unique<RooListProxy>("mPullsN2","mPullsN2",this);
            for (int i = 0; i < this->mNBins; i++) {
                RooRealVar* tempPar = new RooRealVar(Form("parN2%d", i), 
                        Form("parN2%d", i+1), 0);
                tempPar->setConstant(false);
                if (this->mBinStep*(i+1) <= 300)
                    tempPar->setError(1);
                else
                    tempPar->setError(1);
                mParVecN2.push_back(tempPar);
                mPullsN2->add(*tempPar);
            }
            this->addServerList(*mPullsN2);
            this->mCorrelationMatrixN2 = this->SetCorrelationMatrix(*this->mHistNominalN2, *this->mHistShiftN2);
            this->mCovarianceMatrixN2 = this->SetCovarianceMatrix(*this->mCorrelationMatrixN2);
        }

        if (N3) {
            mPullsN3 = std::make_unique<RooListProxy>("mPullsN3","mPullsN3",this);
            for (int i = 0; i < this->mNBins; i++) {
                RooRealVar* tempPar = new RooRealVar(Form("parN3%d", i), 
                        Form("parN3%d", i+1), 0);
                tempPar->setConstant(false);
                if (this->mBinStep*(i+1) <= 300)
                    tempPar->setError(1);
                else
                    tempPar->setError(1);
                mParVecN3.push_back(tempPar);
                mPullsN3->add(*tempPar);
            }
            //this->addServerList(*mPullsN3);
            this->mCorrelationMatrixN3 = this->SetCorrelationMatrix(*this->mHistNominalN3, *this->mHistShiftN3);
            this->mCovarianceMatrixN3 = this->SetCovarianceMatrix(*this->mCorrelationMatrixN3);
        }

        //this->SetToyCorrelationMatrix();
        //this->SetToyCovarianceMatrix();

        //this->SetTestCov();
}
//------------------------------------------------------------------------------
void HighNuFCN::SetEvents(std::string inputFile) {
    TFile tempFile(inputFile.c_str());
    if (!tempFile.IsOpen()) {
        std::cout << "in " << __func__ << std::endl;
        std::cout << "invalid input: " << inputFile << std::endl;
    }
    TTree* tempTree = (TTree*)tempFile.Get("tree");
    float recoNu;          tempTree->SetBranchAddress("recoNu", &recoNu);
    float trueNu;          tempTree->SetBranchAddress("trueNu", &trueNu);
    float reWeight[100];   tempTree->SetBranchAddress("reWeight", &reWeight);
    int numberOfFSNeutron; tempTree->SetBranchAddress("numberOfFSNeutron", 
            &numberOfFSNeutron);
    for (int i = 0; i < tempTree->GetEntries(); ++i) {
        PrintProgress(i, tempTree->GetEntries());
        tempTree->GetEntry(i);

        if (numberOfFSNeutron == 1 && trueNu < 300) continue;

        bool isOk = true;
        Event tempEvent;
        tempEvent.recoNu = recoNu;
        tempEvent.trueNu = trueNu;
        for (int j = 0; j < 10; ++j) {
            if (reWeight[j] < 0) {
                isOk = false;
                break;
            }
            tempEvent.reWeight[j] = reWeight[j];
        }
        tempEvent.numberOfFSNeutron = numberOfFSNeutron;
        if (isOk)
            this->mEvents.push_back(tempEvent);
    }
}
//------------------------------------------------------------------------------
void HighNuFCN::InitializeHistograms() {
    mHistNominalN1NonNeutron = std::make_unique<TH1D>("", 
            "n1 non neutron nominal;reco #nu (MeV)", this->mNBins, 0, this->mNBins * this->mBinStep);
    mHistNominalN1HighNu = std::make_unique<TH1D>("", 
            "n1 high nu nominal;reco #nu (MeV)", this->mNBins, 0, this->mNBins * this->mBinStep);
    mHistNominalN2 = std::make_unique<TH1D>("", 
            "n2 nominal;reco #nu (MeV)", this->mNBins, 0, this->mNBins * this->mBinStep);
    mHistNominalN3 = std::make_unique<TH1D>("", 
            "n3 nominal;reco #nu (MeV)", this->mNBins, 0, this->mNBins * this->mBinStep);

    mHistShiftN1NonNeutron = std::make_unique<TH1D>("", 
            "n1 non neutron reweighted;reco #nu (MeV)", this->mNBins, 0, this->mNBins * this->mBinStep);
    mHistShiftN1HighNu = std::make_unique<TH1D>("", 
            "n1 high nu reweighted;reco #nu (MeV)", this->mNBins, 0, this->mNBins * this->mBinStep);
    mHistShiftN2 = std::make_unique<TH1D>("", 
            "n2 reweighted;reco #nu (MeV)", this->mNBins, 0, this->mNBins * this->mBinStep);
    mHistShiftN3 = std::make_unique<TH1D>("", 
            "n3 reweighted;reco #nu (MeV)", this->mNBins, 0, this->mNBins * this->mBinStep);

    this->SetHistograms();
}
//------------------------------------------------------------------------------
std::unique_ptr<TH1D> HighNuFCN::FillHist(std::string inputFile) {
    std::unique_ptr<TH1D> tempHist = std::make_unique<TH1D>(inputFile.c_str(), 
            "nominal;reco #nu (MeV)", this->mNBins, 0, this->mNBins * this->mBinStep);
            //inputFile.c_str(), this->mNBins, 0, this->mNBins * this->mBinStep);

    TFile tempFile(inputFile.c_str());
    if (!tempFile.IsOpen()) {
        std::cout << "in " << __func__ << std::endl;
        std::cout << "invalid input: " << inputFile << std::endl;
    }
    TTree* tempTree = (TTree*)tempFile.Get("tree");
    float recoNu;          tempTree->SetBranchAddress("recoNu", &recoNu);
    float trueNu;          tempTree->SetBranchAddress("trueNu", &trueNu);
    int numberOfFSNeutron; tempTree->SetBranchAddress("numberOfFSNeutron", 
            &numberOfFSNeutron);

    std::cout << __func__ << std::endl;
    for (int i = 0; i < tempTree->GetEntries(); ++i) {
        PrintProgress(i, tempTree->GetEntries());
        tempTree->GetEntry(i);
        tempHist->Fill(recoNu);
    }
    //for (auto event : this->mEvents) {
    //    PrintProgress(i, tempTree->GetEntries());
    //    tempHist->Fill(recoNu);
    //}
    //tempHist->Scale(1/tempHist->Integral(), "nosw2");

    return std::move(tempHist);
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistGenieNominal(std::string inputFile) {
    std::cout << __func__ << std::endl;
    mHistNominalAll = std::make_unique<TH1D>(inputFile.c_str(), 
            "nominal;reco #nu (MeV)", this->mNBins, 0, this->mNBins * this->mBinStep);

    int i = 0;
    for (auto e : this->mEvents) {
        PrintProgress(i, this->mEvents.size());
        mHistNominalAll->Fill(e.recoNu);
        if (e.numberOfFSNeutron == 1 && e.trueNu < 300)
            mHistNominalN1NonNeutron->Fill(e.recoNu);
        if (e.numberOfFSNeutron == 1 && e.trueNu > 300)
            mHistNominalN1HighNu->Fill(e.recoNu);
        if (e.numberOfFSNeutron == 2)
            mHistNominalN2->Fill(e.recoNu);
        if (e.numberOfFSNeutron > 2)
            mHistNominalN3->Fill(e.recoNu);
        ++i;
    }
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistGenieShift(std::string inputFile) {
    std::cout << __func__ << std::endl;
    mHistShiftAll = std::make_unique<TH1D>("", 
            "reweighted;reco #nu (MeV)", this->mNBins, 0, this->mNBins * this->mBinStep);

    int i = 0;
    for (auto e : this->mEvents) {
        PrintProgress(i, this->mEvents.size());
        double tempAvgReWeight = 0; 
        for (int j = 0; j < 10; ++j) {
            tempAvgReWeight += e.reWeight[j];
        }
        tempAvgReWeight /= 10;
        //std::cout << "tempAvgReWeight: " << tempAvgReWeight << std::endl;
        
        mHistShiftAll->Fill(e.recoNu, tempAvgReWeight);
        if (e.numberOfFSNeutron == 1 && e.trueNu < 300) 
            mHistShiftN1NonNeutron->Fill(e.recoNu, tempAvgReWeight);
        if (e.numberOfFSNeutron == 1 && e.trueNu > 300) 
            mHistShiftN1HighNu->Fill(e.recoNu, tempAvgReWeight);
        if (e.numberOfFSNeutron == 2)
            mHistShiftN2->Fill(e.recoNu, tempAvgReWeight);
        if (e.numberOfFSNeutron > 2)
            mHistShiftN3->Fill(e.recoNu, tempAvgReWeight);
        ++i;
    }
}
//------------------------------------------------------------------------------
void HighNuFCN::SetHistograms() {
    //SetHistGenieNominal("analysis_output_G1801a00000AfterCut.root");
    //SetHistGenieShift("analysis_output_G1802a00000AfterCut.root");

    SetHistGenieNominal("output_10loop_reweightAfterCut.root");
    SetHistGenieShift("output_10loop_reweightAfterCut.root");
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
TH1D HighNuFCN::SamplingEachHistogram(const TH1D& inNominal, const TH1D& inShift) {
    TH1D tempSample("","tempSample", this->mNBins, 0, 
                    this->mBinStep * this->mNBins);

    TH1D tempNominal(inNominal);
    for (int i = 0; i < this->mNBins; ++i) {
        tempNominal.SetBinError(i+1, inNominal.GetBinContent(i+1) - inShift.GetBinContent(i+1));
    }

    TRandom random;
    random.SetSeed(0);
    for (int i = 0; i < this->mNBins; ++i) {
        tempSample.SetBinContent(i+1, //std::max((double)0, 
                random.Gaus(inShift.GetBinContent(i+1),
                            tempNominal.GetBinError(i+1))
                //)
                );
    }
    //tempSample.Scale(
    //        tempNominal.Integral() / tempSample.Integral(),
    //        "nosw2"
    //        );

    return tempSample;
}
//------------------------------------------------------------------------------
std::vector<TH1D> HighNuFCN::SamplingHistograms(int inSamplingNumber, const TH1D& inNominal, const TH1D& inShift) {
    std::cout << __func__ << std::endl;
    /*
    std::vector<TH1D> tempResult(inSamplingNumber);

    mSampleResult = std::make_unique<TH1D>("","sampling result", 
            this->mNBins, 0, this->mBinStep * this->mNBins);

    for (int i = 0; i < inSamplingNumber; ++i) {
        PrintProgress(i, inSamplingNumber);
        tempResult.at(i) = this->SamplingEachHistogram(inNominal, inShift);
        //this->mSampledHist.push_back(this->SamplingEachHistogram(inNominal, inShift));
    }
    */

    std::vector<TH1D> tempResult;

    //TFile tempFile("output_10_reweightAfterCut.root");

    //TTree* tempTree = (TTree*)tempFile.Get("tree");
    //float recoNu;          tempTree->SetBranchAddress("recoNu", &recoNu);
    //float reWeight[100];        tempTree->SetBranchAddress("reWeight", &reWeight);
    //float trueNu;          tempTree->SetBranchAddress("trueNu", &trueNu);
    //int numberOfFSNeutron; tempTree->SetBranchAddress("numberOfFSNeutron", 
    //        &numberOfFSNeutron);

    //for (int j = 0; j < 10; ++j) {
    //    std::cout << "loop: " << j << std::endl;
    //    TH1D tempHist("", "sample", this->mNBins, 0, this->mNBins * this->mBinStep);
    //    for (int i = 0; i < tempTree->GetEntries(); ++i) {
    //        PrintProgress(i, tempTree->GetEntries());
    //        tempTree->GetEntry(i);
    //        if (reWeight[j] == 0 || reWeight[j] < 0)
    //            continue;

    //        if (!N1HighNu && !N1NonNeutron && !N2 && !N3) {
    //            tempHist.Fill(recoNu, reWeight[j]);
    //        }
    //        if (N1HighNu) {
    //            if (numberOfFSNeutron == 1 && trueNu > 300) 
    //                tempHist.Fill(recoNu, reWeight[j]);
    //        }
    //        if (N1NonNeutron) {
    //            if (numberOfFSNeutron == 1 && trueNu < 300) 
    //                tempHist.Fill(recoNu, reWeight[j]);
    //        }
    //        if (N2) {
    //            if (numberOfFSNeutron == 2)
    //                tempHist.Fill(recoNu, reWeight[j]);
    //        }
    //        if (N3) {
    //            if (numberOfFSNeutron > 2)
    //                tempHist.Fill(recoNu, reWeight[j]);
    //        }
    //    }
    //    //tempHist.Scale(this->mHistNominalAll->Integral()/tempHist.Integral(), "nosw2");
    //    tempResult.push_back(tempHist);
    //}
    
    for (int j = 0; j < 10; ++j) {
        std::cout << "loop: " << j << std::endl;
        TH1D tempHist("", "sample", this->mNBins, 0, this->mNBins * this->mBinStep);
        int i = 0;
        for (auto e : this->mEvents) {
            PrintProgress(i, this->mEvents.size());
            ++i;
            if (e.reWeight[j] == 0 || e.reWeight[j] < 0)
                continue;

            if (!N1HighNu && !N1NonNeutron && !N2 && !N3) {
                tempHist.Fill(e.recoNu, e.reWeight[j]);
            }
            if (N1HighNu) {
                if (e.numberOfFSNeutron == 1 && e.trueNu > 300) 
                    tempHist.Fill(e.recoNu, e.reWeight[j]);
            }
            if (N1NonNeutron) {
                if (e.numberOfFSNeutron == 1 && e.trueNu < 300) 
                    tempHist.Fill(e.recoNu, e.reWeight[j]);
            }
            if (N2) {
                if (e.numberOfFSNeutron == 2)
                    tempHist.Fill(e.recoNu, e.reWeight[j]);
            }
            if (N3) {
                if (e.numberOfFSNeutron > 2)
                    tempHist.Fill(e.recoNu, e.reWeight[j]);
            }
        }
        //tempHist.Scale(this->mHistNominalAll->Integral()/tempHist.Integral(), "nosw2");
        tempResult.push_back(tempHist);
    }

    return tempResult;
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
std::unique_ptr<TMatrixD> HighNuFCN::SetCorrelationMatrix(const TH1D& inNominal, const TH1D& inShift) {
    std::cout << __func__ << std::endl;
    TMatrixD tempCorr(this->mNBins, this->mNBins);

    auto tempCorrelationMatrix = std::make_unique<TMatrixD>(this->mNBins, this->mNBins);

    TMatrixD covMat(this->mNBins, this->mNBins);
    TMatrixD XX(this->mNBins, this->mNBins);
    TVectorD X(this->mNBins);

    std::vector<TH1D> tempSampledHist = this->SamplingHistograms(1000, inNominal, inShift);

    for (const auto& h : tempSampledHist) {
        for (int i = 0; i < this->mNBins; ++i) {
            X(i) += h.GetBinContent(i+1);
            for (int j = 0; j < this->mNBins; ++j) {
                XX(i, j) += h.GetBinContent(i+1) * h.GetBinContent(j+1);
            }
        }
    }
    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            double Eij = XX(i, j) / tempSampledHist.size();
            double Ei = X(i) / tempSampledHist.size();
            double Ej = X(j) / tempSampledHist.size();
            covMat(i, j) = Eij - Ei * Ej;
        }
    }

    std::vector<double> fracError;
    for (int i = 0; i < this->mNBins; ++i) {
        double tempDifference = 0;
        tempDifference = inNominal.GetBinContent(i+1)
                             - inShift.GetBinContent(i+1);

        double tempFracError = tempDifference;///(inNominal.GetBinContent(i+1));
        fracError.push_back(tempFracError);
    }

    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            tempCorr(i, j) = covMat(i, j) / (fracError.at(i) * fracError.at(j));
        }
    }
    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            //normalize , use absolute sign(if not doesn't work)
            //(*tempCorrelationMatrix)(i, j) = std::abs
            //    (
            //     tempCorr(i, j)
            //     / std::pow(std::abs(tempCorr(i, i) * tempCorr(j, j)), 0.5));
            (*tempCorrelationMatrix)(i, j) = //std::abs
                (
                 covMat(i, j) / std::pow(covMat(i, i) * covMat(j, j), 0.5)
                 );
        }
    }

    return std::move(tempCorrelationMatrix);
}
//------------------------------------------------------------------------------
std::unique_ptr<TMatrixD> HighNuFCN::SetCovarianceMatrix(const TMatrixD& inCorrMat) {
    std::cout << __func__ << std::endl;
    auto tempCovarianceMatrix = std::make_unique<TMatrixD>(this->mNBins, this->mNBins);
    TMatrixD parError(this->mNBins, this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        if (this->mBinStep*(i+1) <= 300)
            parError(i,i) = 1.;
        else
            parError(i,i) = 1.;
    }

    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            (*tempCovarianceMatrix)(i, j) = inCorrMat(i,j) * (parError(i, i) * parError(j, j)); 
        }
    }

    return std::move(tempCovarianceMatrix);
}
//------------------------------------------------------------------------------
/////////////////////////////////////chi2///////////////////////////////////////
double HighNuFCN::DoMatrixCal(const TVectorD& inVec, TMatrixD inMat) const {
    double tempResult = 0;

    TVectorD vec1(inVec);
    TVectorD vec2(inVec);
    TMatrixD invertMat = inMat.Invert();
    vec1 *= invertMat;
    tempResult += vec2 * vec1;

    return tempResult;
}
double HighNuFCN::PredictionMinusData(const TH1D& inNominal) const {
    TH1D prediction("","prediction", 
            this->mNBins, 0, this->mBinStep * this->mNBins);
    TH1D data("","data", 
            this->mNBins, 0, this->mBinStep * this->mNBins);

    for (int i = 0; i < this->mNBins; ++i) {
        double tempDataBin = 0;
        double tempPredBin = 0;
        if (!N1HighNu && !N1NonNeutron && !N2 && !N3) {
            tempDataBin += this->mHistNominalAll->GetBinContent(i+1);
            tempPredBin += this->mHistNominalAll->GetBinContent(i+1) * (1 + mParVec.at(i)->getValV());
        }
        if (N1HighNu) {
            tempDataBin += this->mHistNominalN1HighNu->GetBinContent(i+1);
            tempPredBin += this->mHistNominalN1HighNu->GetBinContent(i+1) * (1 + mParVecN1HighNu.at(i)->getValV());
        }
        if (N1NonNeutron) {
            tempDataBin += this->mHistNominalN1NonNeutron->GetBinContent(i+1);
            tempPredBin += this->mHistNominalN1NonNeutron->GetBinContent(i+1) * (1 + mParVecN1NonNeutron.at(i)->getValV());
        }
        if (N2) {
            tempDataBin += this->mHistNominalN2->GetBinContent(i+1);
            tempPredBin += this->mHistNominalN2->GetBinContent(i+1) * (1 + mParVecN2.at(i)->getValV());
        }
        if (N3) {
            tempDataBin += this->mHistNominalN3->GetBinContent(i+1);
            tempPredBin += this->mHistNominalN3->GetBinContent(i+1) * (1 + mParVecN3.at(i)->getValV());
        }
        data.SetBinContent(i+1, tempDataBin);
        prediction.SetBinContent(i+1, tempPredBin);
    }

    TVectorD difference(this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        if (this->mBinStep * (i + 1) <= 300)
            difference[i] = 0;
        else
            difference[i] = data.GetBinContent(i+1) 
                            - prediction.GetBinContent(i+1);
        //std::cout << "difference[" << i << "]: ";
        //std::cout << difference[i] << std::endl;
    }
    TMatrixD stat(this->mNBins, this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        for (int j = 0; j < this->mNBins; ++j) {
            if (i == j)
                stat(i, j) = 1.0 * data.GetBinContent(i+1);
            else
                stat(i, j) = 0.;
        }
    }

    return DoMatrixCal(difference, stat);
}
//------------------------------------------------------------------------------
double HighNuFCN::Correlation() const {
    double tempResult = 0;

    if (!N1HighNu && !N1NonNeutron && !N2 && !N3) {
        TVectorD pars(this->mNBins);
        for (int i = 0; i < this->mNBins; ++i) {
            pars[i] = mParVec.at(i)->getValV() - 0;
        }
        tempResult += DoMatrixCal(pars, *this->mCovarianceMatrixAll);
    }
    
    if (N1HighNu) {
        TVectorD parsN1HighNu(this->mNBins);
        for (int i = 0; i < this->mNBins; ++i) {
            parsN1HighNu[i] = mParVecN1HighNu.at(i)->getValV();
        }
        tempResult += DoMatrixCal(parsN1HighNu, *this->mCovarianceMatrixN1HighNu);
    }
    
    if (N1NonNeutron) {
        TVectorD parsN1NonNeutron(this->mNBins);
        for (int i = 0; i < this->mNBins; ++i) {
            parsN1NonNeutron[i] = mParVecN1NonNeutron.at(i)->getValV();
        }
        tempResult += DoMatrixCal(parsN1NonNeutron, *this->mCovarianceMatrixN1NonNeutron);
    }
    
    if (N2) {
        TVectorD parsN2(this->mNBins);
        for (int i = 0; i < this->mNBins; ++i) {
            parsN2[i] = mParVecN2.at(i)->getValV();
        }
        tempResult += DoMatrixCal(parsN2, *this->mCovarianceMatrixN2);
    }
    
    if (N3) {
        TVectorD parsN3(this->mNBins);
        for (int i = 0; i < this->mNBins; ++i) {
            parsN3[i] = mParVecN3.at(i)->getValV();
        }
        tempResult += DoMatrixCal(parsN3, *this->mCovarianceMatrixN3);
    }

    return tempResult;
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
    double P_D = this->PredictionMinusData(*this->mHistNominalN1HighNu);
    double CORR = this->Correlation();
    chi2 += P_D;
    chi2 += CORR;
    //if (!TOY) {
    //    chi2 += this->Correlation();
    //} else {
    //    chi2 += this->CorrelationToyModel();
    //}

    //for (int i = 0; i < this->mNBins; ++i) {
    //    std::cout << "p[" << i << "]: " << mParVec.at(i)->getValV() << std::endl;
    //}
    std::cout << "------------------------------" << std::endl;
    std::cout << "P-D: " << P_D << std::endl;
    std::cout << "P-CV: " << CORR << std::endl;
    std::cout << "chi2: " << chi2 << std::endl;
    
    //chi2 += this->TestChi2();
    //chi2 += this->ExtraPenaltyForParameters();

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
void HighNuFCN::SetParVec(std::vector<double> inVec) {
    if (inVec.size() != this->mParVecN1HighNu.size()) {
        std::cout << "in " << __func__ << ": " << std::endl;
        std::cout << "inVec.size() != this->mParVecN1HighNu.size()" << std::endl;
        return;
    }
    for (int i = 0; i < inVec.size(); ++i) {
        this->mParVecN1HighNu.at(i)->setVal(inVec.at(i));
    }
}
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

    TMatrixD invertCov = this->mTestCov->Invert();
    //invertCov.Invert();
    TCanvas c1;
    invertCov.Draw("colz text");
    c1.SaveAs("inverttestCov.pdf");
}
//------------------------------------------------------------------------------
double HighNuFCN::TestChi2() const {
    TH1D prediction("","1n high nu prediction", 
            this->mNBins, 0, this->mBinStep * this->mNBins);

    for (int i = 0; i < this->mNBins; ++i) {
        prediction.SetBinContent(i+1, this->mHistNominalAll->GetBinContent(i+1) 
                * mParVec.at(i)->getValV());
    }

    TVectorD difference(this->mNBins);
    for (int i = 0; i < this->mNBins; ++i) {
        if (this->mBinStep * (i + 1) <= 300)
            difference[i] = 0;
        else
            difference[i] = mHistNominalAll->GetBinContent(i+1) 
                            - prediction.GetBinContent(i+1);
        //std::cout << "difference[" << i << "]: ";
        //std::cout << difference[i] << std::endl;
    }
    TMatrixD invertCov = this->mTestCov->Invert();
    //for (int i = 0; i < this->mNBins; ++i) {
    //    invertCov(i, i) = 1.0 * this->mHistNominalAll->GetBinContent(i+1);
    //}
    //for (int i = 0; i < this->mNBins; ++i) {
    //    for (int j = 0; j < this->mNBins; ++j) {
    //        if (i == j) continue;
    //        invertCov(i, j) = (*mCorrelationMatrixAll)(i, j) * std::pow(invertCov(i, i), 0.5) * std::pow(invertCov(j, j), 0.5);
    //    }
    //}
    TVectorD mulVec(difference);
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
    TMatrixD invertToyCov = this->mToyCovarianceMatrix->Invert();
    TVectorD mulVec2(pars);
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
TMatrixD HighNuFCN::MyInvert(const TMatrixD& inMatrix) const {
    TMatrixD tempMatrix(inMatrix);
    int _1 = 0;
    int _2 = 1;
    int _3 = 2;
    double d1 = inMatrix(_1,_1) * (inMatrix(_2,_2) * inMatrix(_3,_3) - inMatrix(_2,_3) * inMatrix(_3,_2));
    double d2 = - inMatrix(_1,_2) * (inMatrix(_2,_1) * inMatrix(_3,_3) - inMatrix(_2,_3) * inMatrix(_3,_1));
    double d3 = + inMatrix(_1,_3) * (inMatrix(_2,_1) * inMatrix(_3,_2) - inMatrix(_2,_2) * inMatrix(_3,_1));
    double det = 0;
    //double det = inMatrix(_1,_1) * (inMatrix(_2,_2) * inMatrix(_3,_3) - inMatrix(_2,_3) * inMatrix(_3,_2))
    //           - inMatrix(_1,_2) * (inMatrix(_2,_1) * inMatrix(_3,_3) - inMatrix(_2,_3) * inMatrix(_3,_1))
    //           + inMatrix(_1,_3) * (inMatrix(_2,_1) * inMatrix(_3,_2) - inMatrix(_2,_2) * inMatrix(_3,_1));

    det = d1 + d2 + d3;

    std::cout << "d1: " << d1 << std::endl;
    std::cout << "d2: " << d2 << std::endl;
    std::cout << "d3: " << d3 << std::endl;
    std::cout << "d1 + d2: " << d1 + d2 << std::endl;
    std::cout << "d1 + d2 + d3: " << d1 + d2 + d3 << std::endl;
    std::cout << "det: " << det << std::endl;

    tempMatrix(_1,_1) =  1 * (inMatrix(_2,_2) * inMatrix(_3,_3) - inMatrix(_2,_3) * inMatrix(_3,_2))/det;
    tempMatrix(_1,_2) = -1 * (inMatrix(_2,_1) * inMatrix(_3,_3) - inMatrix(_2,_3) * inMatrix(_3,_1))/det;
    tempMatrix(_1,_3) =  1 * (inMatrix(_2,_1) * inMatrix(_3,_2) - inMatrix(_2,_2) * inMatrix(_3,_1))/det;
    tempMatrix(_2,_1) = -1 * (inMatrix(_1,_2) * inMatrix(_3,_3) - inMatrix(_1,_3) * inMatrix(_3,_2))/det;
    tempMatrix(_2,_2) =  1 * (inMatrix(_1,_1) * inMatrix(_3,_3) - inMatrix(_1,_3) * inMatrix(_3,_1))/det;
    tempMatrix(_2,_3) = -1 * (inMatrix(_1,_1) * inMatrix(_3,_2) - inMatrix(_1,_2) * inMatrix(_3,_1))/det;
    tempMatrix(_3,_1) =  1 * (inMatrix(_1,_2) * inMatrix(_2,_3) - inMatrix(_1,_3) * inMatrix(_2,_2))/det;
    tempMatrix(_3,_2) = -1 * (inMatrix(_1,_1) * inMatrix(_2,_3) - inMatrix(_1,_3) * inMatrix(_2,_1))/det;
    tempMatrix(_3,_3) =  1 * (inMatrix(_1,_1) * inMatrix(_2,_2) - inMatrix(_1,_2) * inMatrix(_2,_1))/det;
    return tempMatrix;
}
