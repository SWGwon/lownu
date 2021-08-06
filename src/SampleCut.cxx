#include "getopt.h"

#include "SampleCut.hxx"

std::string fileName;
//------------------------------------------------------------------------------
void PrintSyntax();
bool ParseArgs(int argc, char* argv[]);
//------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (!ParseArgs(argc, argv)) return 0;

    SampleCut(fileName);

    return 0;
}
//------------------------------------------------------------------------------
bool ParseArgs(int argc, char* argv[]) {
    bool status = false;
    int index;
    int iarg = 0;
    const struct option longopts[] = {
        {"input-file", required_argument, 0, 'i'}, 
        {"tof",        required_argument, 0, '1'}, 
        {"edep",       required_argument, 0, '2'}, 
        {"neighbor-distance",  required_argument, 0, '3'}, 
        {"branch-num",  required_argument, 0, '4'}, 
        {"track-length",     required_argument, 0, '5'}, 
        {"lownu",      required_argument, 0, '6'}, 
        {"max-angle",       required_argument, 0, '7'}, 
        {"max-distance",       required_argument, 0, '8'}, 
        {0,0,0,0}
    };

    while (iarg != -1) {
        iarg = getopt_long(argc, argv, "i:12345678h", longopts, &index);
        switch (iarg) {
            case 'i' : {
                           fileName = optarg;
                           break;
                       }
            case '1' : {
                           tofCut = true;
                           break;
                       }
            case '2' : {
                           eDepCut = true;
                           break;
                       }
            case '3' : {
                           neighborDCut = false;
                           break;
                       }
            case '4' : {
                           numBranchCut = true;
                           break;
                       }
            case '5' : {
                           trackLCut = false;
                           break;
                       }
            case '6' : {
                           lownuCut = false;
                           break;
                       }
            case '7' : {
                           maxACut = true;
                           break;
                       }
            case '8' : {
                           maxDCut = true;
                           break;
                       }
            case 'h' : {
                           PrintSyntax();
                           break;
                       }
        }
    }

    if (!fileName.empty()) {
        status = true;
    }
    if (!status)
        PrintSyntax();
    return status;
}
//------------------------------------------------------------------------------
void PrintSyntax() {
    std::cout << "./SampleCut\n";
    std::cout << "  -i, --input-file filename (REQUIRED)\n";
    std::cout << "  --tof, turn on tof cut, cut value: \n";
    std::cout << "          tof > " << tofCutValue << std::endl;
    std::cout << "  --edep, turn on edep cut, cut value: \n";
    std::cout << "          track eDep > " << eDepCutValueTrack << ", cluster eDep > " << eDepCutValueCluster << std::endl;
    std::cout << "  --neighbor-distance, turn on distance cut, cut value: \n";
    std::cout << "          neighborDistance > " << neighborDCutValue << std::endl;
    std::cout << "  --branch-num, turn on num cut, cut value: \n";
    std::cout << "          number of branches < " << numBranchCutValue << std::endl;
    std::cout << "  --track-length, turn on length cut, cut value: \n";
    std::cout << "          track length > " << trackLCutValue << std::endl;
    std::cout << "  --lownu, turn on lownu cut, cut value: \n";
    std::cout << "          reoc nu < " << lownuCutValue << std::endl;
    std::cout << "  --max-angle, turn on angle cut, cut value: \n";
    std::cout << "          maximum angle < " << maxACutValue << std::endl;
    std::cout << "  --max-distance, turn on distance cut, cut value: \n";
    std::cout << "          maximum distance < " << maxDCutValue << std::endl;
    std::cout << std::endl;

}
//------------------------------------------------------------------------------
void SampleCut(std::string fileName) {
    PrintProgramInfo();

    file = new TFile(fileName.c_str());
    tree = (TTree*)file->Get("tree");

    if (!file->IsOpen() || !tree) {
        PrintBar(barSize);
        std::cout << "invalide input" << std::endl;
        exit(0);
        return;
    }

    outName = fileName;
    for (std::string::iterator it = outName.end(); it != outName.begin(); --it) {
        if (*it == '/') {
            outName.erase(outName.begin(), ++it);
            break;
        }
    }
    for (std::string::iterator it = outName.end(); it != outName.begin(); --it) {
        if (*it == '.') {
            outName.erase(it, outName.end());
            break;
        }
    }
    outName += "AfterCut.root";

    TFile* outFile = new TFile(outName.c_str(), "RECREATE");
    TTree* outTree = tree->CloneTree(0);
    bool isSignal;
    outTree->Branch("isSignal", &isSignal, "isSignal/O");

    SetTree();
    PrintInputInfo();
    PrintCutInfo();

    std::cout << "Applying the cuts" << std::endl;
    for (int i = 0; i < tree->GetEntries(); ++i) {
        PrintProgress(i, tree->GetEntries());

        tree->GetEntry(i);
        isSignal = false;

        if (isnan(recoNu))
            continue;

        bool track = false;
        bool cluster = false;

        if (category == 0 || category == 2)
            track = true;
        if (category == 1 || category == 3)
            cluster = true;

        if ((category == 0 || category == 1) && numberOfFSNeutron == 1 && trueNu < 300)
            isSignal = true;

        if (isSignal)
            totalNumberOfEvents++;

        //0: tof : tof > 0
        if (tofCut && tof < tofCutValue)
            continue;
        if (isSignal)
            sig[0]++;
        else
            bkg[0]++;

        //1: eDep: track + cluster
        if (((track && eDep < eDepCutValueTrack) || (cluster && eDep < eDepCutValueCluster)))
            //if (eDepCut && ((track && eDep < 3700) || (cluster && eDep < 735)))
            continue;
        if (isSignal)
            sig[1]++;
        else
            bkg[1]++;

        //2: neighborDistance > 20
        if (neighborDCut && cluster && neighborDistance < 20)
            continue;
        if (isSignal)
            sig[2]++;
        else
            bkg[2]++;

        //3: numberOfBranches < 1
        if (numBranchCut && track && numberOfBranches > 0)
            continue;
        if (isSignal)
            sig[3]++;
        else
            bkg[3]++;

        //4: trackLength > 20
        if (trackLCut && track && trackLength < 20)
            continue;
        if (isSignal)
            sig[4]++;
        else
            bkg[4]++;

        //5: lownu
        if (lownuCut && recoNu > 200)
            continue;
        if (isSignal)
            sig[5]++;
        else
            bkg[5]++;

        //6: max angle
        if (maxACut && maxAngle > 0.7536)
            continue;
        if (isSignal)
            sig[6]++;
        else
            bkg[6]++;

        //7: max distance
        if (maxDCut && maxDistance > 700)
            continue;
        if (isSignal)
            sig[7]++;
        else {
            bkg[7]++;
        }

        outTree->Fill();
    }
    //std::cout << std::endl;

    PrintResult();

    outFile->Write();
    outFile->Close();
}
//---------------------------------------------------------------------------------------------------------------
void GetCommandLineArgs(int argc, char* argv[])
{
    if (argc < 2) {
        PrintProgramInfo();
        PrintSyntax();
    }
}
//---------------------------------------------------------------------------------------------------------------
void PrintProgramInfo()
{
    //for (int i = 0; i < barSize; ++i) {
    //    std::cout << "=";
    //}
    std::cout << std::setw(barSize) << std::setfill('/') << "" << std::endl;
    std::cout << std::setw((barSize - 10)/2 - 1) << " " << "Sample cut" << " " << std::setw((barSize - 10)/2 ) << std::setfill('/') << std::right<< "/" << std::endl;
    std::cout << std::setw((barSize - 10)/2 - 1) << " " << "2021.08.02" << " " << std::setw((barSize - 10)/2 ) << std::setfill('/') << std::right<< "/" << std::endl;
    std::cout << std::setw((barSize - 11)/2 - 1) << "/" << " Sunwoo Gwon" << " " << std::setw((barSize - 11)/2 ) << std::setfill('/') << std::right<< "/" << std::endl;
    std::cout << std::setw(barSize) << std::setfill('/') << "" << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void SetTree()
{
    tree->SetBranchAddress("category", &category);
    tree->SetBranchAddress("numberOfBranches", &numberOfBranches);
    tree->SetBranchAddress("numberOfFSNeutron", &numberOfFSNeutron);
    tree->SetBranchAddress("trackLength", &trackLength);
    tree->SetBranchAddress("recoNu",&recoNu);
    tree->SetBranchAddress("trueNu",&trueNu);
    tree->SetBranchAddress("tof",&tof);
    tree->SetBranchAddress("maxAngle",&maxAngle);
    tree->SetBranchAddress("maxDistance", &maxDistance);
    tree->SetBranchAddress("neighborDistance", &neighborDistance);
    tree->SetBranchAddress("eDep", &eDep);
}
//---------------------------------------------------------------------------------------------------------------
void PrintInputInfo()
{
    PrintBar(barSize);
    std::cout << "Input info" << std::endl;
    std::cout << "|-total sample size: " << tree->GetEntries() << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void PrintCutInfo()
{
    const int width = 15;
    PrintBar(barSize);
    std::cout << "Cuts info" << std::endl;

    std::cout << std::left << std::setw(width) << std::setfill(' ')  << "|-tofCut "<< "- ";
    if (tofCut) {
        std::cout << "ON: tof > " << tofCutValue << std::endl;
    } else {
        std::cout << "OFF" << std::endl;
    }

    std::cout << std::left << std::setw(width) << std::setfill(' ') << "|-eDepCut " << "- ";
    if (eDepCut) {
        std::cout << "ON: track eDep > " << eDepCutValueTrack << ", cluster eDep > " << eDepCutValueCluster << std::endl;
    } else {
        std::cout << "OFF" << std::endl;
    }

    std::cout.setf(std::ios::left);
    std::cout << std::setw(width) << "|-neighborDCut " << std::setfill(' ') << "- ";
    if (neighborDCut) {
        std::cout << "ON: neighborDistance > " << neighborDCutValue << std::endl;
    } else {
        std::cout << "OFF" << std::endl;
    }
    
    std::cout.setf(std::ios::left);
    std::cout << std::setw(width) << "|-numBranchCut " << std::setfill(' ') << "- ";
    if (numBranchCut) {
        std::cout << "ON: number of branches < " << numBranchCutValue << std::endl;
    } else {
        std::cout << "OFF" << std::endl;
    }
    
    std::cout.setf(std::ios::left);
    std::cout << std::setw(width) << "|-trackLCut " << std::setfill(' ') << "- ";
    if (trackLCut) {
        std::cout << "ON: track length > " << trackLCutValue << std::endl;
    } else {
        std::cout << "OFF" << std::endl;
    }
    
    std::cout.setf(std::ios::left);
    std::cout << std::setw(width) << "|-lownuCut " << std::setfill(' ') << "- ";
    if (lownuCut) {
        std::cout << "ON: reoc nu < " << lownuCutValue << std::endl;
    } else {
        std::cout << "OFF" << std::endl;
    }
    
    std::cout.setf(std::ios::left);
    std::cout << std::setw(width) << "|-maxACut " << std::setfill(' ') << "- ";
    if (maxACut) {
        std::cout << "ON: maximum angle < " << maxACutValue << std::endl;
    } else {
        std::cout << "OFF" << std::endl;
    }
    
    std::cout.setf(std::ios::left);
    std::cout << std::setw(width) << "|-maxDCut " << std::setfill(' ') << "- ";
    if (maxDCut) {
        std::cout << "ON: maximum distance < " << maxDCutValue << std::endl;
    } else {
        std::cout << "OFF" << std::endl;
    }

    PrintBar(barSize);
}
//---------------------------------------------------------------------------------------------------------------
void PrintProgress(int i, int total)
{
    if (i == 0)
        std::cout << "Processing..." << std::endl;
    float progress = (double)(i+1)/total;
    const int barWidth = barSize - 8;
    std::cout << "|-[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();

    if (int(progress * 100.0) >= 100) {
        std::cout << std::endl;
        std::cout << "|-done " << std::endl;
    }
}
//---------------------------------------------------------------------------------------------------------------
void PrintResult()
{
    //for (int j = 0; j < 8; ++j) {
    //    std::cout << "sig[" << j << "]: " << sig[j] << std::endl;
    //    std::cout << "bkg[" << j << "]: " <<  bkg[j] << std::endl;
    //    std::cout << "purity[" << j << "]: " << sig[j] / (sig[j] + bkg[j]) << std::endl;
    //    std::cout << "efficiency[" << j << "]: " << sig[j] / totalNumberOfEvents << std::endl;
    //}

    PrintBar(barSize);
    std::cout << "Result Info" << std::endl;
    std::cout << "|-reduced sample size: " << sig[7] << std::endl;
    std::cout << "|-signal cut efficiency: " << sig[7] / totalNumberOfEvents << std::endl;
    std::cout << "|-purity: " << sig[7] / (sig[7] + bkg[7]) << std::endl;
    std::cout << "|-output: " << outName + "AfterCut.root" << std::endl;
    PrintBar(barSize);
}
//---------------------------------------------------------------------------------------------------------------
void PrintBar(int i)
{
    std::cout << std::setw(i) << std::setfill('-') << "" << std::endl;
}
