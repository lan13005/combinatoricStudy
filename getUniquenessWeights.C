// Since the DSelector runs combo by combo we have to save the values and fill them all.
// This would not be hard but just kind of annoying to implement. It is best to just read in the data and count the number of combos that pass the event. This is what we do here.
//
// Some things to worry about:
// 1. What to do with bestChiSq when two combos have the same chiSq? Currently we just take the first one. Probably doesnt matter since they P4s should very likely be the same right?
// 2. The eventNumbers can be repeated in the MC. So we can have two events with the same event number. Hopefully we avoided this by checking against the previous event number
// 3. Program is very convoluted with a lot of conditions. I know.

bool verbose=false;

void addUTWeightsBranch(string rootFileLoc, string rootFileName, string treeName){
    TFile* dataFile = new TFile((rootFileLoc+rootFileName+".root").c_str());
    TTree* dataTree;
    dataFile->GetObject((treeName).c_str(),dataTree);

    Long64_t nentries = (Long64_t)dataTree->GetEntries();
    ULong64_t event;
    double chiSq;
    double DOF;
    // We will basically count the number of combos that have the same eventNumber.
    dataTree->SetBranchAddress("event",&event);
    dataTree->SetBranchAddress("chiSq",&chiSq);
    dataTree->SetBranchAddress("DOFKinFit",&DOF);

    // clone the tree and add two new branches
    TFile *ut_dataFile = TFile::Open((rootFileLoc+rootFileName+"_UTweights.root").c_str(),"RECREATE"); 
    TTree *ut_dataTree = dataTree->CloneTree(-1,"fast"); 
    double UT_equalWeight;
    double UT_bestChiWeight;
    double UT_probWeight;
    TBranch* UT_equalWeights=ut_dataTree->Branch("UT_equalWeights",&UT_equalWeight,"UT_equalWeights/D");
    TBranch* UT_bestChiWeights=ut_dataTree->Branch("UT_bestChiWeights",&UT_bestChiWeight,"UT_bestChiWeights/D");
    TBranch* UT_probWeights=ut_dataTree->Branch("UT_probWeights",&UT_probWeight,"UT_probWeights/D");
    
    // -------------------------------------
    // Count the number of combos in an event
    // Find the best ChiSq in an event and given a weight of 1, 0's to others
    // Find the pvalues and weight the combo accordingly
    // -------------------------------------
    int count=0;
    ULong64_t previousEvent=-1;
    //nentries=50;
    std::vector<double> pValuesInEvent;
    std::vector<double> chiSqsInEvent;
    double pValue;
    double pValueNormalization=0;
    int idx_bestChiSq;
    double bestChiSq=DBL_MAX;

    //nentries=100;//330681;
    cout << "nentries: " << nentries << endl;
    double countedEntries=0;
    for(Long64_t ientry=0; ientry<nentries; ientry++)
    {
    	dataTree->GetEntry(ientry);
        if (ientry==0) previousEvent=event; // have to make a special case for just the first event to get the right initialization
        if (previousEvent==event){
            ++count;
            chiSqsInEvent.push_back(chiSq); // keeping the chiSq for debugging
            pValue=TMath::Prob(chiSq,DOF);
            pValuesInEvent.push_back(pValue);  
            pValueNormalization += pValue;
            if(chiSq<bestChiSq){
                bestChiSq=chiSq;
                idx_bestChiSq=count-1; // since c++ is zero indexed
            }
        }
        if ( (previousEvent!=event) || (ientry==nentries-1) ) { // if its the the last entry then we will fill it even though the event number can be seen before
            if(verbose) cout << "^--Found new event (or we are at the last entry) Filling all previous combos from previous event" << endl; 
            if(ientry>0){ // dont want to start filling on the first event
                for(int iCount=0; iCount<count; ++iCount){
                    UT_equalWeight=1.0/count;
                    UT_probWeight=pValuesInEvent[iCount]/pValueNormalization;
                    if (iCount==idx_bestChiSq){
                        UT_bestChiWeight=1;
                    }
                    else { UT_bestChiWeight=0; }
                    UT_equalWeights->Fill();
                    UT_probWeights->Fill();
                    UT_bestChiWeights->Fill();
                    ++countedEntries;
                    if (verbose) {
                        cout << "event " << previousEvent << " | chiSq " << chiSqsInEvent[iCount] << " | pValue " << pValuesInEvent[iCount] <<
                            " | equalWeighting " << UT_equalWeight << " | bestChiSqWeight " << UT_bestChiWeight << " | probWeighted " << UT_probWeight << endl;
                    }
                }

                // ---------------
                // RESET EVERYTHIN
                // ---------------
                previousEvent=event;
                count=1; // if we are in this condition then the loop has just found a combo with an event not equal to the previous. So it starts at 1
                bestChiSq=DBL_MAX; 
                pValuesInEvent.clear();
                chiSqsInEvent.clear();

                // ---------------
                // Since we are at a new event we have to save the pValues and chiSqs or else we would have skipped them
                // ---------------
                pValue=TMath::Prob(chiSq,DOF);
                pValueNormalization=pValue; 
                if(verbose)cout << "Reintialized pValueNormalization: " << pValueNormalization << endl << endl;
                chiSqsInEvent.push_back(chiSq);
                if(chiSq<bestChiSq){
                    bestChiSq=chiSq;
                    idx_bestChiSq=count-1; // since c++ is zero indexed
                }
                pValuesInEvent.push_back(pValue); 


                // We need this condition here for the case if the the last combo corresponds to a new event since the section above only fills the combos from the previous event. 
                // If this is true then all the weights are equal to 1 since thats the only choice. 
                if (ientry==nentries-1 && previousEvent!=event) {
                    UT_equalWeight=1;
                    UT_probWeight=1;
                    UT_bestChiWeight=1;
                    UT_equalWeights->Fill();
                    UT_probWeights->Fill();
                    UT_bestChiWeights->Fill();
                    ++countedEntries;
                    if (verbose) {
                        cout << "event " << event << " | chiSq " << chiSq << " | pValue " << pValue <<
                            " | equalWeighting " << 1 << " | bestChiSqWeight " << 1 << " | probWeighted " << 1 << endl;
                    }
                }
            }
        }
        if(verbose) cout << "-event" << event << " chiSq: " << chiSq << endl;
    }

    ut_dataFile->cd();
    ut_dataTree->Write(); 

    cout << "\n--------------" << endl;
    cout << "nentries vs counted entries: " << nentries << ", " << countedEntries << endl;
}


void getUniquenessWeights(){
    string rootFileLoc = "/d/grid15/ln16/pi0eta/092419/";
    string rootFileName = "degALL_a0Test_treeFlat_DSelector";
    string treeName = "degALL_a0Test_tree_flat";
    addUTWeightsBranch(rootFileLoc, rootFileName, treeName);

    //string rootFileLoc = "/d/grid15/ln16/pi0eta/q-values/";
    //string rootFileLoc="/home/lawrence/Desktop/gluex/q-values/";
}
