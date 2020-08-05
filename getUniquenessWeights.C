// Since the DSelector runs combo by combo we have to save the values and fill them all.
// This would not be hard but just kind of annoying to implement. It is best to just read in the data and count the number of combos that pass the event. This is what we do here.
//
// Some things to worry about:
// 1. What to do with bestChiSq when two combos have the same chiSq? Currently we just take the first one. Probably doesnt matter since they P4s should very likely be the same right?
// 2. The eventNumbers can be repeated in the MC. So we can have two events with the same event number. Hopefully we avoided this by checking against the previous event number
// 3. Program is very convoluted with a lot of conditions. I know.

bool verbose=true;

void addUTWeightsBranch(string rootFileLoc, string rootFileName, string treeName){
    TFile* dataFile = new TFile((rootFileLoc+rootFileName+".root").c_str());
    TTree* dataTree;
    dataFile->GetObject((treeName).c_str(),dataTree);

    Long64_t nentries = (Long64_t)dataTree->GetEntries();
    ULong64_t event;
    double chiSq;
    double DOF;
    Int_t numSpect;
    Int_t spectroscopicIDs[20]; // easy way to store as much final state particle IDs we care about. 
    Int_t beamID;
    double Mpi0eta;
    TLorentzVector* beam_p4_meas=NULL;
    double beamEnergy;
    vector<double> beamEnergies;
    // We will basically count the number of combos that have the same eventNumber.
    dataTree->SetBranchAddress("event",&event);
    dataTree->SetBranchAddress("chiSq",&chiSq);
    dataTree->SetBranchAddress("DOFKinFit",&DOF);
    dataTree->SetBranchAddress("numSpect",&numSpect);
    dataTree->SetBranchAddress("spectroscopicID",&spectroscopicIDs);
    dataTree->SetBranchAddress("beamID",&beamID);
    dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
    dataTree->SetBranchAddress("beam_p4_meas",&beam_p4_meas);

    // clone the tree and add two new branches
    TFile *ut_dataFile = TFile::Open((rootFileLoc+rootFileName+"_UTweights.root").c_str(),"RECREATE"); 
    TTree *ut_dataTree = dataTree->CloneTree(-1,"fast"); 
    double UT_equalWeight;
    double UT_spectEqualWeight;
    double UT_bestChiWeight;
    double UT_spectBestChiWeight;
    double UT_probWeight;
    double UT_spectProbWeight;
    string spectroscopicID;
    int UT_uniqueEventNumber=0; // NEED TO UNSCRAMBLE THE EVENT NUMBERS. THERE ARE MULTIPLE EVENTS WITH THE SAME EVENT NUMBER IN THE ROOT TREES. NOT ENTIRELY SURE WHAT THE SOURCE IS. PROBABLY MULTIPLE SIMULATIONS?
    TBranch* UT_equalWeights=ut_dataTree->Branch("UT_equalWeights",&UT_equalWeight,"UT_equalWeights/D");
    TBranch* UT_spectEqualWeights=ut_dataTree->Branch("UT_spectEqualWeights",&UT_spectEqualWeight,"UT_spectEqualWeights/D");
    TBranch* UT_bestChiWeights=ut_dataTree->Branch("UT_bestChiWeights",&UT_bestChiWeight,"UT_bestChiWeights/D");
    TBranch* UT_spectBestChiWeights=ut_dataTree->Branch("UT_spectBestChiWeights",&UT_spectBestChiWeight,"UT_spectBestChiWeights/D");
    TBranch* UT_probWeights=ut_dataTree->Branch("UT_probWeights",&UT_probWeight,"UT_probWeights/D");
    TBranch* UT_spectProbWeights=ut_dataTree->Branch("UT_spectProbWeights",&UT_spectProbWeight,"UT_spectProbWeights/D");
    TBranch* UT_uniqueEventNumbers=ut_dataTree->Branch("UT_uniqueEventNumbers",&UT_uniqueEventNumber,"UT_uniqueEventNumbers/I");
    
    // -------------------------------------
    // Count the number of combos in an event
    // Find the best ChiSq in an event and given a weight of 1, 0's to others
    // Find the pvalues and weight the combo accordingly
    // -------------------------------------
    int count=0;
    int beamCount=0;
    ULong64_t previousEvent=-1;
    //nentries=50;
    std::vector<double> pValuesInEvent;
    std::vector<double> chiSqsInEvent;
    std::vector<string> spectroscopicIDVector;
    std::vector<Int_t> beamIDs;
    std::vector<double> Mpi0etas;
    double pValue;
    vector<double> pValueNormalizations;
    double pValueNormalization;
    int idx_bestChiSq;
    double bestChiSq=DBL_MAX;
    // I will set a maximum that there can only be up to 20 unique spectroscopic combos
    vector<double> bestChiSqs;
    vector<int> idx_bestChiSqs; 
    vector<int> spectCounts;

    //nentries=1000;//330681;
    cout << "nentries: " << nentries << endl;
    double countedEntries=0; 
    map<int,int> alreadySeenBeamID; // grouping into combos with the same beamID, i.e. spectroscopicID vary, and look for the best chiSq
    set<string> alreadySeenSpectCombo; // counting how many spectroscopicID there are
    for(Long64_t ientry=0; ientry<nentries; ientry++)
    {
    	dataTree->GetEntry(ientry);
        beamEnergy=beam_p4_meas->E();
        if(verbose)cout << "(previousEvent "<< previousEvent << ")currentEvent " << event << " EBeam " << beamEnergy << " BeamID " << beamID << " spectroscopicID " << spectroscopicID << " chiSq " << chiSq <<  endl;
        if (ientry==0) previousEvent=event; // have to make a special case for just the first event to get the right initialization
        if (previousEvent==event){
            cout << "-previousEvent same as currentEvent... Filling" << endl;
            ++count; // counts total number of combos in an event

            // --------------------
            // Make a string to hold the spectroscopic (photon) ids so we can do some checks
            spectroscopicID="";
            for (int iSpect=0; iSpect<numSpect; ++iSpect){
                spectroscopicID+=to_string(spectroscopicIDs[iSpect]);
            }

            // --------------------
            // Each beam photon has its own indexer that holds the index of the best chiSq
            if(alreadySeenBeamID.find(beamID)==alreadySeenBeamID.end()){
                alreadySeenBeamID.insert(make_pair(beamID,beamCount));
                bestChiSqs.push_back(DBL_MAX);
                idx_bestChiSqs.push_back(-1); // will overwrite this, just using -1 as a place holder
                ++beamCount;

               spectCounts.push_back(1);  // starts at 1 since we just saw a new beam photon

               pValueNormalizations.push_back(0);
            }
            else{
                // I wont make an explicit check that the spectroscopic combo has not been seen before.
                spectCounts[alreadySeenBeamID[beamID]] += 1; 
            }
            
            // Need to create vectors to hold the informaton for the event. Once we see a new event we will begin filling the branches with these vectors
            Mpi0etas.push_back(Mpi0eta);
            spectroscopicIDVector.push_back(spectroscopicID);
            beamIDs.push_back(beamID);
            beamEnergies.push_back(beamEnergy);
            chiSqsInEvent.push_back(chiSq); // keeping the chiSq for debugging
            pValue=TMath::Prob(chiSq,DOF);
            pValuesInEvent.push_back(pValue);  
            pValueNormalizations[alreadySeenBeamID[beamID]]+= pValue; // need a normalization factor, Sum(p_i) 
            pValueNormalization+=pValue;
            if(chiSq<bestChiSq){
                bestChiSq=chiSq;
                idx_bestChiSq=count-1; // since c++ is zero indexed
            }
            int idxCombosBestChiSq=alreadySeenBeamID[beamID];
            if(chiSq<bestChiSqs[idxCombosBestChiSq]){ // compare chiSq to bestChiSq for the spectroscopic combo
                bestChiSqs[idxCombosBestChiSq]=chiSq;
                idx_bestChiSqs[idxCombosBestChiSq]=count-1;
            }
        }
        
        // --------------------------
        // Once we see a new event we will fill the data for the now previous event.
        // We will be on a new event so we have to include this information in the fill vectors the next event
        // You might notice that if the last event was also a new event then we will not fill it in this loop. We made a special condition for it after this loop
        if (previousEvent!=event || (ientry==(nentries-1))) {
            ++UT_uniqueEventNumber;
            if(verbose) cout << "^--Found new event (or we are at the last entry) Filling all previous combos from previous event" << endl; 
            if(ientry>0){ // dont want to start filling on the first event
                double sumSpectEqualWeights=0;
                double sumBestChiWeights=0;
                for(int iCount=0; iCount<count; ++iCount){
                    UT_equalWeight=1.0/count;
                    UT_spectEqualWeight=1.0/spectCounts[alreadySeenBeamID[beamIDs[iCount]]]; // pretty convoluted but we get the beamID for that combo and get the index for the spectCounts
                    UT_probWeight=pValuesInEvent[iCount]/pValueNormalization;
                    UT_spectProbWeight=pValuesInEvent[iCount]/pValueNormalizations[alreadySeenBeamID[beamIDs[iCount]]];
                    if (iCount==idx_bestChiSq){
                        UT_bestChiWeight=1;
                    }
                    else { UT_bestChiWeight=0; }
                    if ( find(begin(idx_bestChiSqs), end(idx_bestChiSqs), iCount) != end(idx_bestChiSqs) ){
                        UT_spectBestChiWeight=1;
                    }
                    else { UT_spectBestChiWeight=0; }
                    sumSpectEqualWeights+=UT_spectEqualWeight;
                    sumBestChiWeights+=UT_spectBestChiWeight;
                    UT_equalWeights->Fill();
                    UT_probWeights->Fill();
                    UT_spectProbWeights->Fill();
                    UT_bestChiWeights->Fill();
                    UT_spectEqualWeights->Fill();
                    UT_spectBestChiWeights->Fill();
                    UT_uniqueEventNumbers->Fill();
                    ++countedEntries;
                    if (verbose) {
                        cout << "event " << previousEvent << " | Mpi0eta " << Mpi0etas[iCount] << " | beamEnergy " << beamEnergies[iCount] <<
                            " | beamID " << beamIDs[iCount] <<  " | spectroscopicID " << spectroscopicIDVector[iCount] << 
                            " | chiSq " << chiSqsInEvent[iCount] << " | pValue " << pValuesInEvent[iCount] <<
                            " | equalWeight " << UT_equalWeight << " | bestChiSqWeight " << UT_bestChiWeight << " | probWeight " << UT_probWeight << 
                            " | spectEqualWeight " << UT_spectEqualWeight << " | spectBestChiSqWeight " << UT_spectBestChiWeight << 
                            " | spectProbWeight " << UT_spectProbWeight <<
                            " | uniqueEventNumber " << UT_uniqueEventNumber << endl;
                    }
                }
                if (abs(sumBestChiWeights-sumSpectEqualWeights)>0.05){
                    cout << "\nMismatch of summed spect weights: sumBestChiWeights=" << sumBestChiWeights << " sumSpectEqualWeights=" << sumSpectEqualWeights << endl;
                    exit(0);
                    //for(int iCount=0; iCount<count; ++iCount){
                    //    UT_equalWeight=1.0/count;
                    //    UT_spectEqualWeight=1.0/spectCounts[alreadySeenBeamID[beamIDs[iCount]]]; // pretty convoluted but we get the beamID for that combo and get the index for the spectCounts;
                    //    UT_probWeight=pValuesInEvent[iCount]/pValueNormalization;
                    //    if (iCount==idx_bestChiSq){
                    //        UT_bestChiWeight=1;
                    //    }
                    //    else { UT_bestChiWeight=0; }
                    //    if ( find(begin(idx_bestChiSqs), end(idx_bestChiSqs), iCount) != end(idx_bestChiSqs) ){
                    //        UT_spectBestChiWeight=1;
                    //    }
                    //    else { UT_spectBestChiWeight=0; }
                    //    cout << "event " << previousEvent << " | Mpi0eta " << Mpi0etas[iCount] << " | beamID " << beamIDs[iCount] <<  " | spectroscopicID " << spectroscopicIDVector[iCount] << 
                    //        " | chiSq " << chiSqsInEvent[iCount] << " | pValue " << pValuesInEvent[iCount] <<
                    //        " | equalWeight " << UT_equalWeight << " | bestChiSqWeight " << UT_bestChiWeight << " | probWeight " << UT_probWeight << 
                    //        " | spectEqualWeight " << UT_spectEqualWeight << " | spectBestChiSqWeight " << UT_spectBestChiWeight << " | uniqueEventNumber " << UT_uniqueEventNumber << endl;
                    //}
                }

                
                // ---------------
                // RESET EVERYTHIN
                // ---------------
                // We will reset everything now. Will not update previousEvent if we are on the last entry because we will use it to determine if the last event is also new. This requires us to manually fill
                if(ientry!=(nentries-1)){
                    previousEvent=event;
                }
                if(verbose) cout << "-- Resetting everything\n" << endl;
                count=1; // if we are in this condition then the loop has just found a combo with an event not equal to the previous. So it starts at 1
                bestChiSq=DBL_MAX; 
                pValuesInEvent.clear();
                chiSqsInEvent.clear();
                spectroscopicIDVector.clear();
                beamIDs.clear();
                Mpi0etas.clear();
                idx_bestChiSqs.clear();
                bestChiSqs.clear();
                spectCounts.clear();
                pValueNormalizations.clear();
                beamEnergies.clear();

                // ---------------
                // Since we are at a new event we have to save the pValues and chiSqs or else we would have skipped them
                // ---------------
                pValue=TMath::Prob(chiSq,DOF);
                chiSqsInEvent.push_back(chiSq);
                
                // Need to set spectroscopicID for this new event
                spectroscopicID="";
                beamCount=0;
                alreadySeenBeamID.clear();
                alreadySeenSpectCombo.clear();
                for (int iSpect=0; iSpect<numSpect; ++iSpect){
                    spectroscopicID+=to_string(spectroscopicIDs[iSpect]);
                }
                spectroscopicIDVector.push_back(spectroscopicID);
                beamIDs.push_back(beamID);
                Mpi0etas.push_back(Mpi0eta);
                beamEnergy=beam_p4_meas->E();
                beamEnergies.push_back(beamEnergy);

                // Need to set the new bestChiSq and its index into the combos array
                alreadySeenBeamID.insert(make_pair(beamID,beamCount));
                alreadySeenSpectCombo.insert(spectroscopicID);
                spectCounts.push_back(1);
                ++beamCount;
                bestChiSq=chiSq;
                bestChiSqs.push_back(chiSq);
                idx_bestChiSq=count-1; // since c++ is zero indexed
                idx_bestChiSqs.push_back(idx_bestChiSq); 
                
                pValuesInEvent.push_back(pValue); 
                pValueNormalizations.push_back(pValue); // the normalization for this new event, which automatically has a new beam photon, is the current pValue.
                pValueNormalization=pValue;
            }
        }
        // We need this condition here for the case if the the last combo corresponds to a new event since the section above only fills the combos from the previous event. 
        // If this is true then all the weights are equal to 1 since thats the only choice. 
        // Also, since this is the last event we dont have to worry about resetting anything
        if (previousEvent!=event && (ientry==nentries-1)) {
            ++UT_uniqueEventNumber;
            ++countedEntries;
            UT_equalWeight=1;
            UT_probWeight=1;
            UT_bestChiWeight=1;
            UT_spectEqualWeight=1;
            UT_spectBestChiWeight=1;
            UT_spectProbWeight=1;
            UT_equalWeights->Fill();
            UT_probWeights->Fill();
            UT_bestChiWeights->Fill();
            UT_spectEqualWeights->Fill();
            UT_spectBestChiWeights->Fill();
            UT_spectProbWeights->Fill();
            UT_uniqueEventNumbers->Fill();
            if (verbose) {
                cout << "(LastEventIsNew) event " << event << " | Mpi0eta " << Mpi0eta << " | beamEnergy " << beamEnergy << 
                    " | beamID " << beamID <<  " | spectroscopicID " << spectroscopicID << 
                    " | chiSq " << chiSq << " | pValue " << pValue <<
                    " | equalWeight 1 | bestChiSqWeight 1 | probWeight 1" << 
                    " | spectEqualWeight 1 | spectBestChiSqWeight 1" << 
                    " | spectProbWeight 1" <<
                    " | uniqueEventNumber " << UT_uniqueEventNumber << endl;
            }
        }
    }
    ut_dataFile->cd();
    ut_dataTree->Write(); 

    cout << "nentries vs counted entries: " << nentries << ", " << countedEntries << endl;
    cout << "\n--------------" << endl;
}


void getUniquenessWeights(){
    string rootFileLoc = "/d/grid15/ln16/pi0eta/092419/";
    string rootFileName = "degALL_a0a2Test_treeFlat_DSelector";
    string treeName = "degALL_a0a2Test_tree_flat";
    addUTWeightsBranch(rootFileLoc, rootFileName, treeName);

    //string rootFileLoc = "/d/grid15/ln16/pi0eta/q-values/";
    //string rootFileLoc="/home/lawrence/Desktop/gluex/q-values/";
}
