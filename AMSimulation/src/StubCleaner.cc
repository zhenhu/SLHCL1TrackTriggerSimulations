#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/StubCleaner.h"

static const unsigned MIN_NGOODSTUBS = 3;
static const unsigned MAX_NGOODSTUBS = 8;

bool sortByFloat(const std::pair<unsigned, float>& lhs, const std::pair<unsigned, float>& rhs) {
    return lhs.second < rhs.second;
}

template<typename RandomAccessIterator, typename Size, typename T>
void insertSorted(RandomAccessIterator first, Size pos, Size len, const T value) {
    first += len;
    len -= pos;
    while (len>0) {
        *first = std::move(*(first-1));
        --len;
        --first;
    }
    *first = std::move(value);
}


// _____________________________________________________________________________
// Read the input ntuples
int StubCleaner::readFile(TString src) {
    if (src.EndsWith(".root")) {
        if (verbose_)  std::cout << Info() << "Opening " << src << std::endl;
        if (chain_->Add(src) )  // 1 if successful, 0 otherwise
            return 0;

    } else if (src.EndsWith(".txt")) {
        TFileCollection* fc = new TFileCollection("fileinfolist", "", src);
        if (chain_->AddFileInfoList((TCollection*) fc->GetList()) )  // 1 if successful, 0 otherwise
            return 0;
    }

    std::cout << Error() << "Input source should be either a .root file or a .txt file." << std::endl;
    return 1;
}


// _____________________________________________________________________________
// Select one stub per layer, one hit per stub.
// If an event fails selection, all stubs are removed
int StubCleaner::cleanStubs(TString out) {
    if (!out.EndsWith(".root")) {
        std::cout << Error() << "Output filename must be .root" << std::endl;
        return 1;
    }

    // Do not do GetEntries() as it conflicts with CloneTree() when reading from XROOTD
    //long long nentries = chain_->GetEntries();
    //if (nentries <= 0) {
    //    std::cout << Error() << "Input source has zero entry." << std::endl;
    //    return 1;
    //}
    //if (nEvents_ > nentries)
    //    nEvents_ = nentries;

    // For reading
    if (verbose_)  std::cout << Info() << "Reading " << nEvents_ << " events; recreating " << out << " after stub cleaning." << std::endl;

    //std::vector<float> *          vb_x          = 0;
    //std::vector<float> *          vb_y          = 0;
    std::vector<float> *          vb_z          = 0;
    std::vector<float> *          vb_r          = 0;
    std::vector<float> *          vb_eta        = 0;
    std::vector<float> *          vb_phi        = 0;
    std::vector<float> *          vb_coordx     = 0;
    std::vector<float> *          vb_coordy     = 0;
    //std::vector<float> *          vb_roughPt    = 0;
    //std::vector<float> *          vb_trigBend   = 0;
    std::vector<unsigned> *       vb_modId      = 0;
    std::vector<int> *            vb_trkId      = 0;
    std::vector<float> *          vp_pt         = 0;
    std::vector<float> *          vp_eta        = 0;
    std::vector<float> *          vp_phi        = 0;
    std::vector<float> *          vp_vx         = 0;
    std::vector<float> *          vp_vy         = 0;
    std::vector<float> *          vp_vz         = 0;
    std::vector<int> *            vp_charge     = 0;

    chain_->SetBranchStatus("*"                 , 0);
    //chain_->SetBranchStatus("TTStubs_x"         , 1);
    //chain_->SetBranchStatus("TTStubs_y"         , 1);
    chain_->SetBranchStatus("TTStubs_z"         , 1);
    chain_->SetBranchStatus("TTStubs_r"         , 1);
    chain_->SetBranchStatus("TTStubs_eta"       , 1);
    chain_->SetBranchStatus("TTStubs_phi"       , 1);
    chain_->SetBranchStatus("TTStubs_coordx"    , 1);
    chain_->SetBranchStatus("TTStubs_coordy"    , 1);
    //chain_->SetBranchStatus("TTStubs_roughPt"   , 1);
    //chain_->SetBranchStatus("TTStubs_trigBend"  , 1);
    chain_->SetBranchStatus("TTStubs_modId"     , 1);
    chain_->SetBranchStatus("TTStubs_trkId"     , 1);
    chain_->SetBranchStatus("genParts_pt"       , 1);
    chain_->SetBranchStatus("genParts_eta"      , 1);
    chain_->SetBranchStatus("genParts_phi"      , 1);
    chain_->SetBranchStatus("genParts_vx"       , 1);
    chain_->SetBranchStatus("genParts_vy"       , 1);
    chain_->SetBranchStatus("genParts_vz"       , 1);
    chain_->SetBranchStatus("genParts_charge"   , 1);

    //chain_->SetBranchAddress("TTStubs_x"        , &(vb_x));
    //chain_->SetBranchAddress("TTStubs_y"        , &(vb_y));
    chain_->SetBranchAddress("TTStubs_z"        , &(vb_z));
    chain_->SetBranchAddress("TTStubs_r"        , &(vb_r));
    chain_->SetBranchAddress("TTStubs_eta"      , &(vb_eta));
    chain_->SetBranchAddress("TTStubs_phi"      , &(vb_phi));
    chain_->SetBranchAddress("TTStubs_coordx"   , &(vb_coordx));
    chain_->SetBranchAddress("TTStubs_coordy"   , &(vb_coordy));
    //chain_->SetBranchAddress("TTStubs_roughPt"  , &(vb_roughPt));
    //chain_->SetBranchAddress("TTStubs_trigBend" , &(vb_trigBend));
    chain_->SetBranchAddress("TTStubs_modId"    , &(vb_modId));
    chain_->SetBranchAddress("TTStubs_trkId"    , &(vb_trkId));
    chain_->SetBranchAddress("genParts_pt"      , &(vp_pt));
    chain_->SetBranchAddress("genParts_eta"     , &(vp_eta));
    chain_->SetBranchAddress("genParts_phi"     , &(vp_phi));
    chain_->SetBranchAddress("genParts_vx"      , &(vp_vx));
    chain_->SetBranchAddress("genParts_vy"      , &(vp_vy));
    chain_->SetBranchAddress("genParts_vz"      , &(vp_vz));
    chain_->SetBranchAddress("genParts_charge"  , &(vp_charge));

    // Set up TTreeFormula for event selection
    TTreeFormula* ttf_event = new TTreeFormula("ttf_event", eventSelect_, chain_);
    //ttf_event->SetQuickLoad(1);

    // For writing
    TFile* tfile = TFile::Open(out, "RECREATE");
    tfile->mkdir("ntupler")->cd();

    TTree* ttree = (TTree*) chain_->CloneTree(0); // Do not copy the data yet
    // The clone should not delete any shared i/o buffers.
    ResetDeleteBranches(ttree);


    // _________________________________________________________________________
    // Loop over all events

    const int good_trkId    =  1;
    const int unmatch_trkId = -1;

    int curTree = chain_->GetTreeNumber();
    int nPassed = 0, nKept = 0;
    unsigned ievt_step = 0;
    for (long long ievt=0; ievt<nEvents_; ++ievt, ++ievt_step) {
        Long64_t local_entry = chain_->LoadTree(ievt);  // for TChain
        if (local_entry < 0)  break;
        if (chain_->GetTreeNumber() != curTree) {  // for TTreeFormula
            curTree = chain_->GetTreeNumber();
            ttf_event->UpdateFormulaLeaves();
        }
        chain_->GetEntry(ievt);

        const unsigned nstubs = vb_modId->size();
        if (verbose_>1 && ievt_step == 50000) {
            std::cout << Debug() << Form("... Processing event: %7lld, keeping: %7i, passing: %7i", ievt, nKept, nPassed) << std::endl;
            ievt_step -= 50000;
        }
        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # stubs: " << nstubs << std::endl;

        if (!nstubs) {  // skip if no stub
            ++nKept;
            ttree->Fill();
            continue;
        }

        if (nstubs > 100000) {
            std::cout << Error() << "Way too many stubs: " << nstubs << std::endl;
            return 1;
        }

        // _____________________________________________________________________
        // Start cleaning
        bool keep = true;

        // Check event selection
        int ndata = ttf_event->GetNdata();
        if (ndata && filter_)
            keep = (ttf_event->EvalInstance() != 0);

        // Check min # of layers
        bool require = (nstubs >= MIN_NGOODSTUBS);
        if (!require && filter_)
            keep = false;

        // Check sim info
        assert(vp_pt->size() == 1);
        float simPt    = vp_pt->front();
        float simEta   = vp_eta->front();
        float simPhi   = vp_phi->front();
        float simTheta = simPt > 0.0 ? (2.0 * std::atan(std::exp(-simEta)) ) : (simEta >= 0 ? 0 : M_PI);

        // Apply pt, eta, phi requirements
        bool sim = (po.minPt  <= simPt  && simPt  <= po.maxPt  &&
                    po.minEta <= simEta && simEta <= po.maxEta &&
                    po.minPhi <= simPhi && simPhi <= po.maxPhi);
        if (!sim && filter_)
            keep = false;

        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " simPt: " << simPt << " simEta: " << simEta << " simPhi: " << simPhi << " keep? " << keep << std::endl;

        // _____________________________________________________________________
        // Remove multiple stubs in one layer
        std::vector<std::pair<unsigned, float> > vec_index_dZ;
        for (unsigned l=0; (l<nstubs) && keep; ++l) {
            int trkId = vb_trkId->at(l);  // check sim info
            if (trkId == good_trkId || trkId == unmatch_trkId) {  // also catch stubs that fail to find a matched simTrack
                //float stub_eta = vb_eta->at(l);
                //float dEta = std::abs(simEta - stub_eta);
                //if (dEta > 0.8)  // way too far
                //    continue;
                float dZ = std::abs(vb_r->at(l) / std::tan(simTheta) - vb_z->at(l));
                vec_index_dZ.push_back(std::make_pair(l, dZ));
            }
        }

        // Sort: smallest dZ to largest
        std::sort(vec_index_dZ.begin(), vec_index_dZ.end(), sortByFloat);

        // Select only one stub per layer
        std::vector<unsigned> goodLayerStubs(16, 999999);
        {   // scoped to limit the scopes of variables
            id_type moduleId, lay, lay16;
            unsigned l;
            float dZ;
            for (unsigned ll=0; (ll<vec_index_dZ.size()) && keep; ++ll) {
                l  = vec_index_dZ.at(ll).first;
                dZ = vec_index_dZ.at(ll).second;

                moduleId = vb_modId->at(l);
                lay      = decodeLayer(moduleId);
                lay16    = compressLayer(lay);
                //assert(lay16 < 16);

                // For each layer, takes the stub with min dZ to simTrack
                if (goodLayerStubs.at(lay16) == 999999 && dZ < 26.0) {  // CUIDADO: gets rid of stubs due to loopers
                    goodLayerStubs.at(lay16) = l;
                }
            }
        }

        //if (keep && goodLayerStubs.at(0) == 999999)
        //    std::cout << Warning() << "... evt: " << ievt << " no stub in the first layer of the barrel!" << std::endl;
        //if (keep && goodLayerStubs.at(6) != 999999 && goodLayerStubs.at(11) != 999999)
        //    std::cout << Warning() << "... evt: " << ievt << " found stubs in the first layers of both positive and negative endcaps!" << std::endl;


        // _____________________________________________________________________
        // Now make keep-or-ignore decision per stub
        unsigned ngoodstubs = 0;
        id_type moduleId;
        for (unsigned l=0; (l<nstubs) && keep; ++l) {
            bool keepstub = true;

            moduleId = vb_modId->at(l);
            if (verbose_>2)  std::cout << Debug() << "... ... stub: " << l << " moduleId: " << moduleId << " trkId: " << vb_trkId->at(l) << std::endl;

            // Check whether the index l was stored as a good stub
            unsigned count = std::count(goodLayerStubs.begin(), goodLayerStubs.end(), l);
            if (!count)  // do not keep even if filter_ is false
                keepstub = false;

            if (keepstub) {
                // Keep and do something similar to insertion sort
                // First, find the correct position to insert (according to moduleId)
                std::vector<unsigned>::const_iterator pos = vb_modId->begin()+ngoodstubs;
                unsigned ipos = ngoodstubs;
                while (ipos != 0 && moduleId < *(--pos))  // start comparing from the tail
                    --ipos;

                // Insert, keeping only the 'ngoodstubs' elements
                //insertSorted(vb_x->begin()        , ipos, ngoodstubs, vb_x->at(l));
                //insertSorted(vb_y->begin()        , ipos, ngoodstubs, vb_y->at(l));
                insertSorted(vb_z->begin()        , ipos, ngoodstubs, vb_z->at(l));
                insertSorted(vb_r->begin()        , ipos, ngoodstubs, vb_r->at(l));
                insertSorted(vb_eta->begin()      , ipos, ngoodstubs, vb_eta->at(l));
                insertSorted(vb_phi->begin()      , ipos, ngoodstubs, vb_phi->at(l));
                insertSorted(vb_coordx->begin()   , ipos, ngoodstubs, vb_coordx->at(l));
                insertSorted(vb_coordy->begin()   , ipos, ngoodstubs, vb_coordy->at(l));
                //insertSorted(vb_roughPt->begin()  , ipos, ngoodstubs, vb_roughPt->at(l));
                //insertSorted(vb_trigBend->begin() , ipos, ngoodstubs, vb_trigBend->at(l));
                insertSorted(vb_modId->begin()    , ipos, ngoodstubs, vb_modId->at(l));
                insertSorted(vb_trkId->begin()    , ipos, ngoodstubs, vb_trkId->at(l));

                ++ngoodstubs;  // remember to increment
            }
            if (verbose_>2)  std::cout << Debug() << "... ... stub: " << l << " keep? " << keepstub << std::endl;
        }
        assert(ngoodstubs <= nstubs);

        // _____________________________________________________________________
        // Now make keep-or-ignore decision per event

        // Check again min # of layers
        require = (ngoodstubs >= MIN_NGOODSTUBS);
        if (!require && filter_)
            keep = false;

        if (keep)
            ++nPassed;
        else //  do not keep any stub
            ngoodstubs = 0;

        if (keep && ngoodstubs > MAX_NGOODSTUBS) {
            std::cout << Warning() << "... evt: " << ievt << " simPt: " << simPt << " simEta: " << simEta << " simPhi: " << simPhi <<  " ngoodstubs: " << ngoodstubs << std::endl;
        }

        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # good stubs: " << ngoodstubs << " keep? " << keep << std::endl;

        //vb_x        ->resize(ngoodstubs);
        //vb_y        ->resize(ngoodstubs);
        vb_z        ->resize(ngoodstubs);
        vb_r        ->resize(ngoodstubs);
        vb_eta      ->resize(ngoodstubs);
        vb_phi      ->resize(ngoodstubs);
        vb_coordx   ->resize(ngoodstubs);
        vb_coordy   ->resize(ngoodstubs);
        //vb_roughPt  ->resize(ngoodstubs);
        //vb_trigBend ->resize(ngoodstubs);
        vb_modId    ->resize(ngoodstubs);
        vb_trkId    ->resize(ngoodstubs);

        ++nKept;
        ttree->Fill();
    }
    if (verbose_)  std::cout << Info() << "Processed and kept " << nKept << " events, passed " << nPassed << std::endl;

    assert(ttree->GetEntries() == nKept);
    tfile->Write();
    delete ttf_event;
    delete ttree;
    delete tfile;

    return 0;
}


// _____________________________________________________________________________
// Main driver
int StubCleaner::run(TString out, TString src) {
    gROOT->ProcessLine("#include <vector>");  // how is it not loaded?

    int exitcode = 0;
    Timing(1);

    exitcode = readFile(src);
    if (exitcode)  return exitcode;
    Timing();

    exitcode = cleanStubs(out);
    if (exitcode)  return exitcode;
    Timing();

    chain_->Reset();
    return exitcode;
}
