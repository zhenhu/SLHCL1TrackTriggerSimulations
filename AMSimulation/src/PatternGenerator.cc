#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternGenerator.h"

static const unsigned MIN_NGOODSTUBS = 3;
static const unsigned MAX_NTOWERS = 6 * 8;


// _____________________________________________________________________________
// Read and parse the trigger tower csv file
int PatternGenerator::readTriggerTowerFile(TString src) {
    std::vector<unsigned> values;
    if (src.EndsWith(".csv")) {
        if (verbose_)  std::cout << Info() << "Opening " << src << std::endl;

        unsigned i = 0;
        std::string line;
        std::ifstream ifs(src.Data());
        while (std::getline(ifs, line)) {
            std::istringstream iss(line);
            std::string issline;

            while (std::getline(iss, issline, ',')) {  // split by commas
                if (i != 0)  // skip the first line
                    values.push_back(std::stoi(issline));
            }
            ++i;
        }
    }

    if (values.empty()) {
        std::cout << Error() << "Failed to read the trigger tower file: " << src << std::endl;
        return 1;
    }

    triggerTowerMap_.clear();
    unsigned triggerTowerId = 0;
    for (unsigned i=0; i<values.size(); ++i) {
        if (i == 0)  continue;  // skip the first index
        if (values.at(i-1) <= 6 && values.at(i) <= 8) {  // eta_idx, phi_idx
            triggerTowerId = (values.at(i-1)-1) * 8 + (values.at(i)-1);  // mapped to 0-47
            triggerTowerMap_.insert(std::make_pair(triggerTowerId, std::vector<unsigned>()));
        } else if (values.at(i) > 10000) {
            triggerTowerMap_.at(triggerTowerId).push_back(values.at(i));
        }
    }

    if (triggerTowerMap_.size() != MAX_NTOWERS) {
        std::cout << Error() << "Failed to get all trigger towers from the file: " << src << std::endl;
        return 1;
    }

    //for (auto it: triggerTowerMap_)
    //    std::cout << "Tower " << it.first << " has " << it.second.size() << " modules." << std::endl;

    // Prepare the trigger tower reverse map
    triggerTowerReverseMap_.clear();
    //for (auto it: triggerTowerMap_) { // loop over trigger tower map
    //    for (auto it2: it.second) {   // loop over the moduleIds in the tower
    //        triggerTowerReverseMap_[it2].push_back(it.first);
    //    }
    //}
    for (auto it: po.triggerTowers) { // loop over input trigger towers
        const std::vector<unsigned>& moduleIds = triggerTowerMap_[it];
        for (auto it2: moduleIds) {   // loop over the moduleIds in the tower
            triggerTowerReverseMap_[it2].push_back(it);
        }
    }

    //for (auto it: triggerTowerReverseMap_)
    //    std::cout << "Module " << it.first << " is in " << it.second.size() << " towers." << std::endl;

    return 0;
}


// _____________________________________________________________________________
// Read the input ntuples
int PatternGenerator::readFile(TString src) {
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
// Make the patterns
int PatternGenerator::makePatterns_map() {
    long long nentries = chain_->GetEntries();
    if (nentries <= 0) {
        std::cout << Error() << "Input source has zero entry." << std::endl;
        return 1;
    }
    if (nEvents_ > nentries)
        nEvents_ = nentries;

    if (verbose_)  std::cout << Info() << "Reading " << nEvents_ << " events" << std::endl;

    // For reading
    std::vector<float> *          vb_x          = 0;
    std::vector<float> *          vb_y          = 0;
    std::vector<float> *          vb_z          = 0;
    std::vector<float> *          vb_r          = 0;
    std::vector<float> *          vb_phi        = 0;
    std::vector<float> *          vb_coordx     = 0;
    std::vector<float> *          vb_coordy     = 0;
    std::vector<float> *          vb_roughPt    = 0;
    std::vector<unsigned> *       vb_modId      = 0;
    std::vector<unsigned> *       vb_nhits      = 0;
    std::vector<float> *          vb_simPt      = 0;
    std::vector<float> *          vb_simEta     = 0;
    std::vector<float> *          vb_simPhi     = 0;
    std::vector<int> *            vb_trkId      = 0;

    chain_->SetBranchStatus("*"                 , 0);
    chain_->SetBranchStatus("TTStubs_x"         , 1);
    chain_->SetBranchStatus("TTStubs_y"         , 1);
    chain_->SetBranchStatus("TTStubs_z"         , 1);
    chain_->SetBranchStatus("TTStubs_r"         , 1);
    chain_->SetBranchStatus("TTStubs_phi"       , 1);
    chain_->SetBranchStatus("TTStubs_coordx"    , 1);
    chain_->SetBranchStatus("TTStubs_coordy"    , 1);
    chain_->SetBranchStatus("TTStubs_roughPt"   , 1);
    chain_->SetBranchStatus("TTStubs_modId"     , 1);
    chain_->SetBranchStatus("TTStubs_nhits"     , 1);
    chain_->SetBranchStatus("TTStubs_simPt"     , 1);
    chain_->SetBranchStatus("TTStubs_simEta"    , 1);
    chain_->SetBranchStatus("TTStubs_simPhi"    , 1);
    chain_->SetBranchStatus("TTStubs_trkId"     , 1);

    chain_->SetBranchAddress("TTStubs_x"        , &(vb_x));
    chain_->SetBranchAddress("TTStubs_y"        , &(vb_y));
    chain_->SetBranchAddress("TTStubs_z"        , &(vb_z));
    chain_->SetBranchAddress("TTStubs_r"        , &(vb_r));
    chain_->SetBranchAddress("TTStubs_phi"      , &(vb_phi));
    chain_->SetBranchAddress("TTStubs_coordx"   , &(vb_coordx));
    chain_->SetBranchAddress("TTStubs_coordy"   , &(vb_coordy));
    chain_->SetBranchAddress("TTStubs_roughPt"  , &(vb_roughPt));
    chain_->SetBranchAddress("TTStubs_modId"    , &(vb_modId));
    chain_->SetBranchAddress("TTStubs_nhits"    , &(vb_nhits));
    chain_->SetBranchAddress("TTStubs_simPt"    , &(vb_simPt));
    chain_->SetBranchAddress("TTStubs_simEta"   , &(vb_simEta));
    chain_->SetBranchAddress("TTStubs_simPhi"   , &(vb_simPhi));
    chain_->SetBranchAddress("TTStubs_trkId"    , &(vb_trkId));

    allPatterns_map_.clear();
    if ((size_t) nEvents_ >= allPatterns_map_.max_size()) {
        std::cout << Error() << "Number of events is more than the max_size of a std::map." << std::endl;
        return 1;
    }

    // _________________________________________________________________________
    // Loop over all events

    // Containers are declared outside the event loop to reduce memory allocations
    std::vector<id_type> superstripLayers;
    std::vector<addr_type> superstrips;
    std::vector<pattern_type> patterns;

    int nRead = 0, nKept = 0;
    for (long long ievt=0; ievt<nEvents_; ++ievt) {
        Long64_t local_entry = chain_->LoadTree(ievt);  // for TChain
        if (local_entry < 0)  break;
        chain_->GetEntry(ievt);

        unsigned nstubs = vb_modId->size();
        if (verbose_>1 && ievt%1000==0)  std::cout << Debug() << Form("... Processing event: %7lld, keeping: %7i", ievt, nKept) << std::endl;
        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # stubs: " << nstubs << " [# patterns: " << allPatterns_map_.size() << "]" << std::endl;

        if (!nstubs) {  // skip if no stub
            continue;
        }

        // _____________________________________________________________________
        // Start generating patterns
        bool keep = true;

        superstripLayers.clear();
        superstrips.clear();
        patterns.clear();

        // Check min # of layers
        bool require = (nstubs >= MIN_NGOODSTUBS);
        if (!require)
            keep = false;

        // Loop over reconstructed stubs
        for (unsigned l=0; (l<nstubs) && keep; ++l) {
            unsigned moduleId = vb_modId->at(l);
            // If there is a valid hit, but moduleId does not exist in any
            // trigger tower (due to cables, connections, or just too forward),
            // we drop them
            if (po.requireTriggerTower && triggerTowerReverseMap_.find(moduleId) == triggerTowerReverseMap_.end())
                continue;

            unsigned lay = decodeLayer(moduleId);
            unsigned count = std::count(superstripLayers.begin(), superstripLayers.end(), lay);
            if (count != 0) {
                std::cout << Error() << "There should be only one stub in any layer" << std::endl;
                return 1;
            }
            superstripLayers.push_back(lay);

            unsigned nhits = vb_nhits->at(l);
            assert(nhits > 0);
            float coordx = vb_coordx->at(l);
            float coordy = vb_coordy->at(l);

            // Use half-strip unit
            id_type col = halfStripRound(coordy);
            id_type row = halfStripRound(coordx);

            // Find superstrip address
            col = arbiter_ -> subladder(moduleId, col);
            row = arbiter_ -> submodule(moduleId, row);
            addr_type ssId = encodeSuperstripId(moduleId, col, row);

            superstrips.push_back(ssId);

            if (verbose_>2)  std::cout << Debug() << "... ... stub: " << l << " moduleId: " << moduleId << " col: " << col << " row: " << row << " ssId: " << ssId << " trkId: " << vb_trkId->at(l) << " # hits: " << nhits << std::endl;
        }

        // Build the patterns
        if (!superstrips.empty()) {
            std::sort(superstrips.begin(), superstrips.end(), std::less<addr_type>());
            patterns = stitcher_ -> stitch(superstrips);

            // Remove patterns that are not within any trigger tower
            if (po.requireTriggerTower) {
                pattern_type emptyPattern;
                for (unsigned i=0; i<patterns.size(); ++i) {
                    //assert(patterns.at(i) != emptyPattern);
                    if (!isWithinTriggerTower(patterns.at(i)) )
                        patterns.at(i).fill(0);
                }
                patterns.erase(std::remove(patterns.begin(), patterns.end(), emptyPattern));
            }
        }

        if (patterns.empty())
            keep = false;

        if (keep) {
            if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # patterns: " << patterns.size() << std::endl;

            for (unsigned i=0; i<patterns.size(); ++i) {
                ++allPatterns_map_[patterns.at(i)];
                if (verbose_>2)  std::cout << Debug() << "... ... patt: " << i << "  " << patterns.at(i) << std::endl;
            }
            ++nKept;
        }
        ++nRead;
    }
    if (verbose_)  std::cout << Info() << Form("Read: %7i, kept: %7i, # patterns: %7lu", nRead, nKept, allPatterns_map_.size()) << std::endl;

    for (auto it: allPatterns_map_)
        allPatterns_map_pairs_.push_back(it);

    // Clear the map and release memory
    std::map<pattern_type, count_type> emptyMap;
    allPatterns_map_.clear();
    allPatterns_map_.swap(emptyMap);

    // Sort by frequency
    std::sort(allPatterns_map_pairs_.begin(), allPatterns_map_pairs_.end(),
              [=](const std::pair<pattern_type, count_type>& lhs, const std::pair<pattern_type, count_type>& rhs) {
                  return lhs.second > rhs.second;
              } );

    if (verbose_>2) {
        unsigned i=0;
        for (auto it : allPatterns_map_pairs_) {
            std::cout << Debug() << "... patt: " << i << "  " << it.first << "  freq: " << (unsigned) it.second << std::endl;
            ++i;
        }
    }

    if (verbose_)  std::cout << Info() << "Generated " << allPatterns_map_pairs_.size() << " patterns, highest freq: " << (unsigned) allPatterns_map_pairs_.front().second << std::endl;

    return 0;
}


// _____________________________________________________________________________
// Output patterns into a TTree
int PatternGenerator::writePatterns_map(TString out) {
    if (!out.EndsWith(".root")) {
        std::cout << Error() << "Output filename must be .root" << std::endl;
        return 1;
    }

    const long long nentries = allPatterns_map_pairs_.size();
    if (verbose_)  std::cout << Info() << "Recreating " << out << " with " << nentries << " patterns." << std::endl;
    TFile* tfile = TFile::Open(out, "RECREATE");
    TTree* ttree = new TTree(bankName_, "");

    typedef unsigned char unsigned8;
    std::auto_ptr<unsigned8>                  frequency       (new unsigned8(0));
    std::auto_ptr<std::vector<unsigned> >     superstripIds   (new std::vector<unsigned>());

    ttree->Branch("frequency"       , &(*frequency));
    ttree->Branch("superstripIds"   , &(*superstripIds));

    // _________________________________________________________________________
    // Loop over all patterns
    for (long long ievt=0; ievt<nentries; ++ievt) {
        if (verbose_>1 && ievt%10000==0)  std::cout << Debug() << Form("... Writing event: %7lld", ievt) << std::endl;
        superstripIds->clear();

        *frequency = allPatterns_map_pairs_.at(ievt).second;

        // Assume patterns are sorted by frequency
        if (*frequency < minFrequency_)
            break;

        const pattern_type& patt = allPatterns_map_pairs_.at(ievt).first;
        for (unsigned i=0; i<patt.size(); ++i) {
            superstripIds->push_back(patt.at(i));
        }

        ttree->Fill();
    }
    assert(ttree->GetEntries() == nentries);

    tfile->Write();
    delete ttree;
    delete tfile;
    return 0;
}


// _____________________________________________________________________________
// Private functions

bool PatternGenerator::isWithinTriggerTower(const pattern_type& patt) {
    // FIXME: implement this
    return true;
}


// _____________________________________________________________________________
// Main driver
int PatternGenerator::run(TString out, TString src, TString layout) {
    gROOT->ProcessLine("#include <vector>");  // how is it not loaded?

    int exitcode = 0;
    Timing(1);

    if (po.requireTriggerTower) {
        exitcode = readTriggerTowerFile(layout);
        if (exitcode)  return exitcode;
        Timing();
    }

    exitcode = readFile(src);
    if (exitcode)  return exitcode;
    Timing();

    exitcode = makePatterns_map();
    if (exitcode)  return exitcode;
    Timing();

    exitcode = writePatterns_map(out);
    if (exitcode)  return exitcode;
    Timing();

    chain_->Reset();
    return exitcode;
}
