#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitter.h"

#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTRoadReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTTrackReader.h"
#include <bitset>


namespace {
unsigned getPtSegment(float invPt) {  // for PCA
    return (invPt - PCA_MIN_INVPT) / (PCA_MAX_INVPT - PCA_MIN_INVPT) * PCA_NSEGMENTS;
}

unsigned getHitBits(const std::vector<bool>& stubs_bool) {
    unsigned bitset = 0;

    for (unsigned i=0; i<stubs_bool.size(); ++i) {
        bitset |= (stubs_bool.at(i) << i);
    }

    switch (bitset) {
    case 0b111111:  return 0;
    case 0b111110:  return 1;
    case 0b111101:  return 2;
    case 0b111011:  return 3;
    case 0b110111:  return 4;
    case 0b101111:  return 5;
    case 0b011111:  return 6;
    default      :  return 7;
    }
}

// Comparator
bool sortByPt(const TTTrack2& lhs, const TTTrack2& rhs) {
    return lhs.pt() > rhs.pt();
}
}


// _____________________________________________________________________________
// Do track fitting
int TrackFitter::makeTracks(TString src, TString out) {
    if (verbose_)  std::cout << Info() << "Reading " << nEvents_ << " events and fitting tracks." << std::endl;

    // _________________________________________________________________________
    // For reading
    TTRoadReader reader(verbose_);

    if (reader.init(src, prefixRoad_, suffix_)) {
        std::cout << Error() << "Failed to initialize TTRoadReader." << std::endl;
        return 1;
    }

    // _________________________________________________________________________
    // For writing
    TTTrackWriter writer(verbose_);
    if (writer.init(reader.getChain(), out, prefixTrack_, suffix_)) {
        std::cout << Error() << "Failed to initialize TTTrackWriter." << std::endl;
        return 1;
    }

    // _________________________________________________________________________
    // Loop over all events

    // Containers
    std::vector<TTTrack2> tracks;
    tracks.reserve(300);

    // Bookkeepers
    long int nRead = 0, nKept = 0;

    for (long long ievt=0; ievt<nEvents_; ++ievt) {
        if (reader.loadTree(ievt) < 0)  break;
        reader.getEntry(ievt);

        const unsigned nroads = reader.vr_patternRef->size();
        if (verbose_>1 && ievt%100==0)  std::cout << Debug() << Form("... Processing event: %7lld, fitting: %7ld", ievt, nKept) << std::endl;
        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # roads: " << nroads << std::endl;

        if (!nroads) {  // skip if no road
            writer.fill(std::vector<TTTrack2>());
            ++nRead;
            continue;
        }

        tracks.clear();
        int fitstatus = 0;


        // _____________________________________________________________________
        // Track fitters taking fit combinations

        // Loop over the roads
        for (unsigned iroad=0; iroad<nroads; ++iroad) {
            if (iroad >= (unsigned) po_.maxRoads)  break;

            const unsigned patternRef = reader.vr_patternRef->at(iroad);
            if (patternRef >= (unsigned) po_.maxPatterns)  continue;

            // Get combinations of stubRefs
            std::vector<std::vector<unsigned> > stubRefs = reader.vr_stubRefs->at(iroad);
            //for (unsigned ilayer=0; ilayer<stubRefs.size(); ++ilayer) {
            //    if (stubRefs.at(ilayer).size() > (unsigned) po_.maxStubs)
            //        stubRefs.at(ilayer).resize(po_.maxStubs);
            //}

            if (po_.emu != 0) writeFirmwareInput(stubRefs, reader, ievt);

            // const std::vector<std::vector<unsigned> > & combinations = combinationFactory_.combine(stubRefs);

	    std::vector<std::vector<unsigned> > combinations;
	    // The compiler will likely do RVO so the move might not be needed, we prefer to be explicit about it.
	    if (po_.oldCB) combinations = std::move(combinationFactory_.combine(stubRefs));
	    else combinations = std::move(combinationBuilderFactory_->combine(stubRefs));

            for (unsigned icomb=0; icomb<combinations.size(); ++icomb)
                assert(combinations.at(icomb).size() == reader.vr_stubRefs->at(iroad).size());

            if (verbose_>2) {
                std::cout << Debug() << "... ... road: " << iroad << " # combinations: " << combinations.size() << std::endl;
            }

            // Loop over the combinations
            for (unsigned icomb=0; icomb<combinations.size(); ++icomb) {
                if (icomb >= (unsigned) po_.maxCombs)  break;

                // Create and set TTRoadComb
                TTRoadComb acomb;
                acomb.roadRef    = iroad;
                acomb.combRef    = icomb;
                acomb.patternRef = patternRef;
                acomb.ptSegment  = getPtSegment(reader.vr_patternInvPt->at(iroad));
                acomb.stubRefs   = combinations.at(icomb);

                acomb.stubs_r   .clear();
                acomb.stubs_phi .clear();
                acomb.stubs_z   .clear();
                acomb.stubs_bool.clear();

                float conv_r = 0., conv_phi = 0., conv_z = 0.;
                int64_t conv_r_int = 0, conv_phi_int = 0, conv_z_int = 0;
                LocalToGlobal conv_l2g;
                LocalToGlobalInt conv_l2g_int;

                for (unsigned istub=0; istub<acomb.stubRefs.size(); ++istub) {
                    const unsigned stubRef = acomb.stubRefs.at(istub);
                    if (stubRef != CombinationFactory::BAD) {
                        if (po_.emu == 0) {
                            acomb.stubs_r   .push_back(reader.vb_r   ->at(stubRef));
                            acomb.stubs_phi .push_back(reader.vb_phi ->at(stubRef));
                            acomb.stubs_z   .push_back(reader.vb_z   ->at(stubRef));
                            acomb.stubs_bool.push_back(true);
                            acomb.stubs_bitString.push_back("");
                        } else {
                            acomb.stubs_r   .push_back(reader.vb_r   ->at(stubRef));  // full floating-point precision, not fixed-point precision
                            acomb.stubs_phi .push_back(reader.vb_phi ->at(stubRef));  // full floating-point precision, not fixed-point precision
                            acomb.stubs_z   .push_back(reader.vb_z   ->at(stubRef));  // full floating-point precision, not fixed-point precision
                            acomb.stubs_bool.push_back(true);
                            acomb.stubs_bitString.push_back(reader.vb_bitString->at(stubRef));

                            // This block is safe to remove
                            {
                                // Do local-to-global conversion
                                unsigned moduleId = reader.vb_modId   ->at(stubRef);
                                float    strip    = reader.vb_coordx  ->at(stubRef);  // in full-strip unit
                                float    segment  = reader.vb_coordy  ->at(stubRef);  // in full-strip unit
                                float    stub_ds  = reader.vb_trigBend->at(stubRef);  // in full-strip unit
                                l2gmap_ -> convert(moduleId, strip, segment, conv_r, conv_phi, conv_z, conv_l2g);
                                l2gmap_ -> convertInt(moduleId, strip, segment, po_.tower, conv_l2g, conv_r_int, conv_phi_int, conv_z_int, conv_l2g_int);

                                int bend_4b = int(std::round(stub_ds));
                                int strip_7b = (halfStripRound(strip)) >> 4;

                                // Sanity checks
                                if (!acomb.stubs_bitString.empty()) {
                                    size_t index = acomb.stubs_bitString.size() - 1;
                                    assert(acomb.stubs_phi_int(index) == conv_phi_int);
                                    assert(acomb.stubs_r_int(index) == conv_r_int);
                                    assert(acomb.stubs_z_int(index) == conv_z_int);
                                    assert(acomb.bend_int(index) == bend_4b);
                                    assert(acomb.strip_int(index) == strip_7b);
                                }
                                const unsigned superstripId = reader.vr_superstripIds->at(iroad).at(istub);
                                const unsigned superstripId2 = reader.vb_superstripId->at(stubRef);
                                assert(superstripId2 == superstripId);
                            }
                        }
                    } else {
                        acomb.stubs_r   .push_back(0.);
                        acomb.stubs_phi .push_back(0.);
                        acomb.stubs_z   .push_back(0.);
                        acomb.stubs_bool.push_back(false);
                        acomb.stubs_bitString.push_back(emptyFirmwareInputString_);
                    }
                }

                acomb.hitBits = getHitBits(acomb.stubs_bool);

                if (verbose_>2) {
                    std::cout << Debug() << "... ... ... comb: " << icomb << " " << acomb;
                    std::cout << std::endl;
                }

                // _____________________________________________________________
                // Fit
                TTTrack2 atrack;
                fitstatus = fitter_->fit(acomb, atrack);

                if (po_.emu != 0) writeFirmwareOutput(atrack);

                atrack.setTower     (po_.tower);
                atrack.setRoadRef   (acomb.roadRef);
                atrack.setCombRef   (acomb.combRef);
                atrack.setPatternRef(acomb.patternRef);
                atrack.setPtSegment (acomb.ptSegment);
                atrack.setHitBits   (acomb.hitBits);
                atrack.setStubRefs  (acomb.stubRefs);

                if (atrack.chi2Red() < po_.maxChi2)  // reduced chi^2 = chi^2 / ndof
                    tracks.push_back(atrack);

                if (verbose_>2)  std::cout << Debug() << "... ... ... track: " << icomb << " status: " << fitstatus << " reduced chi2: " << atrack.chi2Red() << " invPt: " << atrack.invPt() << " phi0: " << atrack.phi0() << " cottheta: " << atrack.cottheta() << " z0: " << atrack.z0() << std::endl;
            }
        }  // loop over the roads

        std::sort(tracks.begin(), tracks.end(), sortByPt);


        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # tracks: " << tracks.size() << std::endl;
        if (verbose_>3) {
            for (unsigned itrack=0; itrack!=tracks.size(); ++itrack) {
                std::cout << "... ... track: " << itrack << " " << tracks.at(itrack) << std::endl;
            }
        }

        // _____________________________________________________________________
        // Find ghosts

        for (unsigned itrack=0; itrack<tracks.size(); ++itrack) {  // all tracks
            for (unsigned jtrack=0; jtrack<itrack; ++jtrack) {  // only non ghost tracks
                if (tracks.at(jtrack).isGhost())  continue;

                bool isGhost = ghostBuster_.isGhostTrack(tracks.at(jtrack).stubRefs(), tracks.at(itrack).stubRefs());
                if (isGhost) {
                    tracks.at(itrack).setAsGhost();
                }
            }
        }

        if (! tracks.empty())
            ++nKept;

        if (tracks.size() > (unsigned) po_.maxTracks)
            tracks.resize(po_.maxTracks);



	// ---------------------------------------------------------------------
	// Classify tracks as duplicates or not (for duplicate removal)
	// In the algorithm tracking particles are sorted by pT
	// And AM tracks are sorted by logic and pT
	// ---------------------------------------------------------------------
	DuplicateRemoval flagDuplicates;
	flagDuplicates.CheckTracks(tracks, po_.rmDuplicate);

	//----------------------------------------------------------------------
	// Identify and flag duplicates by defining a track-parameter space 
	// inside of which anything is considered to be a single track
	//----------------------------------------------------------------------
	ParameterDuplicateRemoval RemoveParameterDuplicates;
	if(po_.rmParDuplicate) RemoveParameterDuplicates.ReduceTracks(tracks);




        // _____________________________________________________________________        // Track categorization

        if (po_.speedup<1) {
            const unsigned nparts = reader.vp2_primary->size();
            if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # particles: " << nparts << std::endl;

            std::vector<TrackingParticle> trkParts;
            for (unsigned ipart=0; ipart<nparts; ++ipart) {
                bool  primary         = reader.vp2_primary->at(ipart);
                int   simCharge       = reader.vp2_charge->at(ipart);

                if (simCharge!=0 && primary) {
                    float simPt           = reader.vp2_pt->at(ipart);
                    float simEta          = reader.vp2_eta->at(ipart);
                    float simPhi          = reader.vp2_phi->at(ipart);
                    //float simVx           = reader.vp2_vx->at(ipart);
                    //float simVy           = reader.vp2_vy->at(ipart);
                    float simVz           = reader.vp2_vz->at(ipart);
                    int   simCharge       = reader.vp2_charge->at(ipart);
                    int   simPdgId        = reader.vp2_pdgId->at(ipart);

                    float simCotTheta     = std::sinh(simEta);
                    float simChargeOverPt = float(simCharge)/simPt;

                    trkParts.emplace_back(TrackingParticle{  // using POD type constructor
                        (int) ipart,
                        simPdgId,
                        simChargeOverPt,
                        simPhi,
                        simCotTheta,
                        simVz,
                        0.
                    });

                    if (verbose_>3)  std::cout << Debug() << "... ... part: " << ipart << " primary: " << primary << " " << trkParts.back();
                }
            }
            truthAssociator_.associate(trkParts, tracks);
        }

        writer.fill(tracks);
        ++nRead;
    }

    if (nRead == 0) {
        std::cout << Error() << "Failed to read any event." << std::endl;
        return 1;
    }

    if (verbose_)  std::cout << Info() << Form("Read: %7ld, triggered: %7ld", nRead, nKept) << std::endl;


    // _________________________________________________________________________
    // Write histograms

    for (std::map<TString, TH1F *>::const_iterator it=fitter_->histograms_.begin();
         it!=fitter_->histograms_.end(); ++it) {
        if (it->second)  it->second->SetDirectory(gDirectory);
    }

    long long nentries = writer.writeTree();
    assert(nentries == nRead);

    return 0;
}


// _____________________________________________________________________________
// Main driver
int TrackFitter::run() {
    int exitcode = 0;
    Timing(1);

    exitcode = makeTracks(po_.input, po_.output);
    if (exitcode)  return exitcode;
    Timing();

    return exitcode;
}


void TrackFitter::writeFirmwareInput(const std::vector<std::vector<unsigned> > & stubRefs, TTRoadReader & reader,
                                     const long long & ievt)
{
    // Keep track of the event id
    firmwareInputFile_ << std::setfill('0') << std::setw(5) << ievt << " ";
    unsigned ilayer = 0;
    for (auto layer : stubRefs) {
        // std::cout << "layer " << ilayer << ": ";
        unsigned numStubs = 0;
        for (auto stubRef = layer.rbegin(); stubRef != layer.rend(); ++stubRef) {
            // std::cout << reader.vb_bitString->at(*stubRef) << " ";
            firmwareInputFile_ << reader.vb_bitString->at(*stubRef);
            ++numStubs;
        }
        // Fill up missing stubs with zeros
        for (unsigned missingStubs = numStubs; missingStubs < 4; ++missingStubs) {
            // std::cout << emptyFirmwareInputString_ << " ";
            firmwareInputFile_ << emptyFirmwareInputString_;
        }
        ++ilayer;
        // std::cout << std::endl;
    }
    firmwareInputFile_ << std::endl;
}


void TrackFitter::writeFirmwareOutput(const TTTrack2 & atrack)
{
    std::vector<int64_t> parsInt(atrack.parsInt());
    std::vector<int64_t> chi2TermsInt(atrack.chi2TermsInt());
    if (parsInt.size() == 4) {
        // std::cout << parsInt.at(0) << " ";
        // std::cout << parsInt.at(1) << " ";
        // std::cout << chi2TermsInt.at(0) << " ";
        // std::cout << chi2TermsInt.at(1) << " ";
        // std::cout << chi2TermsInt.at(2) << " ";
        // std::cout << chi2TermsInt.at(3) << " ";
        // std::cout << parsInt.at(2) << " ";
        // std::cout << parsInt.at(3) << " ";
        // std::cout << chi2TermsInt.at(4) << " ";
        // std::cout << chi2TermsInt.at(5) << " ";
        // std::cout << chi2TermsInt.at(6) << " ";
        // std::cout << chi2TermsInt.at(7) << " ";
        // std::cout << std::endl;
        firmwareOutputFile_ << parsInt.at(0) << " ";
        firmwareOutputFile_ << parsInt.at(1) << " ";
        firmwareOutputFile_ << chi2TermsInt.at(0) << " ";
        firmwareOutputFile_ << chi2TermsInt.at(1) << " ";
        firmwareOutputFile_ << chi2TermsInt.at(2) << " ";
        firmwareOutputFile_ << chi2TermsInt.at(3) << " ";
        firmwareOutputFile_ << parsInt.at(2) << " ";
        firmwareOutputFile_ << parsInt.at(3) << " ";
        firmwareOutputFile_ << chi2TermsInt.at(4) << " ";
        firmwareOutputFile_ << chi2TermsInt.at(5) << " ";
        firmwareOutputFile_ << chi2TermsInt.at(6) << " ";
        firmwareOutputFile_ << chi2TermsInt.at(7) << " ";
        firmwareOutputFile_ << std::endl;
    }
}
