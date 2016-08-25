#ifndef AMSimulation_TrackFitter_h_
#define AMSimulation_TrackFitter_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ProgramOption.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/LocalToGlobalMap.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitterAlgoPCA.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitterAlgoATF.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitterAlgoLTF.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/CombinationFactory.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/CombinationBuilderFactory.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/GhostBuster.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/DuplicateRemoval.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ParameterDuplicateRemoval.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/MCTruthAssociator.h"
using namespace slhcl1tt;


class TrackFitter {

  public:
    // Constructor
    TrackFitter(const ProgramOption& po)
    : po_(po),
      nEvents_(po.maxEvents), verbose_(po.verbose),
    prefixRoad_("AMTTRoads_"), prefixTrack_("AMTTTracks_"), suffix_(""),
    combinationBuilderFactory_(std::make_shared<CombinationBuilderFactory>(po_.FiveOfSix)),
      emptyFirmwareInputString_(66, '0') {

        l2gmap_ = new LocalToGlobalMap();
        l2gmap_->read(po_.datadir);

        // Decide the track fitter to use
        fitter_ = 0;
        if (po.algo == "PCA4" || po.algo == "PCA5") {
            fitter_ = new TrackFitterAlgoPCA(po);
        } else if (po.algo == "ATF4") {
            fitter_ = new TrackFitterAlgoATF(false);
        } else if (po.algo == "ATF5") {
            fitter_ = new TrackFitterAlgoATF(true);
        } else if (po.algo == "LTF") {
            fitter_ = new TrackFitterAlgoLTF(po);
        } else {
            throw std::invalid_argument("unknown track fitter algo.");
        }

        if (po_.emu != 0) {
            firmwareInputFile_.open("FirmwareInputRoads.txt");
            if (!firmwareInputFile_) {
                std::cout << "error opening FirmwareInputRoads.txt" << std::endl;
                return;
            }
            firmwareOutputFile_.open("EmulatorOutput.txt");
            if (!firmwareOutputFile_) {
                std::cout << "error opening FirmwareOutput.txt" << std::endl;
                return;
            }
        }
    }

    // Destructor
    ~TrackFitter() {
        if (l2gmap_)  delete l2gmap_;
        if (fitter_)  delete fitter_;
        firmwareInputFile_.close();
        firmwareOutputFile_.close();
    }

    // Main driver
    int run();

    void writeFirmwareInput(const std::vector<std::vector<unsigned> > & stubRefs, TTRoadReader & reader,
                            const long long & ievt);
    void writeFirmwareOutput(const TTTrack2 & atrack);

  private:
    // Member functions
    int makeTracks(TString src, TString out);

    // Program options
    const ProgramOption po_;
    long long nEvents_;
    int verbose_;

    // Configurations
    const TString prefixRoad_;
    const TString prefixTrack_;
    const TString suffix_;

    // Operators
    LocalToGlobalMap  * l2gmap_;

    // Track fitter
    TrackFitterAlgoBase *  fitter_;

    // Combination factory
    CombinationFactory combinationFactory_;

    // SCB and ACB combination factory
    std::shared_ptr<CombinationBuilderFactory> combinationBuilderFactory_;

    // Ghost buster
    GhostBuster ghostBuster_;

    // MC truth associator
    MCTruthAssociator truthAssociator_;

    std::ofstream firmwareInputFile_;
    std::ofstream firmwareOutputFile_;
    std::string emptyFirmwareInputString_;
};

#endif
