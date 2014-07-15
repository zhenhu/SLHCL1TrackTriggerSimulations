#ifndef AMSimulation_StubCleaner_h_
#define AMSimulation_StubCleaner_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/Helper.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HelperMath.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternBankOption.h"
using namespace slhcl1tt;

#include "TFile.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TTreeFormula.h"
#include "TString.h"


// SETTINGS: none
// INPUT   : TTree with moduleId, hitId, sim info
// OUTPUT  : TTree with moduleId, hitId, sim info

class StubCleaner {
  public:
    // Constructor
    StubCleaner(PatternBankOption option)
    : po(option),
      filter_(true),
      nEvents_(999999999), verbose_(1) {

        chain_ = new TChain("ntupler/tree");

        eventSelect_ = "(1)";  // always on
    }

    // Destructor
    ~StubCleaner() {}


    // Setters
    void setFilter(bool b=true)     { filter_ = b; }
    void setNEvents(int n)          { if (n != -1)  nEvents_ = std::max(0, n); }
    void setVerbosity(int n)        { verbose_ = n; }

    // Getters
    // none

    // Functions
    int readFile(TString src);

    int cleanStubs(TString out);

    // Main driver
    int run(TString out, TString src);


  public:
    // Configurations
    const PatternBankOption po;

  private:
    // Program options
    bool filter_;
    int nEvents_;
    int verbose_;

    // Event selection
    TString eventSelect_;

    // Containers
    TChain * chain_;
};

#endif
