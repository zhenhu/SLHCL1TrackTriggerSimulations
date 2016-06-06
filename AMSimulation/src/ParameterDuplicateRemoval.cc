#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ParameterDuplicateRemoval.h"
using namespace slhcl1tt;

void ParameterDuplicateRemoval::ReduceTracks(std::vector<TTTrack2>& Tracks){
  std::vector<TrackCloud> Clouds; //container for merged tracks
  for(unsigned t=0; t<Tracks.size(); ++t){
    if(Clouds.size()==0){ //start a new cloud, if there isn't any, yet
      TrackCloud element;
      Clouds.push_back(element);
      Clouds[0].Add(t,Tracks.at(t).stubRefs(),Tracks.at(t).chi2()/Tracks.at(t).ndof(),Tracks.at(t).eta(),Tracks.at(t).phi0());
    }
    else{
      bool FoundCloud=false;
      for(unsigned c=0; c<Clouds.size(); ++c){
	FoundCloud=Clouds[c].Add(t,Tracks.at(t).stubRefs(),Tracks.at(t).chi2()/Tracks.at(t).ndof(),Tracks.at(t).eta(),Tracks.at(t).phi0());
	if(FoundCloud) break;
      }//end cloud loop
      if(!FoundCloud){
	TrackCloud element;
	Clouds.push_back(element);
	Clouds[Clouds.size()-1].Add(t,Tracks.at(t).stubRefs(),Tracks.at(t).chi2()/Tracks.at(t).ndof(),Tracks.at(t).eta(),Tracks.at(t).phi0());
      }
    }
  }//end track loop

  //new track collection of unique tracks
  std::vector<TTTrack2> UniqueTracks;
  for(unsigned c=0; c<Clouds.size(); ++c){
    UniqueTracks.push_back(Tracks.at(Clouds[c].BestTrack));
  }//end cloud loop 2
  
  Tracks.resize(UniqueTracks.size());
  Tracks=UniqueTracks;
  assert(Tracks.size() == UniqueTracks.size());
}
