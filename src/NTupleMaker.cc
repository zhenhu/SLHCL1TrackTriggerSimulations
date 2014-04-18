#include "SLHCL1TrackTriggerSimulations/NTupleTools/interface/NTupleMaker.h"

#include "FWCore/Framework/interface/ConstProductRegistry.h"
#include "FWCore/Framework/interface/ProductSelector.h"
#include "FWCore/Framework/interface/ProductSelectorRules.h"
#include "DataFormats/Provenance/interface/Selections.h"
#include "Math/LorentzVector.h"
#include "Math/Vector3D.h"

#include <map>
#include "boost/foreach.hpp"
#include <TBranch.h>
#include <TLorentzVector.h>

void NTupleMaker::
analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  BOOST_FOREACH( BranchConnector* connector, connectors)
    connector->connect(iEvent);
  tree->Fill();
}

template <class T>
void NTupleMaker::TypedBranchConnector<T>::
connect(const edm::Event& iEvent) {
  edm::Handle<T> handle_;
  iEvent.getByLabel(ml, pin, handle_);
  object_ = *handle_;
}

template <class T> 
NTupleMaker::TypedBranchConnector<T>::
TypedBranchConnector(edm::BranchDescription const* desc,
                     std::string t,
                     TTree * tree)
  :  ml( desc->moduleLabel() ),
     pin( desc->productInstanceName() )
{
  object_ptr_ = &object_;
  std::string s=pin+t;
  if(t!="")  { tree->Branch(pin.c_str(),  object_ptr_, s.c_str() );}  //raw type
  else       { tree->Branch(pin.c_str(), &object_ptr_            );}  //vector<type>
}

void NTupleMaker::
beginJob() {
  tree = fs->make<TTree>("tree", "");

  typedef std::map<std::string,       bool> mapStringBool;
  typedef std::map<std::string,        int> mapStringInt;
  typedef std::map<std::string,std::string> mapStringString;
  typedef std::map<std::string,std::vector<float> > mapStringDoubles;
  typedef std::vector<std::vector<float> > vectorVectorFloats;
  typedef std::vector<std::vector<int> > vectorVectorInts;

  std::map<std::string, LEAFTYPE> leafmap;
  leafmap["bool"]      = BOOL;       leafmap["bools"]     = BOOL_V;
  leafmap["short int"] = SHORT;      leafmap["shorts"]    = SHORT_V;
  leafmap["ushort int"]= U_SHORT;    leafmap["ushorts"]   = U_SHORT_V;
  leafmap["int"]       = INT;        leafmap["ints"]      = INT_V;
  leafmap["uint"]      = U_INT;      leafmap["uints"]     = U_INT_V;
  leafmap["float"]     = FLOAT;      leafmap["floats"]    = FLOAT_V;
  leafmap["double"]    = DOUBLE;     leafmap["doubles"]   = DOUBLE_V;
  leafmap["lint"]      = LONG;       leafmap["longs"]     = LONG_V;
  leafmap["ulint"]     = U_LONG;     leafmap["ulongs"]    = U_LONG_V;
  leafmap["StringStringstdmap"] = STRING_STRING_M;
  leafmap["Stringboolstdmap"  ] = STRING_BOOL_M;
  leafmap["Stringintstdmap"   ] = STRING_INT_M;
  leafmap["Stringfloatsstdmap"] = STRING_FLOAT_V_M;
  leafmap["floatss"] = FLOAT_V_V;
  leafmap["intss"] = INT_V_V;
  // leafmap[""] = LORENTZ_V_V;
  
  //
  leafmap["String"]     = STRING;     leafmap["Strings"]    = STRING_V;

  edm::Service<edm::ConstProductRegistry> reg;
  edm::Selections allBranches = reg->allBranchDescriptions();
  edm::ProductSelectorRules groupSelectorRules_(pset, "outputCommands", "NTupleMaker");
  edm::ProductSelector groupSelector_;
  groupSelector_.initialize(groupSelectorRules_, allBranches);

  std::set<std::string> branchnames;

  BOOST_FOREACH( const edm::Selections::value_type& selection, allBranches) {
    if(groupSelector_.selected(*selection)) {

      //Check for duplicate branch names
      if (branchnames.find( selection->productInstanceName()) != branchnames.end() ) {
        throw edm::Exception(edm::errors::Configuration)
          << "More than one branch named: "
          << selection->productInstanceName() << std::endl
          << "Exception thrown from NTupleMaker::beginJob" << std::endl;
      }
      else {
        branchnames.insert( selection->productInstanceName() );
      }

      //Create NTupleMaker branch
      switch(leafmap.find( selection->friendlyClassName() )->second) {
      case BOOL            :  connectors.push_back( new TypedBranchConnector                      <bool>         (selection, "/O", tree) ); break;
      case BOOL_V          :  connectors.push_back( new TypedBranchConnector<std::vector          <bool> >       (selection,   "", tree) ); break;
      case INT             :  connectors.push_back( new TypedBranchConnector                       <int>         (selection, "/I", tree) ); break;
      case INT_V           :  connectors.push_back( new TypedBranchConnector<std::vector           <int> >       (selection,   "", tree) ); break;
      case U_INT           :  connectors.push_back( new TypedBranchConnector              <unsigned int>         (selection, "/i", tree) ); break;
      case U_INT_V         :  connectors.push_back( new TypedBranchConnector<std::vector  <unsigned int> >       (selection,   "", tree) ); break;
      case SHORT           :  connectors.push_back( new TypedBranchConnector                     <short>         (selection, "/S", tree) ); break;
      case SHORT_V         :  connectors.push_back( new TypedBranchConnector<std::vector         <short> >       (selection,   "", tree) ); break;
      case U_SHORT         :  connectors.push_back( new TypedBranchConnector            <unsigned short>         (selection, "/s", tree) ); break;
      case U_SHORT_V       :  connectors.push_back( new TypedBranchConnector<std::vector<unsigned short> >       (selection,   "", tree) ); break;
      case FLOAT           :  connectors.push_back( new TypedBranchConnector                     <float>         (selection, "/F", tree) ); break;
      case FLOAT_V         :  connectors.push_back( new TypedBranchConnector<std::vector         <float> >       (selection,   "", tree) ); break;
      case DOUBLE          :  connectors.push_back( new TypedBranchConnector                    <double>         (selection, "/D", tree) ); break;
      case DOUBLE_V        :  connectors.push_back( new TypedBranchConnector<std::vector        <double> >       (selection,   "", tree) ); break;
      case LONG            :  connectors.push_back( new TypedBranchConnector                      <long>         (selection, "/L", tree) ); break;
      case LONG_V          :  connectors.push_back( new TypedBranchConnector<std::vector          <long> >       (selection,   "", tree) ); break;
      case U_LONG          :  connectors.push_back( new TypedBranchConnector             <unsigned long>         (selection, "/l", tree) ); break;
      case U_LONG_V        :  connectors.push_back( new TypedBranchConnector<std::vector <unsigned long> >       (selection,   "", tree) ); break;
        //	           										         
      case STRING          :  connectors.push_back( new TypedBranchConnector             <std::string  >         (selection,   "", tree) ); break;
      case STRING_V        :  connectors.push_back( new TypedBranchConnector<std::vector <std::string  > >       (selection,   "", tree) ); break;
	
      case STRING_INT_M    :  connectors.push_back( new TypedBranchConnector<mapStringInt>        (selection,   "", tree) ); break;
      case STRING_BOOL_M   :  connectors.push_back( new TypedBranchConnector<mapStringBool>       (selection,   "", tree) ); break;
      case STRING_STRING_M :  connectors.push_back( new TypedBranchConnector<mapStringString>     (selection,   "", tree) ); break;
      case STRING_FLOAT_V_M:  connectors.push_back( new TypedBranchConnector<mapStringDoubles>    (selection,   "", tree) ); break;
      case FLOAT_V_V       :  connectors.push_back( new TypedBranchConnector<vectorVectorFloats>  (selection,   "", tree) ); break;
      case INT_V_V         :  connectors.push_back( new TypedBranchConnector<vectorVectorInts>    (selection,   "", tree) ); break;

      default:
        {
          std::string leafstring = "";
          typedef std::pair<std::string, LEAFTYPE> pair_t;
          BOOST_FOREACH( const pair_t& leaf, leafmap)
            leafstring+= "\t" + leaf.first + "\n";

          throw edm::Exception(edm::errors::Configuration)
            << "class NTupleMaker does not handle leaves of type " << selection->className() << " like\n"
            <<   selection->friendlyClassName()   << "_"
            <<   selection->moduleLabel()         << "_"
            <<   selection->productInstanceName() << "_"
            <<   selection->processName()         << std::endl
            << "Valid leaf types are (friendlyClassName):\n"
            <<   leafstring
            << "Exception thrown from NTupleMaker::beginJob\n";
        }
      }
    }
  }
}

