#ifndef ttbarTreeWriter_h
#define ttbarTreeWriter_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
// ROOT stuff
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"

using namespace lcio;
using namespace marlin;


/** 
 *  <h4>Output</h4> 
 *  A ROOT Tree containing information of desired Particles and MCParticles for semi-leptonic ttbar event
 * 
 * @author Philippe Doublet, LAL modified by Jeremy ROUENE
 * @version $Id: ttbarTreeWriter.h,v 1.0 2012/05/30 $ 
 *
 */

class ttbarTreeWriter : public Processor 
{
  
 public:
  
  virtual Processor*  newProcessor() { return new ttbarTreeWriter; }
  
  
  ttbarTreeWriter();

  virtual void init();

  virtual void processRunHeader( LCRunHeader* run );

  virtual void processEvent( LCEvent * evt ); 
    
  virtual void check( LCEvent * evt ); 

  virtual void end();
  
  
 protected:
  

  int _nRun;
  int _nEvt;

enum {
    nJETS=6
};
  TFile* ROOTfile;
  std::string _hfilename;
  TTree* _tree;

  std::string  _clustersCollection;
  std::string  _tracksCollection;
  std::string  _pfosCollection;
  std::string  _mcCollection;
  std::string  _myPfosCollection;
  std::string  _myJetsCollection;
  std::string  _myLeptonCollection;

  TString *_Process;
  float _xsec;

  bool _isttsl;
  bool _iszwwsl;


  /*  Reconstructed variables  */

  TLorentzVector *_jets[nJETS];
  float _jetsBtag[nJETS];
  float _jetsBCtag[nJETS];
  float _jetsCtag[nJETS];
  int _jetsTrueFlavour[nJETS];
  int _jetsTrueCharge[nJETS];
  float _yplus;
  float _yminus;
  bool _isGoodBs;
  bool _isGoodComb;

  TLorentzVector *_myLepton;
  int _nbOfLep;
  int _lepPDG;
  float _lepCharge;
  bool _isGoodLep;
  bool _isLepFromTau;
  TLorentzVector *_totalWOlepton;
  TLorentzVector *_Whad;
  int _selectedBjet;
  TLorentzVector *_top_had;

  int _numberOfPFOs;
  int _numberOfTracks;
  int _numberOfClusters;
  float _pfosFoxWolfram[7];
  float _pfosSphericity;
  float _pfosAplanarity;
  float _pfosOblateness;
  float _pfosThrust;

  /*  Monte Carlo variables   */

  TLorentzVector *_MCtotal;
  
  TLorentzVector *_MCphoton[2];
  TLorentzVector *_MCf[6];
  int _MCfpdg[6];
  
  TLorentzVector *_MCb;
  TLorentzVector *_MCbbar;

  TLorentzVector *_MClepton;
  int _MClepton_PDG;
  TLorentzVector *_MCneutrino;
  int _MCneutrino_PDG;
  TLorentzVector *_MCq[2];
  int _MCqpdg[2];

  TLorentzVector *_MCW_jet;
  TLorentzVector *_MCW_lepton;

  TLorentzVector *_MCtop_jet;
  int _MCtop_jet_Q;

  TLorentzVector *_MCtop_lepton;
  int _MCtop_lepton_Q;
  
  TLorentzVector *_MCbbbar;

};

#endif
