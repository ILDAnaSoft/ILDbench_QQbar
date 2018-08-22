#include "ttbarTreeWriter.h"
#include <iostream>
#include <math.h>
#include <sstream>

#include <LCIOSTLTypes.h>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/PIDHandler.h>
#include <EVENT/MCParticle.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

ttbarTreeWriter attbarTreeWriter;

ttbarTreeWriter::ttbarTreeWriter() : Processor("ttbarTreeWriter")
{
  
  // processor description

  _description = "Analyse ttbar to semileptonic mode at Reconstruction and MC level";
  
  
  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "ROOTFileName",
                              "Output ROOT File Name",
                              _hfilename,
                              std::string("ttbarTreeWriter.root") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "PandoraCollection" ,
			   "Output collection given by Pandora" ,
			   _pfosCollection ,
			   std::string("PandoraPFOs") );

  registerInputCollection( LCIO::MCPARTICLE,
			   "MCCollection" ,
			   "MC collection (MCParticle for Mokka, MCParticlesSkimmed for DST)" ,
			   _mcCollection ,
			   std::string("MCParticle") );
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "MyPFOsCollection" ,
                           "Name of my collection"  ,
                           _myPfosCollection ,
                           std::string("PhD_PFOs") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "4Jets" ,
                           "Name of the 4 Jets collection"  ,
                           _myJetsCollection ,
                           std::string("PhD_Durham4Jets") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "Lepton" ,
                           "Name of the lepton collection"  ,
                           _myLeptonCollection ,
                           std::string("PhD_SelectedLepton") );

}


void ttbarTreeWriter::init()
{ 

  std::cout<< " ---------  ttbar Tree Writer Processor Start ------  " << std::endl;  
  _nRun = 0;
  _nEvt = 0;

  ROOTfile = new TFile(_hfilename.c_str(), "RECREATE", _hfilename.c_str());
  
  _tree = new TTree("tree", "tree");
  
  _tree->Branch("runNumber", &_nRun, "runNumber/I");
  _tree->Branch("eventNumber", &_nEvt, "eventNumber/I");
  _tree->Branch("crossSection", &_xsec, "crossSection/F");

  _Process = new TString();
  _tree->Branch("Process","TString",&_Process,16000,0);

  _tree->Branch("isTTbarSL", &_isttsl, "isTTbarSL/O");
  _tree->Branch("isZWWSL", &_iszwwsl, "isZWWSL/O");


  //  Reconstructed variables 

  //  for(int i = 0 ; i < 4 ; i++)
  for(int i = 0 ; i < nJETS ; i++)
    {
      _jets[i] = new TLorentzVector();
      std::stringstream name;
      name << "jet" << i;
      _tree->Branch(name.str().c_str(),"TLorentzVector",&_jets[i],16000,0);
      _tree->Branch((name.str()+"_btag").c_str(), &_jetsBtag[i], (name.str()+"_btag/F").c_str());
      _tree->Branch((name.str()+"_bctag").c_str(), &_jetsBCtag[i], (name.str()+"_bctag/F").c_str());
      _tree->Branch((name.str()+"_ctag").c_str(), &_jetsCtag[i], (name.str()+"_ctag/F").c_str());
      _tree->Branch((name.str()+"_TrueFlavour").c_str(), &_jetsTrueFlavour[i], (name.str()+"_TrueFlavour/I").c_str());
      _tree->Branch((name.str()+"_TrueCharge").c_str(), &_jetsTrueCharge[i], (name.str()+"_TrueCharge/I").c_str());
    }

  _tree->Branch("YPlus", &_yplus, "YPlus/F");
  _tree->Branch("YMinus", &_yminus, "YMinus/F");
  _tree->Branch("isGoodBjets", &_isGoodBs, "isGoodBjets/O");
  _tree->Branch("isGoodCombination", &_isGoodComb, "isGoodCombination/O");

  _myLepton = new TLorentzVector();
  _tree->Branch("lepton","TLorentzVector",&_myLepton,16000,0);
  _tree->Branch("leptonPDG", &_lepPDG, "leptonPDG/I");
  _tree->Branch("leptonCharge", &_lepCharge, "leptonCharge/F");
  _tree->Branch("NbCandidates", &_nbOfLep, "NbCandidates/I");
  _tree->Branch("isGoodLepton", &_isGoodLep, "isGoodLepton/O");
  _tree->Branch("isLepFromTau", &_isLepFromTau, "isLepFromTau/O");
  _totalWOlepton = new TLorentzVector();
  _tree->Branch("allWOlepton","TLorentzVector",&_totalWOlepton,16000,0);
  _Whad = new TLorentzVector();
  _tree->Branch("Whad","TLorentzVector",&_Whad,16000,0);
  _tree->Branch("SelectedBjet", &_selectedBjet, "SelectedBjet/I");
  _top_had = new TLorentzVector();
  _tree->Branch("top_had","TLorentzVector",&_top_had,16000,0);

  
  _tree->Branch("numberOfPFOs", &_numberOfPFOs, "numberOfPFOs/I");
  _tree->Branch("pfosFoxWolfram", _pfosFoxWolfram, "pfosFoxWolfram[7]/F");
  _tree->Branch("pfosSphericity", &_pfosSphericity, "pfosSphericity/F");
  _tree->Branch("pfosAplanarity", &_pfosAplanarity, "pfosAplanarity/F");
  _tree->Branch("pfosOblateness", &_pfosOblateness, "pfosOblateness/F");
  _tree->Branch("pfosThrust", &_pfosThrust, "pfosThrust/F");


  // Monte Carlo variables 

  _MCtotal = new TLorentzVector();
  _tree->Branch("MCtotal","TLorentzVector",&_MCtotal,16000,0);
  
  _MCphoton[0] = new TLorentzVector();
  _tree->Branch("MCphoton1","TLorentzVector",&_MCphoton[0],16000,0);
  _MCphoton[1] = new TLorentzVector();
  _tree->Branch("MCphoton2","TLorentzVector",&_MCphoton[1],16000,0);
  
  
  for(int i = 0 ; i < nJETS ; i++)
    {
      _MCf[i] = new TLorentzVector();
      std::stringstream name;
      name << "MCfermion" << i;
      _tree->Branch(name.str().c_str(),"TLorentzVector",&_MCf[i],16000,0);
      name << "_PDG";
      _tree->Branch(name.str().c_str(), &_MCfpdg[i], (name.str()+"/I").c_str());
    }
  
  _MCb = new TLorentzVector();
  _tree->Branch("MCb","TLorentzVector",&_MCb,16000,0);

  _MCbbar = new TLorentzVector();
  _tree->Branch("MCbbar","TLorentzVector",&_MCbbar,16000,0);
  
  _MClepton = new TLorentzVector();
  _tree->Branch("MClepton","TLorentzVector",&_MClepton,16000,0);
  _tree->Branch("MClepton_PDG", &_MClepton_PDG, "MClepton_PDG/I");
  
  _MCneutrino = new TLorentzVector();
  _tree->Branch("MCneutrino","TLorentzVector",&_MCneutrino,16000,0);
  _tree->Branch("MCneutrino_PDG", &_MCneutrino_PDG, "MCneutrino_PDG/I");

  _MCq[0] = new TLorentzVector();
  _tree->Branch("MCq1","TLorentzVector",&_MCq[0],16000,0);
  _tree->Branch("MCq1_PDG", &_MCqpdg[0], "MCq1_PDG/I");  

  _MCq[1] = new TLorentzVector();
  _tree->Branch("MCq2","TLorentzVector",&_MCq[1],16000,0);
  _tree->Branch("MCq2_PDG", &_MCqpdg[1], "MCq2_PDG/I");  
  
  _MCW_jet = new TLorentzVector();
  _tree->Branch("MCW_jet","TLorentzVector",&_MCW_jet,16000,0);

  _MCW_lepton = new TLorentzVector();
  _tree->Branch("MCW_lepton","TLorentzVector",&_MCW_lepton,16000,0);

  _MCtop_jet = new TLorentzVector();
  _tree->Branch("MCtop_jet","TLorentzVector",&_MCtop_jet,16000,0);
  _tree->Branch("MCtop_jet_sign", &_MCtop_jet_Q, "MCtop_jet_sign/I");

  _MCtop_lepton = new TLorentzVector();
  _tree->Branch("MCtop_lepton","TLorentzVector",&_MCtop_lepton,16000,0);
  _tree->Branch("MCtop_lepton_sign", &_MCtop_lepton_Q, "MCtop_lepton_sign/I");
 
  _MCbbbar = new TLorentzVector();
  _tree->Branch("MCbbbar","TLorentzVector",&_MCbbbar,16000,0);
  
  
}


void ttbarTreeWriter::processRunHeader(LCRunHeader* run)
{ 
  _nRun++;
} 


void ttbarTreeWriter::processEvent(LCEvent * evt) 
{ 
  
  _nEvt ++;
  ROOTfile->cd();
  
  try {

    // Get Process name and cross section
    *_Process = evt->getParameters().getStringVal("Process");
    _xsec = evt->getParameters().getFloatVal("CrossSection_fb");
    
    _isttsl = false;
    _iszwwsl = false;
 
    // Get information on the MCParticles
    LCCollection* mccol = evt->getCollection(_mcCollection);
    std:: cout << " MC Collection contains " << mccol->getNumberOfElements() << " Elements. " << std::endl;
    if(mccol->getNumberOfElements() > 8)
      {
	// Initialise
	_MCtotal->SetPxPyPzE(0.,0.,0.,0.);
	_MCtop_jet->SetPxPyPzE(0.,0.,0.,0.);
	_MCtop_jet_Q = 0;
	_MCtop_lepton->SetPxPyPzE(0.,0.,0.,0.);
	_MCtop_lepton_Q = 0;
	_MCphoton[0]->SetPxPyPzE(0.,0.,0.,0.);
	_MCphoton[1]->SetPxPyPzE(0.,0.,0.,0.);
	_MCq[0]->SetPxPyPzE(0.,0.,0.,0.) ; _MCqpdg[0] = 0;
	_MCq[1]->SetPxPyPzE(0.,0.,0.,0.) ; _MCqpdg[1] = 0;
	for(int i = 0 ; i < 6 ; i++)
	  { 
	    _MCf[i]->SetPxPyPzE(0.,0.,0.,0.);
	    _MCfpdg[i] = 0; 
	  }
	_MCW_jet->SetPxPyPzE(0.,0.,0.,0.);
	_MCW_lepton->SetPxPyPzE(0.,0.,0.,0.);
	_MCb->SetPxPyPzE(0.,0.,0.,0.);
	_MCbbar->SetPxPyPzE(0.,0.,0.,0.);
	_MClepton->SetPxPyPzE(0.,0.,0.,0.);
	_MClepton_PDG = 0;
	_MCneutrino->SetPxPyPzE(0.,0.,0.,0.);
	_MCneutrino_PDG = 0;
	    
	// 2 photons + the 6 leptons
	for(int i = 0 ; i < 8 ; i++)
	  {
	    MCParticle *mcpart = dynamic_cast<MCParticle*>(mccol->getElementAt(i));
	    TLorentzVector mcVec(TVector3(mcpart->getMomentum()),mcpart->getEnergy());
	
	    *_MCtotal += mcVec;
	
	    if(i < 2) *_MCphoton[i] = mcVec;
	
	    else 
	      {
		*_MCf[i-2] = mcVec;
		_MCfpdg[i-2] = mcpart->getPDG();
	      }
	  }
	    
	// get b quarks, lepton and neutrino information
	for(int i = 0 ; i < 6 ; i++)
	  {
	    if(_MCfpdg[i] == 5) *_MCb = *_MCf[i]; // b quark has charge -1/3

	    else if(_MCfpdg[i] == -5) *_MCbbar = *_MCf[i]; // bbar quark has charge +1/3
	
	    else if(abs(_MCfpdg[i]) == 11 || abs(_MCfpdg[i]) == 13 || abs(_MCfpdg[i]) == 15 ) // lepton
	      {
		*_MClepton = *_MCf[i];
		_MClepton_PDG = _MCfpdg[i];
	      }
	
	    else if(abs(_MCfpdg[i]) == 12 || abs(_MCfpdg[i]) == 14 || abs(_MCfpdg[i]) == 16 ) // neutrino
	      {
		*_MCneutrino = *_MCf[i];
		_MCneutrino_PDG = _MCfpdg[i];
	      }

	    else if(abs(_MCfpdg[i]) <= 4 && abs(_MCfpdg[i])%2 == 0) // first non-b quark
	      {
		*_MCq[0] = *_MCf[i];
		_MCqpdg[0] = _MCfpdg[i];
	      }

	    else if(abs(_MCfpdg[i]) <= 4 && abs(_MCfpdg[i])%2 == 1) // second non-b quark
	      {
		*_MCq[1] = *_MCf[i];
		_MCqpdg[1] = _MCfpdg[i];
	      }
	  }

	// Reconstruct the b-bbar system
	*_MCbbbar = *_MCb + *_MCbbar;
    
	// reconstruct the two W's
	*_MCW_jet = *_MCq[0] + *_MCq[1];
	*_MCW_lepton = *_MClepton + *_MCneutrino;

	// reconstruct tops
	if(_MClepton_PDG > 0) // i.e. lepton charge = -1
	  {
	    *_MCtop_lepton = *_MCbbar + *_MCW_lepton;
	    _MCtop_lepton_Q = -1;
	    *_MCtop_jet = *_MCb + *_MCW_jet;
	    _MCtop_jet_Q = +1;
	  }

	else if(_MClepton_PDG < 0) // i.e. lepton charge = +1
	  {
	    *_MCtop_lepton = *_MCb + *_MCW_lepton;
	    _MCtop_lepton_Q = +1;
	    *_MCtop_jet = *_MCbbar + *_MCW_jet;
	    _MCtop_jet_Q = -1;
	  }
	    
	if((TMath::Abs(_MCtop_lepton->M()-174)<5*1.51 || TMath::Abs(_MCtop_jet->M()-174)<5*1.51 || TMath::Sqrt(pow(_MCtop_jet->M()-174,2)+pow(_MCtop_lepton->M()-174,2))<(15*1.51)) || TMath::Abs(_MCbbbar->M()-91)>25)
	  {
	    _isttsl = true; 
	  }
	else _iszwwsl = true;
      }
      
    

    // Write PFOs information
    _numberOfPFOs = 0;
    LCCollection* col = evt->getCollection(_pfosCollection);
    _numberOfPFOs = col->getNumberOfElements();
    //    std:: cout << " Number of PFOs " << _numberOfPFOs << std::endl; 
    for(int i = 0 ; i < 7 ; i++)
      {
	_pfosFoxWolfram[i] = 0.;
	std::stringstream name;
	name << "FoxWolfram_moment(" << i << ")";
	_pfosFoxWolfram[i] = col->getParameters().getFloatVal(name.str());
      }
    _pfosSphericity = col->getParameters().getFloatVal("sphericity");
    _pfosAplanarity = col->getParameters().getFloatVal("aplanarity");
    _pfosOblateness = col->getParameters().getFloatVal("Oblateness");
    _pfosThrust = col->getParameters().getFloatVal("principleThrustValue");

    // Write myPFOs information
    col = evt->getCollection(_myPfosCollection);
    // Calculate total 4-momentum (WO lepton)
    _totalWOlepton->SetPxPyPzE(0.,0.,0.,0.);
    for(int i = 0 ; i < col->getNumberOfElements() ; i++)
      {
	ReconstructedParticle *a_part = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
	TLorentzVector Partvec(TVector3(a_part->getMomentum()),a_part->getEnergy());
	*_totalWOlepton += Partvec;
      }

    // Initialise and write : LEPTON
    _myLepton->SetPxPyPzE(0.,0.,0.,0.);
    LCCollection* leptonCol = evt->getCollection(_myLeptonCollection);
    _nbOfLep = 0;
    _nbOfLep = leptonCol->getParameters().getIntVal("NumberOfLeptonsCandidates");
    _lepPDG = 0 ;
    _lepPDG = leptonCol->getParameters().getIntVal("SelectedLeptonPDG");
    _lepCharge = 0.;
    _isGoodLep = false;
    _isLepFromTau = false;
    if(leptonCol->getParameters().getIntVal("isGoodLeptonFound") == 1) _isGoodLep = true;
    if(leptonCol->getParameters().getIntVal("isLeptonFromTau") == 1) _isLepFromTau = true;
    for(int i = 0 ; i < (int) leptonCol->getNumberOfElements() ; i++)
      {
	ReconstructedParticle *mylepton = dynamic_cast<ReconstructedParticle*>(leptonCol->getElementAt(i));
	TLorentzVector lepVec(TVector3(mylepton->getMomentum()),mylepton->getEnergy());
	_lepCharge += mylepton->getCharge();
	*_myLepton += lepVec;
      }

    // Initialise and write : JETS
    LCCollection* jetCol = evt->getCollection(_myJetsCollection);
    //    for(int i = 0 ; i < 4 ; i++)
    for(int i = 0 ; i < nJETS ; i++)
      {
	_jets[i]->SetPxPyPzE(0.,0.,0.,0.);
	_jetsBtag[i] = 0.;
	_jetsCtag[i] = 0.;
	_jetsBCtag[i] = 0.;
	_jetsTrueFlavour[i] = 0;
	_jetsTrueCharge[i] = 0;
      }
    _yplus = 0.;
    _yminus = 0.;
    _isGoodBs = false;
    _isGoodComb = false;
    
    if(jetCol->getNumberOfElements() == 6)
      {
	_yplus = jetCol->getParameters().getFloatVal("YPlus");
	_yminus = jetCol->getParameters().getFloatVal("YMinus");
	//get the PIDHandler fot this jet collection
	PIDHandler pidh(jetCol);
	//get the algorithm id for the Flavour Tag info
	int alid = pidh.getAlgorithmID("lcfiplus");
	int ibtag = pidh.getParameterIndex(alid,"BTag");
	int ictag = pidh.getParameterIndex(alid,"CTag");
	int ibctag = pidh.getParameterIndex(alid,"BCTag");
	//get the algorithm id for the MC truth info
	int mcid = pidh.getAlgorithmID("MCTruth");
	// intermediate variables
	TLorentzVector jets[4];
	float bTag[4], cTag[4], cTagBBack[4];
	int jetType[4], jetCharge[4];
	
	for(int i = 0 ; i < nJETS ; i++)
	  {
	    ReconstructedParticle *a_jet = dynamic_cast<ReconstructedParticle*>(jetCol->getElementAt(i));
	    TLorentzVector Jet4vec(TVector3(a_jet->getMomentum()),a_jet->getEnergy());
	    jets[i] = Jet4vec;
	    //get the Particle id object containing the Flavour Tag info
	    const ParticleID& pid = pidh.getParticleID(a_jet,alid);
	    //get the actual flavour tags
	    bTag[i] = pid.getParameters()[ibtag];
	    cTag[i] = pid.getParameters()[ictag];
	    cTagBBack[i] = pid.getParameters()[ibctag];
	    //get the Particle id object containing the MC info
	    const ParticleID& mcpid = pidh.getParticleID(a_jet,mcid);
	    //get the parameters for the MC info
	    FloatVec mcparams = mcpid.getParameters();
	    //get the true jet flavour
	    jetType[i] = (int) mcparams[pidh.getParameterIndex(mcid,"TrueJetFlavour")];	
	    jetCharge[i] = (int) mcparams[pidh.getParameterIndex(mcid,"TruePartonCharge")];
	  }
	

	// sort jets by btag
	int sortPositions[4];
	float tempBtag[4];
	for(int i = 0 ; i < nJETS ; i++) 
	  {
	    tempBtag[i] = bTag[i];
	    sortPositions[i] = i;
	  }
	for(int i = 0 ; i < nJETS ; i++)
	  {
	    int highestBtagPosition = i;
	    for(int j = i ; j < 4 ; j++)
	      {
		if(tempBtag[j] > tempBtag[highestBtagPosition]) highestBtagPosition = j;
	      }
	    float temp = tempBtag[i];
	    tempBtag[i] = tempBtag[highestBtagPosition];
	    tempBtag[highestBtagPosition] = temp;
	    int tempInt = sortPositions[i];
	    sortPositions[i] = sortPositions[highestBtagPosition];
	    sortPositions[highestBtagPosition] = tempInt;
	  }
	
	// write sorted jets in the tree
	for(int i = 0 ; i < nJETS ; i++)
	  {
	    *_jets[i] = jets[sortPositions[i]];
	    _jetsBtag[i] = bTag[sortPositions[i]];
	    _jetsCtag[i] = cTag[sortPositions[i]];
	    _jetsBCtag[i] = cTagBBack[sortPositions[i]];
	    _jetsTrueFlavour[i] = jetType[sortPositions[i]];
	    _jetsTrueCharge[i] = jetCharge[sortPositions[i]];
	  }
  
	if(_jetsTrueFlavour[0] == 5 && _jetsTrueFlavour[1] == 5) _isGoodBs = true;

	// Reconstruct the top
	*_Whad = *_jets[2] + *_jets[3];
	_selectedBjet = -1;
	_top_had->SetPxPyPzE(0.,0.,0.,0.);
	TLorentzVector top1 = *_Whad + *_jets[0];
	TLorentzVector top2 = *_Whad + *_jets[1];
	Double_t bpstar1 = top1.Gamma()*_jets[1]->Energy()*(1 - top1.Beta()*TMath::Cos(top1.Angle( - _jets[1]->Vect())));
	Double_t bpstar2 = top2.Gamma()*_jets[0]->Energy()*(1 - top2.Beta()*TMath::Cos(top2.Angle( - _jets[0]->Vect())));

	if(pow(top1.M()-174,2)/pow(6.3,2) + pow(top1.E()-250,2)/pow(8.0,2) + pow(bpstar1-69,2)/pow(10.0,2) < pow(top2.M()-174,2)/pow(6.3,2) + pow(top2.E()-250,2)/pow(8.0,2) + pow(bpstar2-69,2)/pow(10.0,2))
	  {
	    *_top_had = top1;
	    _selectedBjet = 0;
	    if(_isGoodBs && _jetsTrueCharge[0] == _MCtop_jet_Q) _isGoodComb = true;
	  }
	else
	  {
	    *_top_had = top2;
	    _selectedBjet = 1;
	    if(_isGoodBs && _jetsTrueCharge[1] == _MCtop_jet_Q) _isGoodComb = true;
	  }
      }
    _tree->Fill();
  }
  
  catch(DataNotAvailableException &e) { }
  
}



void ttbarTreeWriter::check(LCEvent * evt)
{ 
  
}


void ttbarTreeWriter::end()
{ 
  std::cout<< " ---------  ttbar Tree Writer Processor End ------  " << std::endl;   
  ROOTfile->Write();
  ROOTfile->Close();
  delete ROOTfile;
}
