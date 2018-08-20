#ifndef ttbar_h
#define ttbar_h 1

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
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"

using namespace lcio;
using namespace marlin;


/** 
 *  <h4>Output</h4> 
 *  A ROOT Tree containing information of desired Particles and MCParticles for semi-leptonic ttbar event
 * 
 * @author Philippe Doublet, LAL modified by Jeremy ROUENE
 * @version $Id: ttbar.h,v 1.0 2012/05/30 $ 
 *
 */

class ttbar : public Processor 
{
  
 public:
  virtual Processor*  newProcessor() { return new ttbar; }
  ttbar();
  virtual void init();
  virtual void processRunHeader( LCRunHeader* run );
  virtual void processEvent( LCEvent * evt ); 
  virtual void check( LCEvent * evt ); 
  virtual void end();
  
 protected:
  int _nRun;
  int _nEvt;
  int _eventnumber;
  int _runnumber;
  enum{
    nJETS=6  
};
  TFile* ROOTfile;
  std::string _hfilename;
  TTree* _tree;

  std::string  _myJetsCollection;
  std::string  _myvtxJetsCollection;
  std::string  _myrelJetsCollection;
  std::string  _mcCollection;

  TString *_Process;
  float _xsec;

  bool _isttsl;
  bool _iszwwsl;

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
  int nbqs;
  int ncqs;
  int nws;
  int nwhs;
  int nwls;

  int nqs;
  int nbs;

  bool isTop;
  bool ismyTop;
  bool isTop1;
  bool isTop2;


 bool B0Meson;
 bool BpMeson;
 bool BmMeson;


 int nbjets;

  int ntops;

    float WMass;
    float TopMass;

   float testop1Mass;
    float testop2Mass;
    float testop3Mass;
    float testop4Mass;

    float costheta;






  TLorentzVector *_jets[nJETS];
  TLorentzVector jets[nJETS];
  TLorentzVector bs[2];

  TLorentzVector mcbVec;
  TLorentzVector mcabVec;


  TLorentzVector mccVec;
  TLorentzVector mcacVec;

  TLorentzVector mctopVec;
  TLorentzVector mcatopVec;

  TLorentzVector B0MesonVec;
  TLorentzVector BmMesonVec;
  TLorentzVector BpMesonVec;
  TLorentzVector vtxpVec;


  TVector3 mcbmom;
  TVector3 mcabmom;
  //  TLorentzVector jets34;
  //  TLorentzVector jets35;
  //  TLorentzVector jets36;
  //  TLorentzVector jets45;
  //  TLorentzVector jets46;
  //  TLorentzVector jets56;
  //  TLorentzVector top1;
  //  TLorentzVector top2;
  //   TLorentzVector t1t2;

  float nforward;
  float nback;


  float _jetsBtag[nJETS];
  float _jetsCtag[nJETS];
  float _jetsUtag[nJETS];
  float _jetsDtag[nJETS];
  float _jetsBCtag[nJETS];
  float _jetsCharge[nJETS];
  float _jetsMass[nJETS];
  float _jetsMass34;
  float _jetsMass35;
  float _jetsMass36;
  float _jetsMass45;
  float _jetsMass46;
  float _jetsMass56;

  float _W1mass;
  float _W2mass;
  float _t1mass;
  float _t2mass;
  float _t1t2mass;



  float _WW1mass;
  float _WW2mass;
  float _Wt1mass;
  float _Wt2mass;


  float _Wtx1mass;
  float _Wtx2mass;

  float diffWt1;
  float diffWt2;
  float diffWt3;
  float diffWt4;

  float qptsumoptf1v;
  float qptsumoptf2v;


  float qmcb;
  float qmcab;
  
  int _jetsTrueFlavour[nJETS];
  int _jetsTrueCharge[nJETS];

  float bTag[nJETS];    
  float cTag[nJETS];    
  float uTag[nJETS];    
  float dTag[nJETS];
  float	cTagBBack[nJETS];
  float jetCharge[nJETS];    
  float jetMass[nJETS];    


  float xbTag[4];    

  float b1charge[4];
  float testpartu2q[6];
  float damtag[6];

  float bjetsq[2];
  float bkTag[2];

  float b1jettag;
  float b2jettag;


  float b1jetvq;
  float b2jetvq;
  float mcbenergx;
  float mcbchargex;


  TH1F* mtforchi;
  TH1F* etforchi;
  TH1F* pbforchi;

  TH1F* mcbenergy;

  TH1F* cost1t2;

  TH1F* cosbmcbrec;
  TH1F* coscmcbrec;

  TH2F* coscmcbrecvsbmc2;
  TH2F* coscmcbrecvsbmc22;
  TH2F* coscmcbrecvsbmc32;
  TH2F* coscmcbrecvsbmc42;
  TH2F* coscmcbrecvsbmcc2;
  TH2F* coscmcbrecvsbmc2c2;
  TH2F* coscmcbrecvsbmc3c2;
  TH2F* coscmcbrecvsbmc4c2;


  TH2F* coscmcbrecvsbmc1;
  TH2F* coscmcbrecvsbmc21;
  TH2F* coscmcbrecvsbmc31;
  TH2F* coscmcbrecvsbmc41;
  TH2F* coscmcbrecvsbmcc1;
  TH2F* coscmcbrecvsbmc2c1;
  TH2F* coscmcbrecvsbmc3c1;
  TH2F* coscmcbrecvsbmc4c1;



  TH2F* cosbreccmcbtag2;
  TH2F* cosbreccmcbtagc2;

  TH2F* cosbreccmcctag2;
  TH2F* cosbreccmcctagc2;

  TH2F* cosbreccmcbtag1;
  TH2F* cosbreccmcbtagc1;

  TH2F* cosbreccmcctag1;
  TH2F* cosbreccmcctagc1;


  TH1F* cosb1mcb2mc;
  TH1F* cosb1mccmc;
  TH1F* cosb2mccmc;

  TH1F* cosb1mcb2mccut;
  TH1F* cosb1mccmccut;
  TH1F* cosb2mccmccut;

  TH2F* cosb1mcb2mcb1mccmc;
  TH2F* cosb1mcb2mcb2mccmc;


  //  TH1F* cosb1mcb2mc;


  TH1F* btagsusp2; 
  TH1F* btagsusp1; 

  TH1F* costmctrecdiff;
  TH1F* costmctrecdiffbc;

  TH1F* costmctrecdiff1;
  TH1F* costmctrecdiff2;
  TH1F* costmctrecdiff3;
  TH1F* costmctrecdiff4;
  TH1F* costmctrecdiff5;


  TH1F* b1jetvq_b;
  TH1F* b2jetvq_b;

  TH1F* b1jetvq_ab;
  TH1F* b2jetvq_ab;

  TH1F* QSum_mcbp;
  TH1F* QSum_mcabp;

  TH1F* QPtSum_mcbp;
  TH1F* QPtSum_mcabp;


  TH1F* gud_bjbmcangle;
  TH1F* bad_bjbmcangle;
  //  TH1F* QPtSum_mcabp;





  float top1Mass;
  float top2Mass;
  float ttbarMass;


  float jetMass34;    
  float jetMass35;    
  float jetMass36;    
  float jetMass45;    
  float jetMass46;    
  float jetMass56;    


  float difff3456;
  float difff3546;
  float difff3645;


  float vtxpartmom;    
  float vtxpartmomsum;    
  float bbmass; 

  
  //  float bpmass;
  //  float bmmass;

  float qptsum;
  float qptsumb1;
  float qptsumb2;

  float  mccostheta ;
  float mccosinustheta;
  //  float  mccostheta ;

  int jetType[nJETS]; 

  //jetCharge[nJETS];

  float _costhetat;
  float _costhetatbar;

  float _costhetatq;
  float _costhetatbarq;



    TLorentzVector MCTopf1;
    TLorentzVector MCTopf2;
    TLorentzVector MCTopf1b;
    TLorentzVector MCTopf2b;



  TH1F*  gudbjmom;
  TH1F*  badbjmom;


  TH1F*  gudbjtag;
  TH1F*  badbjtag;


  TH1F* h_vtxmom;
  TH1F* h_vtxpartmom;
  TH1F* h_vtxpartmomsum;
  
  TH1F* h_vtxcharge;


  TH1F* h_btag;
  TH1F* h_btag1;
  TH1F* h_btag2;

  TH1F* h_btagz;
  TH1F* h_btagz1;

  TH1F* h_costheta;
  TH1F* h_costheta1;
  TH1F* h_costheta2;
  TH1F* h_costhetaavg;
  TH2F* h_costheta2d;
  TH1F* h_Chargesum;
  TH1F* h_Chargesum2;
  TH1F* h_Chargesum3;




  TH1F* h_mccosthetatop;
  TH1F* h_mccosthetatopbar;

  TH1F* h_costhetatop;
  TH1F* h_costhetatopbar;

  TH1F* h_btagnew;
  TH1F* h_ctag;
  TH1F* h_top1mass;
  TH1F* h_top2mass;
  TH1F* h_top1massuncut;
  TH1F* h_top2massuncut;
  TH1F* h_wtop1mass;
  TH1F* h_wtop2mass;
  TH1F* h_w1mass;
  TH1F* h_w2mass;
  TH1F* h_top12mass;

  TH1F* h_BmesonQ;
  TH1F* h_BmmesonQ;
  TH1F* h_BpmesonQ;

  TH1F* h_BmesonQnoB;
  TH1F* h_BmmesonQnoB;
  TH1F* h_BpmesonQnoB;


  TH1F* h_BmpmesonQ;


  //  TH1F* h_BpartonQ;
  TH1F* h_bplusartonQ;
  TH1F* h_bminuspartonQ;

  TH1F* h_bplusartonQ2;
  TH1F* h_bminuspartonQ2;
  //  TH1F* h_bpmpartonQ;


  TH2F* h_Bangle;
  TH2F* h_Bpangle;
  TH2F* h_Bmangle;



  TH1F* topangle;
  TH1F* topangle_bc;
  TH1F* topangle_mc;
  TH1F* topangle_fmc;
  TH1F* topangle_fmc2;
  TH1F* topangle_fmc3;

  TH1F* topangle_diff;

  TH1F* topangle_c1bc;
  TH1F* topangle_c1;
  TH1F* topangle_c2bc;
  TH1F* topangle_c2;

  TH1F* b1jetvq_c1;
  TH1F* b1jetvq_c2;
  TH1F* b2jetvq_c1;
  TH1F* b2jetvq_c2;

  TH1F* b1jetvq_c1bc;
  TH1F* b1jetvq_c2bc;
  TH1F* b2jetvq_c1bc;
  TH1F* b2jetvq_c2bc;



  TH1F* topquarkangle;
  TH1F* topbarquarkangle;

  TH1F* topbarquarkangleq_c;

  TH1F* topquarkangleq;
  TH1F* topbarquarkangleq;

  TH1F* topangleq;
  TH1F* topangleq_c;


  TH1F* topquarkangleqpt;
  TH1F* topbarquarkangleqpt;

  //  TH1F* h_bplusartonQ;


  TH2F* qmcb_b;
  TH2F* qmcab_ab;

  TH1F* QSum;


  TH1F* Wrong_Q;

  TH1F* QSum2;
  TH1F* QSum3;
  TH1F* QSum4;

  TH1F* QPtSum;
  TH1F* QPtSumb;
  TH1F* QPtSumbbar;
  TH1F* QPtSum0V;
  TH1F* QPtSum1V;
  TH1F* QPtSum2V;

  TH1F* QPtSum1V_S;
  TH1F* QPtSum2V_S;

  TH1F* QPtSumVtx1;
  TH1F* QPtSumVtx2;

  TH1F* PtSumVtx1;
  TH1F* PtSumVtx2;



  TH1F* QPtSumb_C;
  TH1F* QPtSumbbar_C;

  TH2F* h_top1top2mass;
  TH2F* h_top1top2mass_sel;


  TH1F* NMCTops;
  TH1F* NMCBqs;
  TH1F* NMCWs;
  TH1F* NMCWfs;
  TH1F* NMCWfmass;

  TH1F* NMCTopfmass;
  TH1F* NMCTopfmass_sum;
  TH1F* NMCTopftmass;

  TH1F* NMCW1theta;
  TH1F* NMCW2theta;



  TH2F* btagvstheta1;
  TH2F* btagvstheta2;
  TH2F* btagvstheta3;

  TH2F* btagvstheta1c;
  TH2F* btagvstheta2c;

  TH2F* btagvsen1;
  TH2F* btagvsen2;

  //  TH2F* btagvstheta3;


  TH2F* NMCTop12fmass;
  TH2F* NMCTop12fmass_sel;
  TH2F* NMCTop12fmass_sel1;
  TH2F* NMCTop12fmass_sel2;

  TH2F* NMCW12fenergy;
  TH2F* NMCW1Etheta;
  TH2F* NMCW2Etheta;
  TH2F* recW12fenergy;
  
  TH2F* NMCW12fmass;
  TH1F* NMCbbfmass;
  TH1F* NMCbbfmasszww;
  TH1F* NMCbbfmasstt;
 
  TH1F* NRecTops;

  TH1F* NVtx;
  TH1F* NEvts;

  TH1F* Topanglecomp1;
  TH1F* Topanglecomp2;
  TH1F* Topanglecomp3;
  TH1F* Topanglecomp4;
  TH1F* Topanglecomp1c;
  TH1F* Topanglecomp2c;
  TH1F* Topanglecomp3c;
  TH1F* Topanglecomp4c;

  TH1F* Topmasseff;
  TH1F* Topmasseff1;
  TH1F* Topmasseff2;
  TH1F* Topmasseff3;
  TH1F* Topmasseff4;
  TH1F* Topmasseff5;

  TH1F* Topmasseff1t;
  TH1F* Topmasseff2t;
  TH1F* Topmasseff3t;
  TH1F* Topmasseff4t;
  TH1F* Topmasseff5t;

  TH1F* Topmasseffchi200;
  TH1F* Topmasseffchi100;
  TH1F* Topmasseffchi50;
  TH1F* Topmasseffchi40;
  TH1F* Topmasseffchi30;

  TH1F* Topmasseffchi200bt;
  TH1F* Topmasseffchi100bt;
  TH1F* Topmasseffchi50bt;
  TH1F* Topmasseffchi40bt;
  TH1F* Topmasseffchi30bt;

  TH1F* Topmasseffchi200t;
  TH1F* Topmasseffchi100t;
  TH1F* Topmasseffchi50t;
  TH1F* Topmasseffchi40t;
  TH1F* Topmasseffchi30t;


  TH2F* topmasschia;
  TH2F* topmasschib;

  TH2F* topmasschiaE;
  TH2F* topmasschibE;

  TH2F* topmasschiaP;
  TH2F* topmasschibP;

  TH1F* b1mom;
  TH1F* b1momsel;

  TH1F* b1en;
  TH1F* b2en;

  TH1F* w1en;
  TH1F* w2en;

  TH2F* top1momvsmass;

  TH1F* h_mass34;
  TH1F* h_mass35;
  TH1F* h_mass36;

  TH1F* h_bbmass;
  TH1F* h_bbmass_tag;
  TH1F* h_bbmass_notag;
  TH1F* h_bbmass_sel;

  TH1F* h_difff3456;
  TH1F* h_difff3546;
  TH1F* h_difff3645;

      float difftx1, difftx2;


};

#endif
