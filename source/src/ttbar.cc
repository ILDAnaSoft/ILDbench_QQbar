#include "ttbar.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <TStyle.h>
#include <LCIOSTLTypes.h>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/PIDHandler.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ParticleID.h>
#include <TH1.h>
#include <TH2.h>
#include <EVENT/LCRelation.h>
#include <EVENT/Vertex.h>
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;
using namespace std;

ttbar attbar;

ttbar::ttbar() : Processor("ttbar")
{
  _description = "Analyse ttbar to fully hadronic mode at Reconstruction and MC level";
  
  registerProcessorParameter( "ROOTFileName","Output ROOT File Name",
                              _hfilename,std::string("ttbar.root"));

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, "RefinedJets" ,
                           "Name of the Jets collection"  ,
                           _myJetsCollection, std::string("myJets"));

  registerInputCollection( LCIO::LCRELATION, "RefinedJets_rel" ,
                           "Name of the Vertex Jets relation"  ,
                           _myrelJetsCollection, std::string("myrelJets"));
  
  registerInputCollection( LCIO::MCPARTICLE, "MCCollection" ,
                           "MC collection (MCParticle for Mokka, MCParticlesSkimmed for DST)" ,
                           _mcCollection, std::string("MCParticle"));
}


void ttbar::init()
{ 
  mccostheta=-999;
  nforward=0;
  nback =0;
  _nRun = 0;
  _nEvt = 0;
  _eventnumber=0;
  ntops=0;
  nbqs=0; 
  ncqs=0; 
  nws=0; 
  nwhs=0; 
  nwls=0; 
  nqs=0; 
 
 ROOTfile = new TFile(_hfilename.c_str(), "RECREATE", _hfilename.c_str());

     h_btag = new TH1F ("BTag", "BTag", 100, 0,1 );
     h_btag1 = new TH1F ("BTag2", "BTag2", 100, 0,1 );
     h_btag2 = new TH1F ("BTag3", "BTag3", 100, 0,1 ); 
     h_btagz = new TH1F ("BTagZ", "BTagZ", 100, 0,1 );
     h_btagz1 = new TH1F ("BTagZ2", "BTagZ2", 100, 0,1 );

     b1en = new TH1F ("b1energy", "b1energy", 100, 0,200 );
     b2en = new TH1F ("b2energy", "b2energy", 100, 0,200 );

     w1en = new TH1F ("w1energy", "w1energy", 100, 0,300 );
     w2en = new TH1F ("w2energy", "w2energy", 100, 0,300 );

     h_ctag = new TH1F ("CTag", "CTag", 100, 0,1 );
     h_top1mass = new TH1F ("Top1Mass", "Top1Mass", 100, 0,350);
     h_top2mass = new TH1F ("Top2Mass", "Top2Mass", 100, 0,350 );

     mtforchi = new TH1F ("mt", "mt", 100, 0,300);
     etforchi = new TH1F ("et", "et", 100, 0,350); 
     pbforchi = new TH1F ("pb", "pb", 100, 0,200); 

     cosbmcbrec=new TH1F ("CosbMCvsbRec", "CosbMCvsbRec", 100, -1,1 );
     coscmcbrec=new TH1F ("CoscMCvsbRec", "CoscMCvsbRec", 100, -1,1 );

     costmctrecdiff=new TH1F ("CostMCtRecdiff", "CostMCtRecdiff", 100, -1,1 );
     costmctrecdiffbc=new TH1F ("CostMCtRecdiffbc", "CostMCtRecdiffbc", 100, -1,1 ); 

     coscmcbrecvsbmc2=new TH2F ("CoscMCbRecbMC2D2", "CoscMCvsbRecbMC2D2", 100, -1,1, 100, -1,1 );
     coscmcbrecvsbmc22=new TH2F ("CoscMCbRecbMC2Dn2", "CoscMCvsbRecbMC2Dn2", 100, -1,1, 100, -1,1 );
     coscmcbrecvsbmc32=new TH2F ("CoscMCbRecbMC2Dnt2", "CoscMCvsbRecbMC2Dnt2", 100, -1,1, 100, -1,1 );
     coscmcbrecvsbmc42=new TH2F ("CoscMCbRecbMC2Dntt2", "CoscMCvsbRecbMC2Dntt2", 100, -1,1, 100, -1,1 );

     coscmcbrecvsbmcc2=new TH2F ("CoscMCbRecbMC2Dc2", "CoscMCvsbRecbMC2Dc2", 100, -1,1, 100, -1,1 );
     coscmcbrecvsbmc2c2=new TH2F ("CoscMCbRecbMC2Dnc2", "CoscMCvsbRecbMC2Dnc2", 100, -1,1, 100, -1,1 );
     coscmcbrecvsbmc3c2=new TH2F ("CoscMCbRecbMC2Dntc2", "CoscMCvsbRecbMC2Dntc2", 100, -1,1, 100, -1,1 );
     coscmcbrecvsbmc4c2=new TH2F ("CoscMCbRecbMC2Dnttc2", "CoscMCvsbRecbMC2Dnttc2", 100, -1,1, 100, -1,1 );


     coscmcbrecvsbmc1=new TH2F ("CoscMCbRecbMC2D1", "CoscMCvsbRecbMC2D1", 100, -1,1, 100, -1,1 );
     coscmcbrecvsbmc21=new TH2F ("CoscMCbRecbMC2Dn1", "CoscMCvsbRecbMC2Dn1", 100, -1,1, 100, -1,1 );
     coscmcbrecvsbmc31=new TH2F ("CoscMCbRecbMC2Dnt1", "CoscMCvsbRecbMC2Dnt1", 100, -1,1, 100, -1,1 );
     coscmcbrecvsbmc41=new TH2F ("CoscMCbRecbMC2Dntt1", "CoscMCvsbRecbMC2Dntt1", 100, -1,1, 100, -1,1 );

     coscmcbrecvsbmcc1=new TH2F ("CoscMCbRecbMC2Dc1", "CoscMCvsbRecbMC2Dc1", 100, -1,1, 100, -1,1 );
     coscmcbrecvsbmc2c1=new TH2F ("CoscMCbRecbMC2Dnc1", "CoscMCvsbRecbMC2Dnc1", 100, -1,1, 100, -1,1 );
     coscmcbrecvsbmc3c1=new TH2F ("CoscMCbRecbMC2Dntc1", "CoscMCvsbRecbMC2Dntc1", 100, -1,1, 100, -1,1 );
     coscmcbrecvsbmc4c1=new TH2F ("CoscMCbRecbMC2Dnttc1", "CoscMCvsbRecbMC2Dnttc1", 100, -1,1, 100, -1,1 );


     
     cosbreccmcbtag2 =new TH2F ("CoscMCbRecBtag2", "CoscMCvsbRecBtag2", 100, 0,1, 100, -1,1 );
     cosbreccmcbtagc2 =new TH2F ("CoscMCbRecBtagc2", "CoscMCvsbRecBtagc2", 100, 0,1, 100, -1,1 );

     cosbreccmcctag2 =new TH2F ("CoscMCbRecCtag2", "CoscMCvsbRecCtag2", 100, 0,1, 100, -1,1 );
     cosbreccmcctagc2 =new TH2F ("CoscMCbRecCtagc2", "CoscMCvsbRecCtagc2", 100, 0,1, 100, -1,1 );

     cosbreccmcbtag1 =new TH2F ("CoscMCbRecBtag1", "CoscMCvsbRecBtag1", 100, 0,1, 100, -1,1 );
     cosbreccmcbtagc1 =new TH2F ("CoscMCbRecBtagc1", "CoscMCvsbRecBtagc1", 100, 0,1, 100, -1,1 );

     cosbreccmcctag1 =new TH2F ("CoscMCbRecCtag1", "CoscMCvsbRecCtag1", 100, 0,1, 100, -1,1 );
     cosbreccmcctagc1 =new TH2F ("CoscMCbRecCtagc1", "CoscMCvsbRecCtagc1", 100, 0,1, 100, -1,1 );
 

     btagsusp2=new TH1F ("BTagSusp2", "BTagSusp2", 100, 0,1 ); 
     btagsusp1=new TH1F ("BTagSusp1", "BTagSusp1", 100, 0,1 ); 
    

     cosb1mcb2mc =new TH1F ("cosb1mcb2mc", "cosb1mcb2mc", 100, -1,1 );
     cosb1mccmc =new TH1F ("cosb1mccmc", "cosb1mccmc", 100, -1,1 );
     cosb2mccmc =new TH1F ("cosb2mccmc", "cosb2mccmc", 100, -1,1 );

     cosb1mcb2mccut =new TH1F ("cosb1mcb2mccut", "cosb1mcb2mccut", 100, -1,1 );
     cosb1mccmccut =new TH1F ("cosb1mccmccut", "cosb1mccmccut", 100, -1,1 );
     cosb2mccmccut =new TH1F ("cosb2mccmccut", "cosb2mccmccut", 100, -1,1 );

     cosb1mcb2mcb1mccmc  =new TH2F ("Cosb1mcb2mcb1mccms", "Cosb1mcb2mcb1mccmc", 100, -1,1, 100, -1,1 );
     cosb1mcb2mcb2mccmc  =new TH2F ("Cosb1mcb2mcb2mccms", "Cosb1mcb2mcb2mccmc", 100, -1,1, 100, -1,1 );

     costmctrecdiff1=new TH1F ("CostMCtRecdiff1", "CostMCtRecdiff1", 100, -1,1 );
     costmctrecdiff2=new TH1F ("CostMCtRecdiff2", "CostMCtRecdiff2", 100, -1,1 );
     costmctrecdiff3=new TH1F ("CostMCtRecdiff3", "CostMCtRecdiff3", 100, -1,1 );
     costmctrecdiff4=new TH1F ("CostMCtRecdiff4", "CostMCtRecdiff4", 100, -1,1 );
     costmctrecdiff5=new TH1F ("CostMCtRecdiff5", "CostMCtRecdiff5", 100, -1,1 );

     h_costheta = new TH1F ("CosTheta", "Costheta", 100, -1,1.1);
     h_costheta1 = new TH1F ("CosTheta1", "Costheta1", 100, -1,1.1);
     h_costheta2 = new TH1F ("CosTheta2", "Costheta2", 100, -1,1.1);
     h_costhetaavg = new TH1F ("CosThetaavg", "Costhetaavg", 100, -1,1.1);
     h_costheta2d= new TH2F ("CosTheta2d", "Costheta2d", 100, -1,1.1, 100 ,-1,1.1);

     h_Chargesum=  new TH1F ("CargeSum", "ChargeSum", 10, -5,5);
     h_Chargesum2=  new TH1F ("CargeSumVtx", "ChargeSumVtx", 10, -5,5);
     h_Chargesum3=  new TH1F ("CargeSumVtx_rp", "ChargeSumVtx_rp", 10, -5,5);
     h_BmesonQ = new TH1F ("B0_Mesons_q", "B0_Mesons_q", 10, -5,5);
     h_BmmesonQ = new TH1F ("Bm_Mesons_q", "Bm_Mesons_q", 10, -5,5);
     h_BpmesonQ = new TH1F ("Bp_Mesons_q", "Bp_Mesons_q", 10, -5,5);

     h_BmesonQnoB = new TH1F ("B0_Mesons_qnoB", "B0_Mesons_qnoB", 10, -5,5);
     h_BmmesonQnoB = new TH1F ("Bm_Mesons_qnoB", "Bm_Mesons_qnoB", 10, -5,5);
     h_BpmesonQnoB = new TH1F ("Bp_Mesons_qnoB", "Bp_Mesons_qnoB", 10, -5,5);

     h_BmpmesonQ = new TH1F ("Bp_Bm_Mesons_q", "Bp_Bm_Mesons_q", 10, -5,5);
     
     h_bplusartonQ = new TH1F ("bp_parton_q", "bp_parton_q", 10, -5,5);
     h_bminuspartonQ = new TH1F ("bm_parton_q", "bm_parton_q", 10, -5,5);

     h_bplusartonQ2 = new TH1F ("bp_parton_q2", "bp_parton_q2", 10, -5,5);
     h_bminuspartonQ2 = new TH1F ("bm_parton_q2", "bm_parton_q2", 10, -5,5);


     h_wtop1mass = new TH1F ("WTop1Mass", "WTop1Mass", 100, 0,350);
     h_wtop2mass = new TH1F ("WTop2Mass", "WTop2Mass", 100, 0,350);
     h_w1mass = new TH1F ("W1Mass", "W1Mass", 100, 0,210);
     h_w2mass = new TH1F ("W2Mass", "W2Mass", 100, 0,210);
     h_top1massuncut= new TH1F ("Top1MassUncut", "Top1MassUncut", 100, 0,400);
     h_top2massuncut= new TH1F ("Top2MassUncut", "Top2MassUncut", 100, 0,400);
     h_top12mass= new TH1F ("Top12Mass", "Top12Mass", 100, 0,600);
     h_top1top2mass= new TH2F ("Top1Top2Mass", "Top1Top2Mass", 100, 0,400, 100, 0,400); 
     h_top1top2mass_sel= new TH2F ("Top1Top2Mass_sel", "Top1Top2Mass_sel", 100, 0,400, 100, 0,400); 

     h_mccosthetatop =  new TH1F ("MCCosThetaT", "MCCosthetaT", 100, -1,1);

     h_mccosthetatopbar =  new TH1F ("MCCosThetaTbar", "MCCosthetaTbar", 100, -1,1);    

     h_vtxcharge =  new TH1F ("Vertex_Charge", "Vertex_Charge", 100, -3,3);    

     h_costhetatop=  new TH1F ("CosThetaTop", "CosthetaTop", 100, -1,1);
     h_costhetatopbar =  new TH1F ("CosThetaTopBar", "CosthetaTopBar", 100, -1,1);


     h_Bangle =  new TH2F ("B0angle", "B0angle", 10, -5,5, 1000, -1,1.1);
     h_Bpangle =  new TH2F ("Bpangle", "Bpangle", 10, -5,5, 1000, -1,1.1);
     h_Bmangle =  new TH2F ("Bmangle", "Bmangle", 10, -5,5, 1000, -1,1.1);

     topangle =  new TH1F ("TopAngle", "TopAngle", 100, -1,1);
     topangle_bc =  new TH1F ("TopAngle_bc", "TopAngle_bc", 100, -1,1);
     topangle_mc =  new TH1F ("TopAngle_mc", "TopAngle_mc", 100, -1,1);
     topangle_fmc =  new TH1F ("TopAngle_fmc", "TopAngle_fmc", 100, -1,1);

     topangle_fmc2 =  new TH1F ("TopAngle_fmc2", "TopAngle_fmc2", 100, -1,1);
     topangle_fmc3 =  new TH1F ("TopAngle_fmc3", "TopAngle_fmc3", 100, -1,1);

     topangle_diff =  new TH1F ("TopAngle_diff", "TopAngle_diff", 100, -1,1);
     topangle_c1 = new TH1F ("TopAngle_c1", "TopAngle_c1", 100, -1,1); 
     topangle_c1bc = new TH1F ("TopAngle_c1bc", "TopAngle_c1bc", 100, -1,1); 
     topangle_c2 = new TH1F ("TopAngle_c2", "TopAngle_c2", 100, -1,1); 
     topangle_c2bc = new TH1F ("TopAngle_c2bc", "TopAngle_c2bc", 100, -1,1); 
     
     topangleq =  new TH1F ("TopAngle_q", "TopAngle_q", 100, -1,1);
     topangleq_c =  new TH1F ("TopAngle_qc", "TopAngle_qc", 100, -1,1 );

     Topanglecomp1  =  new TH1F ("Topanglecomp1","Topanglecomp1", 100, -1,1);
     Topanglecomp2  =  new TH1F ("Topanglecomp2","Topanglecomp2", 100, -1,1);
     Topanglecomp3  =  new TH1F ("Topanglecomp3","Topanglecomp3", 100, -1,1);
     Topanglecomp4  =  new TH1F ("Topanglecomp4","Topanglecomp4", 100, -1,1);
     Topanglecomp1c  =  new TH1F ("Topanglecomp1c","Topanglecomp1c", 100, -1,1);
     Topanglecomp2c  =  new TH1F ("Topanglecomp2c","Topanglecomp2c", 100, -1,1);
     Topanglecomp3c  =  new TH1F ("Topanglecomp3c","Topanglecomp3c", 100, -1,1);    
     Topanglecomp4c  =  new TH1F ("Topanglecomp4c","Topanglecomp4c", 100, -1,1);    


     Topmasseff = new TH1F ("TopMassEff", "TopMassEff", 500, 0,550);
     Topmasseff1 = new TH1F ("TopMassEff1", "TopMassEff1", 500, 0,550);
     Topmasseff2 = new TH1F ("TopMassEff2", "TopMassEff2", 500, 0,550);
     Topmasseff3 = new TH1F ("TopMassEff3", "TopMassEff3", 500, 0,550);
     Topmasseff4 = new TH1F ("TopMassEff4", "TopMassEff4", 500, 0,550);
     Topmasseff5 = new TH1F ("TopMassEff5", "TopMassEff5", 500, 0,550);

     Topmasseff1t = new TH1F ("TopMassEff1t", "TopMassEff1t", 500, 0,550);
     Topmasseff2t = new TH1F ("TopMassEff2t", "TopMassEff2t", 500, 0,550);
     Topmasseff3t = new TH1F ("TopMassEff3t", "TopMassEff3t", 500, 0,550);
     Topmasseff4t = new TH1F ("TopMassEff4t", "TopMassEff4t", 500, 0,550);
     Topmasseff5t = new TH1F ("TopMassEff5t", "TopMassEff5t", 500, 0,550);

     Topmasseffchi200 = new TH1F ("TopMassEffchi200", "TopMassEffchi200", 500, 0,550);
     Topmasseffchi100 = new TH1F ("TopMassEffchi100", "TopMassEffchi100", 500, 0,550);
     Topmasseffchi50  = new TH1F ("TopMassEffchi50",  "TopMassEffchi50",  500, 0,550);
     Topmasseffchi40  = new TH1F ("TopMassEffchi40",  "TopMassEffchi40",  500, 0,550);
     Topmasseffchi30  = new TH1F ("TopMassEffchi30",  "TopMassEffchi30",  500, 0,550);


     Topmasseffchi200bt = new TH1F ("TopMassEffchi200bt", "TopMassEffchi200bt", 500, 0,550);
     Topmasseffchi100bt = new TH1F ("TopMassEffchi100bt", "TopMassEffchi100bt", 500, 0,550);
     Topmasseffchi50bt = new TH1F ("TopMassEffchi50bt",  "TopMassEffchi50bt",  500, 0,550);
     Topmasseffchi40bt  = new TH1F ("TopMassEffchi40bt",  "TopMassEffchi40bt",  500, 0,550);
     Topmasseffchi30bt  = new TH1F ("TopMassEffchi30bt",  "TopMassEffchi30bt",  500, 0,550);

     Topmasseffchi200t = new TH1F ("TopMassEffchi200t", "TopMassEffchi200t", 500, 0,550);
     Topmasseffchi100t = new TH1F ("TopMassEffchi100t", "TopMassEffchi100t", 500, 0,550);
     Topmasseffchi50t = new TH1F ("TopMassEffchi50t",  "TopMassEffchi50t",  500, 0,550);
     Topmasseffchi40t  = new TH1F ("TopMassEffchi40t",  "TopMassEffchi40t",  500, 0,550);
     Topmasseffchi30t  = new TH1F ("TopMassEffchi30t",  "TopMassEffchi30t",  500, 0,550);
     

 
     topmasschia=  new TH2F ("TopMassChia", "TopMassChia", 100, 0,400, 100,0,500 );
     topmasschib=  new TH2F ("TopMassChib", "TopMassChib", 100, 0,400, 100,0,1000 );

     topmasschiaE=  new TH2F ("TopMassChiaE", "TopMassChiaE", 100, 0,400, 100,0,500 );
     topmasschibE=  new TH2F ("TopMassChibE", "TopMassChibE", 100, 0,400, 100,0,1000 );

     topmasschiaP=  new TH2F ("TopMassChiaP", "TopMassChiaP", 100, 0,400, 100,0,500 );
     topmasschibP=  new TH2F ("TopMassChibP", "TopMassChibP", 100, 0,400, 100,0,1000 );

     topquarkangle =  new TH1F ("TopQuarkAngle", "TopQuarkAngle", 100, -1,1);
     topbarquarkangle = new TH1F ("TopBarQuarkAngle", "TopBarQuarkAngle", 100, -1,1 );

     topquarkangleq =  new TH1F ("TopQuarkAngle_Q", "TopQuarkAngle_Q", 100, -1,1 );
     topbarquarkangleq = new TH1F ("TopBarQuarkAngle_Q", "TopBarQuarkAngle_Q", 100, -1,1);


     topbarquarkangleq_c = new TH1F ("TopBarQuarkAngle_QC", "TopBarQuarkAngle_QC", 100, -1,1);

     topquarkangleqpt =  new TH1F ("TopQuarkAngle_QP", "TopQuarkAngle_QP", 100, -1,1 );
     topbarquarkangleqpt = new TH1F ("TopBarQuarkAngle_QP", "TopBarQuarkAngle_QP", 100, -1,1);
 
     cost1t2 = new TH1F ("Angle_t1t2", "Angle_t1t2", 100, -1,1);
     qmcb_b =  new TH2F ("qmcb_b", "qmcb_b", 10, -5,5,  10, -5,5);
     qmcab_ab =   new TH2F ("qmcab_ab","qmcab_ab", 10, -5,5,  10, -5,5);

     Wrong_Q= new TH1F ("WrongQ", "WrongQ", 100, -1,1);

     QSum= new TH1F ("QSum_sorted", "QSum_sorted", 10, -5,5);
     QSum2= new TH1F ("QSum_sorted2", "QSum_sorted2", 10, -5,5);
     QSum3= new TH1F ("QSum_sorted3", "QSum_sorted3", 10, -5,5);
     QSum4= new TH1F ("QSum_sorted4", "QSum_sorted4", 10, -5,5);

     h_vtxmom =  new TH1F ("VtxPt", "VtxPt", 100, 0,100);
     h_vtxpartmom =  new TH1F ("VtxPartPt", "VtxPartPt", 100, 0,40);
     h_vtxpartmomsum =  new TH1F ("VtxPartPtSum", "VtxPartPtSum", 100, 0,100);

     QPtSum  =  new TH1F ("QPtSum", "QPtSum", 100, -1.1, 1.1);
     QPtSumb  =  new TH1F ("QPtSumb", "QPtSumb", 100, -1.1, 1.1);
     QPtSumbbar  =  new TH1F ("QPtSumbbar", "QPtSumbar", 100, -1.1, 1.1);
     
     NMCTops   =  new TH1F ("NMCTops", "NMCTops", 10, 0  , 10);
     NMCBqs   =  new TH1F ("NMCBqs", "NMCBqs", 10, 0  , 10);
     NMCWs   =  new TH1F ("NMCWs", "NMCWs", 10, 0  , 10);
     NMCWfs   =  new TH1F ("NMCWfs", "NMCWfs", 10, -5  , 5);

     NMCWfmass   =  new TH1F ("NMCWfmass", "NMCWfmass", 160, 0 , 160);

     NMCTopfmass   =  new TH1F ("NMCTopfmass", "NMCTopfmass", 500, 0 , 400);

     NMCTopfmass_sum   =  new TH1F ("NMCTopfmass_sum", "NMCTopfmass_sum", 500, 0 , 500);

     NMCTopftmass   =  new TH1F ("NMCTopftmass", "NMCTopftmass", 500, 0 , 550);

     NMCTop12fmass   =  new TH2F ("NMCTop12fmass", "NMCTop12fmass", 500, 0 , 400, 500,0,400);
     NMCTop12fmass_sel   =  new TH2F ("NMCTop12fmass_sel", "NMCTop12fmass_sel", 500, 0 , 400, 500,0,400);
     NMCTop12fmass_sel1   =  new TH2F ("NMCTop12fmass_sel1", "NMCTop12fmass_sel1", 500, 0 , 400, 500,0,400);
     NMCTop12fmass_sel2   =  new TH2F ("NMCTop12fmass_sel2", "NMCTop12fmass_sel2", 500, 0 , 400, 500,0,400);

     NMCW12fmass   =  new TH2F ("NMCW12fmass", "NMCW12fmass", 160, 0 , 160, 160,0,160);
     NMCW12fenergy   =  new TH2F ("NMCW12fenergy", "NMCW12fenergy", 300, 0 , 400, 300,0,400);

     NMCW1theta   =  new TH1F ("NMCW1theta", "NMCW1theta", 100, -1 , 1);
     NMCW2theta   =  new TH1F ("NMCW2theta", "NMCW2theta", 100, -1 , 1);

     NMCW1Etheta   =  new TH2F ("NMCW1Etheta", "NMCW1Etheta", 100, -1 , 1, 300,0,400);
     NMCW2Etheta   =  new TH2F ("NMCW2Etheta", "NMCW2Etheta", 100, -1 , 1, 300,0,400);
     
     recW12fenergy   =  new TH2F ("recW12fenergy", "recW12fenergy", 300, 0 , 400, 300,0,400);

     NMCbbfmass   =  new TH1F ("NMCbbfmass", "NMCbbfmass", 500, 0 , 400);
     NMCbbfmasszww   =  new TH1F ("NMCbbfmasszww", "NMCbbfmasszww", 500, 0 , 400);
     NMCbbfmasstt   =  new TH1F ("NMCbbfmasstt", "NMCbbfmasstt", 500, 0 , 400);

     NRecTops   =  new TH1F ("NRecTops", "NRecTops", 10, 0  , 10);

     NVtx =   new TH1F ("NVtx", "NVtx", 10, 0  , 10);
     NEvts =   new TH1F ("NEvts", "NEvts", 10, 0,4);

     QPtSum0V  =  new TH1F ("QPtSum0V", "QPtSum0V", 100, -1.1, 1.1);
     QPtSum1V  =  new TH1F ("QPtSum1V", "QPtSum1V", 100, -1.1, 1.1);
     QPtSum2V  =  new TH1F ("QPtSum2V", "QPtSum2V", 100, -1.1, 1.1);


     QPtSum1V_S  =  new TH1F ("QPtSum1V_S", "QPtSum1V_S", 100, -1.1, 1.1);
     QPtSum2V_S  =  new TH1F ("QPtSum2V_S", "QPtSum2V_S", 100, -1.1, 1.1);

     QPtSumVtx1  =  new TH1F ("QPtSumVtx1", "QPtSumVtx1", 100, -1.1, 1.1);
     QPtSumVtx2  =  new TH1F ("QPtSumVtx2", "QPtSumVtx2", 100, -1.1, 1.1);

     PtSumVtx1  =  new TH1F ("PtSumVtx1", "PtSumVtx1", 100, 0,100);
     PtSumVtx2  =  new TH1F ("PtSumVtx2", "PtSumVtx2", 100, 0,100);

     QSum_mcbp = new TH1F ("QSum_mcbp", "QSum_mcbp", 10, -5,5);
     QSum_mcabp = new TH1F ("QSum_mcabp", "QSum_mcabp", 10, -5,5);

     b1jetvq_b = new TH1F ("b1jetvq_b ", "b1jetvq_b", 10, -5,5);
     b2jetvq_b =  new TH1F ("b2jetvq_b", "b2jetvq_b", 10, -5,5);

     b1jetvq_ab = new TH1F ("b1jetvq_ab ", "b1jetvq_ab", 10, -5,5);
     b2jetvq_ab =  new TH1F ("b2jetvq_ab", "b2jetvq_ab", 10, -5,5);

     b1jetvq_c1bc = new TH1F ("b1jetvq_c1bc ", "b1jetvq_c1bc", 10, -5,5);
     b2jetvq_c1bc =  new TH1F ("b2jetvq_c1bc", "b2jetvq_c1bc", 10, -5,5);

     b1jetvq_c2bc = new TH1F ("b1jetvq_c2bc ", "b1jetvq_c2bc", 10, -5,5);
     b2jetvq_c2bc =  new TH1F ("b2jetvq_c2bc", "b2jetvq_c2bc", 10, -5,5);

     b1jetvq_c1 = new TH1F ("b1jetvq_c1 ", "b1jetvq_c1", 10, -5,5);
     b2jetvq_c1 =  new TH1F ("b2jetvq_c1", "b2jetvq_c1", 10, -5,5);

     b1jetvq_c2 = new TH1F ("b1jetvq_c2 ", "b1jetvq_c2", 10, -5,5);
     b2jetvq_c2 =  new TH1F ("b2jetvq_c2", "b2jetvq_c2", 10, -5,5);

     QPtSumb_C  =  new TH1F ("QPtSumb_C", "QPtSumb_C", 100, -1.1, 1.1);
     QPtSumbbar_C  =  new TH1F ("QPtSumbbar_C", "QPtSumbbar_C", 100, -1.1, 1.1);
     
     QPtSum_mcbp  =  new TH1F ("QPtSum_mcbp", "QPtSum_mcbp", 100, -1.1, 1.1);
     QPtSum_mcabp  =  new TH1F ("QPtSum_mcabp", "QPtSum_mcabp", 100, -1.1, 1.1);

     gud_bjbmcangle  =  new TH1F ("gud_bjbmcangle", "gud_bjbmcangle", 100, -1, 1); 
     bad_bjbmcangle  =  new TH1F ("bad_bjbmcangle", "bad_bjbmcangle", 100, -1, 1); 

     gudbjmom =  new TH1F ("gud_bjmom", "bad_bjmom", 100, 0, 140);  
     badbjmom =  new TH1F ("bad_bjmom", "bad_bjmom", 100, 0, 140); 

     gudbjtag =  new TH1F ("gud_bjtag", "bad_bjtag", 100, 0, 1);  
     badbjtag =  new TH1F ("bad_bjtag", "bad_bjtag", 100, 0, 1); 

     b1mom =  new TH1F ("b1mom", "b1mom", 100, 0, 140);  
     b1momsel =  new TH1F ("b1momsel", "b1momsel", 100, 0, 140);  
     top1momvsmass =  new TH2F ("t1pvm", "t1pvm", 100, 100, 300, 100, 0,160);


     btagvstheta1 =  new TH2F ("BTagvsCosTheta1", "BTagvsCosTheta1", 100, -1, 1, 100, 0,1);
     btagvstheta2 =  new TH2F ("BTagvsCosTheta2", "BTagvsCosTheta2", 100, -1, 1, 100, 0,1);
     btagvstheta3 =  new TH2F ("BTagvsCosTheta3", "BTagvsCosTheta3", 100, -1, 1, 100, 0,1);


     btagvsen1 =  new TH2F ("BTagvsEn1", "BTagvsCosEn1", 100, 0, 250, 100, 0,1);
     btagvsen2 =  new TH2F ("BTagvsEn2", "BTagvsCosEn2", 100, 0, 250, 100, 0,1);

     mcbenergy  = new TH1F ("mcbenergy", "mcbenergy", 100, 0, 250);  
     
     h_mass34 = new TH1F ("mass34", "mass34", 100, 0, 400);  
     h_mass35 = new TH1F ("mass35", "mass35", 100, 0, 400);  
     h_mass36 = new TH1F ("mass36", "mass36", 100, 0, 400);  
     
     h_bbmass = new TH1F ("bbmass", "bbmass", 100, 0, 400);       
     h_bbmass_tag = new TH1F ("bbmass_tag", "bbmass_tag", 100, 0, 400);  
     h_bbmass_notag = new TH1F ("bbmass_notag", "bbmass_notag", 100, 0, 400);  
     h_bbmass_sel = new TH1F ("bbmass_sel", "bbmass_sel", 100, 0, 400);  
     
     h_difff3456 = new TH1F ("difff3456", "difff3456", 100, 0, 500);  
     h_difff3546 = new TH1F ("difff3546", "difff3546", 100, 0, 500);  
     h_difff3645 = new TH1F ("difff3645", "difff3645", 100, 0, 500);  
     
     _tree = new TTree("tree", "tree");
     
     
     for(int i = 0 ; i < nJETS ; i++){
       _jets[i] = new TLorentzVector();
       std::stringstream name;
       name << "jet" << i;
       _tree->Branch(name.str().c_str(),"TLorentzVector",&_jets[i],16000,0);
       _tree->Branch((name.str()+"_btag").c_str(), &_jetsBtag[i], (name.str()+"_btag/F").c_str());
       _tree->Branch((name.str()+"_ctag").c_str(), &_jetsCtag[i], (name.str()+"_ctag/F").c_str());
       _tree->Branch((name.str()+"_bctag").c_str(), &_jetsBCtag[i], (name.str()+"_bctag/F").c_str());
       _tree->Branch((name.str()+"_charge").c_str(), &_jetsCharge[i], (name.str()+"_charge/F").c_str());
       _tree->Branch((name.str()+"_mass").c_str(), &_jetsMass[i], (name.str()+"_mass/F").c_str());
     }
     
     _tree->Branch("runNumber", &_nRun, "runNumber/I");
     _tree->Branch("eventNumber", &_nEvt, "eventNumber/I");
     _tree->Branch("crossSection", &_xsec, "crossSection/F");
     
     _Process = new TString();
     _tree->Branch("Process","TString",&_Process,16000,0);
     
     for(int i = 0 ; i < 6 ; i++)
       {
	 _MCf[i] = new TLorentzVector();
	 std::stringstream name;
	 name << "MCfermion" << i;
	 _tree->Branch(name.str().c_str(),"TLorentzVector",&_MCf[i],16000,0);
	 name << "_PDG";
	 _tree->Branch(name.str().c_str(), &_MCfpdg[i], (name.str()+"/I").c_str());
       }
     
     _tree->Branch("_mass34", &_jetsMass34, "_mass34/F");
     _tree->Branch("_mass35", &_jetsMass35, "_mass35/F");
     _tree->Branch("_mass36", &_jetsMass36, "_mass36/F");
     _tree->Branch("_mass45", &_jetsMass45, "_mass45/F");
     _tree->Branch("_mass46", &_jetsMass46, "_mass46/F");
     _tree->Branch("_mass56", &_jetsMass56, "_mass56/F");
     _tree->Branch("_t1mass", &_t1mass, "_t1mass/F");
     _tree->Branch("_t2mass", &_t2mass, "_t2mass/F");
     _tree->Branch("_t1t2mass", &_t1t2mass, "_t1t2mass/F");
     _tree->Branch("CosTheta", &costheta, "costheta/F");
     _tree->Branch("_W1mass", &_W1mass, "_W1mass/F");
     _tree->Branch("_W2mass", &_W2mass, "_W2mass/F");
     _tree->Branch("_WW1mass", &_WW1mass, "_WW1mass/F");
     _tree->Branch("_WW2mass", &_WW2mass, "_WW2mass/F");
     _tree->Branch("_Wt1mass", &_Wt1mass, "_Wt1mass/F");
     _tree->Branch("_Wt2mass", &_Wt2mass, "_Wt2mass/F");
     _tree->Branch("_Wtx1mass", &_Wtx1mass, "_Wxt1mass/F");
     _tree->Branch("_Wtx2mass", &_Wtx2mass, "_Wtx2mass/F");
     _tree->Branch("_costhetat", &_costhetat ,"_costhetat/F");
     _tree->Branch("_costhetatbar", &_costhetatbar ,"_costhetatbar/F");
     _tree->Branch("_costhetatq", &_costhetatq ,"_costhetatq/F");
     _tree->Branch("_costhetatbarq", &_costhetatbarq ,"_costhetatbarq/F");
     
      _tree->Branch("difff3456", &difff3456 ,"difff3456/F");
      _tree->Branch("difff3546", &difff3546 ,"difff3546/F");
      _tree->Branch("difff3645", &difff3645 ,"difff3645/F");
      
            
      _tree->Branch("difftx1", &difftx1 ,"difftx1/F");
      _tree->Branch("difftx2", &difftx2 ,"difftx2/F");
      
      _tree->Branch("mcbenergx", &mcbenergx ,"mcbenergx/F");
      _tree->Branch("mcbchargex", &mcbchargex ,"mcbchargex/F");
      
}


void ttbar::processRunHeader(LCRunHeader* run)
{ 
  _nRun++;
} 


void ttbar::processEvent(LCEvent * evt) 
{ 
  _nEvt ++;

  cout << "\n------------------------------------------------" << endl;
  cout << "| ******  Run " << _runnumber << " Event Number " << _eventnumber << " ******** |" <<  endl;  
  cout << "------------------------------------------------" << endl;
  ROOTfile->cd();
  
  try {
    bool isttbar=true;
    bool iszww=false;
    qmcb = 0;
    qmcab = 0;    
    isTop=false;
    ismyTop=false;
    isTop1=false;
    isTop2=false;
    B0Meson=true;
    BpMeson=true;
    BmMeson=true;
    int  nmctops=0;
    int  nmccqs=0;
    int  nmcbqs=0;
    int  nmcws=0;
    int  nmcwhs=0;
    int  nmcwls=0;
    mcbenergx=0;
    mcbchargex=0;
    float npartons;
    npartons=0;
    float chargeW1=0;;
    float chargeW2=0;;
    float charge1=0;
    float charge2=0;
    float charge3=0;
    float charge4=0;
    
    nbs=0; 

    TLorentzVector system;
    TLorentzVector MCWf1;
    TLorentzVector MCWf2;
    TLorentzVector MCbf1;
    TLorentzVector MCbf2;
    
    TLorentzVector MCbbf;
    
    TLorentzVector q1;
    TLorentzVector q2;
    TLorentzVector q3;
    TLorentzVector q4;
    TLorentzVector q5;
    TLorentzVector q6;

    int top=0;
    // Get Process name and cross section
    *_Process = evt->getParameters().getStringVal("Process");
    _xsec = evt->getParameters().getFloatVal("CrossSection_fb");
    _eventnumber=evt->getEventNumber();
    _runnumber=evt->getRunNumber();
    LCCollection* mccol = evt->getCollection(_mcCollection);

    for(int i = 0 ; i <mccol->getNumberOfElements(); i++)
      {

	MCParticle *mcpart1 = dynamic_cast<MCParticle*>(mccol->getElementAt(i));
	int PDG=mcpart1->getPDG();
	TVector3 MCPMomentum = TVector3(mcpart1->getMomentum());

        TLorentzVector mctVec(TVector3(mcpart1->getMomentum()),mcpart1->getEnergy());
	MCParticleVec mcpdau = mcpart1->getDaughters();

	if(abs(PDG)==511||abs(PDG)==10511||abs(PDG)==513||abs(PDG)==10513||abs(PDG)==20513||abs(PDG)==515||abs(PDG)==531||abs(PDG)==10531||abs(PDG)==533||abs(PDG)==10533||abs(PDG)==20533 || abs(PDG)==535){
	  TLorentzVector B0Vec(TVector3(mcpart1->getMomentum()),mcpart1->getEnergy());
	  B0MesonVec = B0Vec;
	  B0Meson=false;	
	}
	
	if(abs(PDG)==521||abs(PDG)==10521||abs(PDG)==523||abs(PDG)==10523||abs(PDG)==20523||abs(PDG)==525||abs(PDG)==541||abs(PDG)==10541||abs(PDG)==543||abs(PDG)==10543||abs(PDG)==545 || abs(PDG)==20543 ){
	  if(PDG>0){
	    TLorentzVector BpVec(TVector3(mcpart1->getMomentum()),mcpart1->getEnergy());
	    BpMesonVec =  BpVec ; 
	    BpMeson=false; 
	  }else	if(PDG<0){
	    TLorentzVector BmVec(TVector3(mcpart1->getMomentum()),mcpart1->getEnergy());
	    BmMesonVec= BmVec;
	    BmMeson=false; 
	  }
	  MCParticleVec BPar = mcpart1->getParents();
	}
	if( abs(PDG)<5 ){
	  nqs++;
	  npartons++;
	  
	  if(npartons<5) {
	    //    cout << " ************** Quarks ************ " << endl;
	    //	    cout << PDG << " : " << mcpart1->getCharge()<< " : " <<  mcpart1->getEnergy()<< " : " << npartons <<" : " << MCPMomentum.Mag()<< " : " <<MCPEnergyx <<  endl; 
	    if (PDG==4 ){  //c quarks
	      nmccqs++;
	      ncqs++;
	      TLorentzVector xmccVec(TVector3(mcpart1->getMomentum()),mcpart1->getEnergy()) ;
	      mccVec=xmccVec;
	    } else if (PDG==-4 ){
	      TLorentzVector xmcacVec(TVector3(mcpart1->getMomentum()),mcpart1->getEnergy()) ;
	      mcacVec=xmcacVec;
	      nmccqs++;
	      ncqs++;
	    }
	    if(npartons==1){
	    charge1=mcpart1->getCharge();	 
	    q1=mctVec;
	    }
	    if(npartons==2){
	    charge2=mcpart1->getCharge();	 
	    q2=mctVec;
	    }
	    if(npartons==3){
	    charge3=mcpart1->getCharge();	 
	    q3=mctVec;
	    }
	    if(npartons==4){
	    charge4=mcpart1->getCharge();	 
	    q4=mctVec;
	    }
	  }
	}
	chargeW1=charge1+charge2;	// 1,2 for bbcyyc & bbuyyu 1,4 for yycyyu yyuyyc 
	chargeW2=charge3+charge4;	 // 3,4 for bbxxxx and 2,3 yyxxxx
	if(fabs(chargeW1)>0)  MCWf1 = q1+q2;
	if(fabs(chargeW2)>0)  MCWf2 = q3+q4;
	
	if(abs(PDG)==5){
	  nbs++;
	  if(nbs<3){
	    if(PDG==5) {
	      q5=mctVec;
	      MCbf1=mctVec;
	    }
	    if(PDG==-5) {
	      q6=mctVec;
	      MCbf2=mctVec;
	    }
	  }
	}  
	
	
	MCbbf=q5+q6;
	system=MCbbf+ MCWf1+ MCWf2;	

	//        dift = fabs(MCTopf1.M()-MCTopf2.M());

	if(chargeW1==1){
	  if(MCWf1.E()>0) MCTopf1= MCWf1;
	  if(MCWf2.E()>0) MCTopf2= MCWf2;
	}else if(chargeW2==1){
	  if(MCWf2.E()>0)	  MCTopf1= MCWf2;
	  if(MCWf1.E()>0)	  MCTopf2= MCWf1;

	}

	if(MCbbf.M()<88 || MCbbf.M()>94){
	  MCTopf1b= MCTopf1+q5;
	  MCTopf2b= MCTopf2+q6;	
	}	
	MCTopf1= MCTopf1+q5;
	MCTopf2= MCTopf2+q6;

	if(MCWf1.E()>210||MCWf2.E()>210) isttbar=false;
	if (abs(PDG)==5 ){
	  nmcbqs++;
	  nbqs++;
	  mcbenergx=mcpart1->getEnergy();
	  if(nmcbqs<3)   mcbenergy->Fill(mcpart1->getEnergy());
	  mcbchargex=mcpart1->getCharge();
	  if(nmcbqs<3) NMCBqs->Fill(nmcbqs);
	}


	if (abs(PDG)==24 ){
	  nws++; 
	  nmcws++;
	  NMCWs->Fill(nmcws);
	  MCParticleVec wdau = mcpart1->getDaughters();
	  for (u_int i=0; i< wdau.size();i++){
	    MCParticle *wdaus = wdau[i];
	    int wdaupdg =  wdaus->getPDG() ;
	    if( abs(wdaupdg)<6){
	      nmcwhs++;
	      nwhs++;	 
	    }
	    else if( abs(wdaupdg)==11 || abs(wdaupdg)==13 ||  abs(wdaupdg)==15  ){
	      nmcwls++;
	      nwls++;	 
	    }
	  }	
	}
	if(abs(PDG)==6){
	  top++;
	  nmctops++;
	  if(nmctops<3) NMCTops->Fill(nmctops);
	  TVector3 TopMom(mcpart1->getMomentum()[0],mcpart1->getMomentum()[1],mcpart1->getMomentum()[2]);
	  mccostheta=TopMom.CosTheta();
	  ntops++;
	  if(PDG==6 && nmctops<3 ){
	    h_mccosthetatop->Fill(mccostheta);
	    topangle_mc->Fill(mccostheta);
	  }else if(PDG== -6 &&  nmctops < 3 ){
	    h_mccosthetatopbar->Fill(mccostheta);
	    mccosinustheta=-1*mccostheta;
	    topangle_mc->Fill(mccosinustheta);
	  }
	  MCParticleVec tdau = mcpart1->getDaughters();
	  for (u_int i=0; i< tdau.size();i++){
	    MCParticle *tdaus = tdau[i];
	    int tdauPDG = tdaus->getPDG();
	    TLorentzVector topdauVec(TVector3(tdaus->getMomentum()),tdaus->getEnergy());
	  
	    if (tdauPDG==5 ){
	      //	      cout << " *** MC b jet !!!"  << " : " << tdaus->getEnergy() << endl;
	      qmcb=tdaus->getCharge();
	      mcbVec=topdauVec;
	      /*	      cout << " Angle B0 & b quark " << cos(B0MesonVec.Angle(mcbVec.Vect())) << endl;
	      cout << " Angle Bp & b quark " << cos(BpMesonVec.Angle(mcbVec.Vect())) << endl;
	      cout << " Angle Bm & b quark " << cos(BmMesonVec.Angle(mcbVec.Vect())) << endl;
	      cout<<tdauPDG << " E : "<<topdauVec.E()<< " = "<< mcbVec.E() << " ,Px : " << topdauVec.X() << " =" << mcbVec.X() << endl;
	      */  } else   if (tdauPDG==-5 ){
	      //	      nmcbqs++;
	      //	      nbqs++;
	      qmcab=tdaus->getCharge();
	      mcabVec=topdauVec;
	      /*cout << " *** MC bbar jet !!!" << " : " << tdaus->getEnergy() << endl;
	      cout << " Angle B0 & bbar quark " << cos(B0MesonVec.Angle(mcabVec.Vect())) << endl;
	      cout << " Angle Bp & bbar quark " << cos(BpMesonVec.Angle(mcabVec.Vect())) << endl;
	      cout << " Angle Bm & bbar quark " << cos(BmMesonVec.Angle(mcabVec.Vect())) << endl;
	      cout<<tdauPDG << " E : "<<topdauVec.E()<< " = "<< mcabVec.E() << " ,Px : " << topdauVec.X() << " =" << mcabVec.X() << endl;*/
	    }
	  }
	}
      }

    /*cout << "W1 " << chargeW1 << " W2 " << chargeW2 << endl;	  
	cout << "Mass W 1 " << MCWf1.M() << " Mass W2 " << MCWf2.M() << endl;	  
	cout << " Tops Mass " << MCTopf1.M()<< " : " << MCTopf2.M() << endl;*/
	mctopVec = MCTopf1;
	mcatopVec = MCTopf2;

	NMCWfs->Fill(chargeW1);
	NMCWfs->Fill(chargeW2);
	if(MCWf1.M()>3)	NMCWfmass->Fill(MCWf1.M());
	if(MCWf2.M()>3)	NMCWfmass->Fill(MCWf2.M());
	NMCW12fenergy->Fill(MCWf1.E(),MCWf2.E());   
	NMCTopfmass->Fill(MCTopf1.M());
	NMCTopfmass->Fill(MCTopf2.M());
	NMCTopfmass_sum->Fill(MCTopf1.M()+MCTopf2.M());

	if(chargeW1==1 && MCTopf1.M()<120 )	NMCW1Etheta->Fill(MCWf1.CosTheta(),MCWf1.E());
	if(chargeW2==-1 && MCTopf2.M()<120 )	NMCW2Etheta->Fill(MCWf2.CosTheta(),MCWf2.E());

	if(MCWf1.E()>200 )	NMCW1theta->Fill(MCWf1.CosTheta());
	if(MCWf2.E()>200 )	NMCW2theta->Fill(MCWf2.CosTheta());


	NMCTopftmass->Fill(system.M());
	if(MCbbf.M()>88 && MCbbf.M()<94 ){
	  NMCTopfmass_sum->Fill(MCbbf.M());
	  isttbar=false;
	}

	if(MCWf1.M()>3 && MCWf2.M()>3)	NMCW12fmass->Fill(MCWf1.M(),MCWf2.M() );
	NMCTop12fmass->Fill(MCTopf1.M(),MCTopf2.M() );
	if( MCTopf1.M()>170 &&  MCTopf1.M()<178 &&  MCTopf2.M()>170 &&  MCTopf2.M()<178) NMCTop12fmass_sel->Fill(MCTopf1.M(),MCTopf2.M() );
	if( MCTopf1.M()>170 &&  MCTopf1.M()<178) NMCTop12fmass_sel1->Fill(MCTopf1.M(),MCTopf2.M() );
	if( MCTopf2.M()>170 &&  MCTopf2.M()<178) NMCTop12fmass_sel2->Fill(MCTopf1.M(),MCTopf2.M() );
	if(MCbbf.M()>4)NMCbbfmass->Fill(MCbbf.M());
	
	float tempcostheta = -1 * (MCTopf2.CosTheta());
	if(MCTopf1.M()>120)	topangle_fmc->Fill(MCTopf1.CosTheta());	 
	if(MCTopf2.M()>120)	topangle_fmc->Fill(tempcostheta);	 

	if(MCTopf1.M()>5 || MCTopf1.M()>5 ) ismyTop=true;
	



	if(isttbar){
	  topangle_fmc2->Fill(MCTopf1.CosTheta());	 
	  topangle_fmc2->Fill(tempcostheta);	 
	}
	if(MCTopf1b.M()>120)topangle_fmc3->Fill(MCTopf1b.CosTheta());	 
	if(MCTopf2b.M()>120)topangle_fmc3->Fill(-1*MCTopf2b.CosTheta());	 
 
	/*    cout <<" N MC Tops = " <<  nmctops << endl;
    cout <<" N MC b quarks = " <<  nmcbqs << endl;
    cout <<" N MC c quarks = " <<  nmccqs << endl;
    cout << top << endl;
    cout  << " __________________________________________________________________________________________________ " << endl;*/    
    if (top > 1 )  isTop=true;
    if( (MCbbf.M()>88) && (MCbbf.M()<94) & (!isTop)){
      iszww=true;   
    }
    if(iszww ){
      NMCbbfmasszww->Fill(MCbbf.M());
    }else if ( !iszww && isTop ){
      NMCbbfmasstt->Fill(MCbbf.M());
    }
    

  }catch(DataNotAvailableException &e) {
    cout << "WARNGING can not find collection with name " << _mcCollection << endl;
  }
  
  
  try {
    LCCollection* jetscol = evt->getCollection(_myJetsCollection);
    for(int i = 0 ; i < nJETS ; i++){
      _jetsBtag[i] = 0.;
    }
    
    TLorentzVector jets34(0,0,0,0);
    TLorentzVector jets35(0,0,0,0);
    TLorentzVector jets36(0,0,0,0);
    TLorentzVector jets45(0,0,0,0);
    TLorentzVector jets46(0,0,0,0);
    TLorentzVector jets56(0,0,0,0);
    
    //    float  bpmass=0;
    //    float  bmmass=0;
      
      WMass=80.4;
      TopMass=174.2;
      _t1t2mass = 0;
      _t1mass = 0;
      _t2mass = 0;
      _jetsBtag[0]=0;
      _jetsBtag[1]=0;

      _WW1mass=0;
      _WW2mass=0;
      _Wt1mass=0;
      _Wt2mass=0;
      _Wtx1mass=0;
      _Wtx2mass=0;

      TLorentzVector W1(0,0,0,0);
      TLorentzVector W2(0,0,0,0);

      TLorentzVector WW1(0,0,0,0);
      TLorentzVector WW2(0,0,0,0);
    
      TLorentzVector b1used(0,0,0,0);
      TLorentzVector b2used(0,0,0,0);

      TVector3 b1usedmom(0,0,0);
      TVector3 b2usedmom(0,0,0);
    
      TLorentzVector testtop1(0,0,0,0);    
      TLorentzVector testtop2(0,0,0,0);    
      TLorentzVector testtop3(0,0,0,0);    
      TLorentzVector testtop4(0,0,0,0);    


      TLorentzVector Wtop1(0,0,0,0);    
      TLorentzVector Wtop2(0,0,0,0);    
      TLorentzVector Wtop3(0,0,0,0);    
      TLorentzVector Wtop4(0,0,0,0);    

      TLorentzVector top1(0,0,0,0);
      TLorentzVector top2(0,0,0,0);
      TLorentzVector t1t2(0,0,0,0);

      TLorentzVector topx1(0,0,0,0);
      TLorentzVector topx2(0,0,0,0);
      TLorentzVector t1tx2(0,0,0,0);


      TLorentzVector topxx1(0,0,0,0);
      TLorentzVector topxx2(0,0,0,0);

      TLorentzVector bb(0,0,0,0);

      bbmass=0;

      qptsumoptf1v=-10;
      qptsumoptf2v=-10;
      nbjets=0;
      b1jettag=0;
      b2jettag=0;
      b1jetvq=-10;
      b2jetvq=-10; 
      qptsumb1=-10;
      qptsumb2=-10;
      vtxpartmom=0;
      vtxpartmomsum=0;

      //      float testmass;
      float ptf=-10;
      for(int u=0; u<6;++u){
	damtag[u]=0;
	testpartu2q[u]=-5;
      }

      PIDHandler pidh( jetscol );
      int algo = pidh.getAlgorithmID( "lcfiplus" );
      int ibtag = pidh.getParameterIndex(algo, "BTag");
      int ictag = pidh.getParameterIndex(algo, "CTag");
      int ibctag = pidh.getParameterIndex(algo,"BCTag");
      //      cout << " :::::: List of Pointers in Refined Jets :::::: " << endl;
      bkTag[0]=0;
      bkTag[1]=0;

      for(int i=0; i < jetscol->getNumberOfElements(); i++) {
	ReconstructedParticle *part = 	  dynamic_cast<ReconstructedParticle*>( jetscol->getElementAt( i ) );
	//	cout << part  << " " ;	cout << endl;
	TLorentzVector Jet4vec(TVector3(part->getMomentum()),part->getEnergy());
	jets[i] = Jet4vec;
	const ParticleID &pid = pidh.getParticleID(part, algo);
	jetCharge[i]=part->getCharge();
	jetMass[i]=Jet4vec.M();
	//	testmass = jetMass[i] ;
	bTag[i] = pid.getParameters()[ibtag];
	float qsum2=0;
	float vtxmom;
      	if(bTag[i]>0.5) {
	  if(nbjets==0){
	    b1jettag=bTag[i] ;
	  }else if(nbjets==1){
	    b2jettag=bTag[i] ;	  
	  }
	  nbjets++;
	  //  cout << " Tags today ---> " << b1jettag << " :" << b2jettag << endl;
	  //	  cout << "BJets with btag " << bTag[i] << " & Charge = " << jetCharge[i] << " & Pointer " << part << endl;
	  LCCollection* jetsrelcol = evt->getCollection(_myrelJetsCollection);
	  //	  float Chi2[10];
	  int nvtx=0;
	  for(int k=0; k < jetsrelcol->getNumberOfElements(); k++) {
	    LCRelation * rel =  dynamic_cast<LCRelation*>( jetsrelcol->getElementAt( k ) );
	    ReconstructedParticle *testpartu = dynamic_cast<ReconstructedParticle*> (rel->getFrom());
	    Vertex *testvertex = dynamic_cast <Vertex*> (rel->getTo());
	    ReconstructedParticle *testpartu2 = dynamic_cast<ReconstructedParticle*> (testvertex->getAssociatedParticle());
	    TVector3 tempvtxmom(testpartu2->getMomentum()[0],testpartu2->getMomentum()[1],testpartu2->getMomentum()[2]);
	    vtxmom=tempvtxmom.Mag();
	    //	    Chi2[k]=testvertex->getChi2();
	    if(testpartu==part){
	      nvtx++;
	      if(nvtx>1){
		//		cout << " 2nd Vertex " << endl;
		//		cout << nvtx << " :  First Vertex had  : " <<  qptsum << " / " << ptf << " = " << qptsum/ptf << endl;
		QPtSumVtx1->Fill(qptsum/ptf);	
		if(fabs(qptsum/ptf)==1) PtSumVtx1->Fill(ptf);
	      }
	      h_vtxmom->Fill(vtxmom);	 
	      //	      cout <<"Matching Pointers :" << testpartu << "  = "  << part << endl;
	      //	      cout << " Vertex particle from LC Relation Charge = " << testpartu2->getCharge() << endl; 
	      testpartu2q[k]=testpartu2->getCharge();
	      if(nbjets==1 && nvtx==1){
		//		cout << " First bjet " << bTag[i]<< " : " << b1jettag <<  " Charge = " << testpartu2->getCharge() <<endl;
		b1jetvq=testpartu2->getCharge();
	      }
	      if(nbjets==2 && nvtx==1 ){
		//		cout << " 2nd bjet " << bTag[i]  << " : " << b2jettag <<  " Charge = " << testpartu2->getCharge()  <<   endl;
		b2jetvq=testpartu2->getCharge();
		qptsumb2=qptsum/ptf;
	      }

	      h_Chargesum2->Fill(testpartu2->getCharge());
	      ReconstructedParticleVec testpartu2vec=testpartu2->getParticles();
	      TVector3 tempvtxpartmomsum(0,0,0);
	      qptsum=0;
	      ptf=0;
	      //	      cout << " Q Pt Init " << qptsum << endl;
	      for (u_int n=0; n< testpartu2vec.size();n++){
		ReconstructedParticle *testpartu3 = testpartu2vec[n];
		TLorentzVector tempvtxpart (TVector3(testpartu3->getMomentum()),testpartu3->getEnergy());
		TVector3 tempvtxpartmom(testpartu3->getMomentum()[0],testpartu3->getMomentum()[1],testpartu3->getMomentum()[2]);
		vtxpartmom= tempvtxpartmom.Mag();
		float tempcharge=testpartu3->getCharge();
		float tempqptpro = tempcharge*vtxpartmom;
		if(fabs(tempcharge) <4)	qptsum = qptsum + tempqptpro;
		h_vtxpartmom->Fill(vtxpartmom);	 
		tempvtxpartmomsum+=tempvtxpartmom;
	      }
	      ptf = tempvtxpartmomsum.Mag();	 
	      h_vtxpartmomsum->Fill(tempvtxpartmomsum.Mag());
	      cout << " Q * Pt Test " << qptsum << " / "  << ptf << " = " << qptsum/ptf << endl;
	      qsum2 = qsum2 + (testpartu2->getCharge()) ; // qsum2+=testpartu2->getCharge();
	      TLorentzVector vtxpVec1(TVector3(testpartu2->getMomentum()),testpartu2->getEnergy());
	      vtxpVec=vtxpVec1;
	      
	      float A0 = cos(vtxpVec.Angle(B0MesonVec.Vect()));
	      float Am = cos(vtxpVec.Angle(BmMesonVec.Vect()));
	      float Ap = cos(vtxpVec.Angle(BpMesonVec.Vect()));
	      // v1.Angle(v2.Vect())
	      cout << " A0 : " << A0 << ", Am : " << Am << ", Ap : " << Ap << endl; 
	      float Aqp = cos(vtxpVec.Angle(mcbVec.Vect()));
	      float Aqm = cos(vtxpVec.Angle(mcabVec.Vect()));
	      cout << " Aqp = " << Aqp << ", Aqm = " << Aqm << endl;
	      if(Aqp > Aqm && Aqp> 0.7){
		if(testpartu2q[k]!=-5  &&  nvtx==1  )		h_bplusartonQ->Fill(testpartu2q[k]);
		//		cout << " b+ parton angle with jets " << Aqp << endl;
		//		bpmass=testmass;	   
		//		cout << " Jets used " << testmass<<"= ? "<<jetMass[i]<< endl;
		QPtSumb->Fill(qptsum/(tempvtxpartmomsum.Mag()));
	      }		  
	      if(Aqm > Aqp && Aqm> 0.7  &&  nvtx==1  ) {
		if(testpartu2q[k]!=-5)		h_bminuspartonQ->Fill(testpartu2q[k]);
		//		cout << " b- parton angle with jets " << Aqm << endl;		  
		//		bmmass= testmass;
		//		cout << " Jets used " << testmass << "= ? "<< jetMass[i]<< endl;		   
		QPtSumbbar->Fill(qptsum/(tempvtxpartmomsum.Mag()));	      
	      }
	      
	      if(testpartu2q[k]!=-5)	      h_Bangle->Fill(testpartu2q[k],A0);
	      if(testpartu2q[k]!=-5)	      h_Bpangle->Fill(testpartu2q[k],Ap);
	      if(testpartu2q[k]!=-5)	      h_Bmangle->Fill(testpartu2q[k],Am);
	      


	      if(A0> Ap && A0>Am && A0>0.7){
		if(testpartu2q[k]!=-5 &&  nvtx==1 )		h_BmesonQ->Fill(testpartu2q[k]);		      
		//		cout << " A0 Taken "  << endl;
		if(B0Meson && testpartu2q[k]!=-5) h_BmesonQnoB->Fill(testpartu2q[k]);
	      }
	      if(Ap> A0 && Ap>Am && Ap>0.7){
		if(testpartu2q[k]!=-5  &&  nvtx==1  )		h_BpmesonQ->Fill(testpartu2q[k]);		      
		//		cout << " Ap taken " << endl;
	    	if(BpMeson && testpartu2q[k]!=-5) h_BpmesonQnoB->Fill(testpartu2q[k]);
	      }		
	      if(Am> A0 && Am>Ap && Am>0.7){
		if(testpartu2q[k]!=-5 &&  nvtx==1 )		h_BmmesonQ->Fill(testpartu2q[k]);		      
		//		cout << "Am taken " << endl;
		if(BmMeson && testpartu2q[k]!=-5) h_BmmesonQnoB->Fill(testpartu2q[k]);	      
	      }		
	      if(nvtx>1){
		//		cout << " 2nd Vertex " << endl;
		//		cout << nvtx << " :  Second Vertex has  : " <<  qptsum << " / " << ptf << " = " << qptsum/ptf << endl;
		QPtSumVtx2->Fill(qptsum/ptf);
		if(fabs(qptsum/ptf)==1) PtSumVtx2->Fill(ptf);
		//		qptsumoptf2v=qptsum/ptf;	
	      }
	      if(nbjets==1 && nvtx==1){
		qptsumb1=qptsum/ptf;
	      }else if(nbjets==2 && nvtx==1){
		qptsumb2=qptsum/ptf;	      
	      }
	      
	    }
	  }
	  cout << " Vertices = " << nvtx << ", qptsum/ ptf = " << qptsum/ptf<< endl;     
	  if(nvtx==0){
	    cout << nvtx  << " Vertices " << endl;
	    cout << " Q*Pt Sum for 0 Vertex case : " << qptsum/ptf << endl;	  
	  }


	  NVtx->Fill(nvtx);
	  h_Chargesum3->Fill(qsum2);      //btagvsen1->Fill(b1.E(), _jetsBtag[0]);
	  if(nvtx>0)  QPtSum->Fill(qptsum/ptf);
	  if(nvtx==0)QPtSum0V->Fill(qptsum/ptf);	
	  if(nvtx==1){
	    QPtSum1V->Fill(qptsum/ptf);
	    qptsumoptf1v=qptsum/ptf;
	  }	 
	  if(nvtx==2){
	    QPtSum2V->Fill(qptsum/ptf);
	    qptsumoptf2v=qptsum/ptf;	  
	  }	
	}

	cTag[i] = pid.getParameters()[ictag];
	cTagBBack[i] = pid.getParameters()[ibctag];
	if(bTag[i]>0.5){
	  TrackVec track = part->getTracks();
	}
      }
      /*      cout << " UnSorted Charge & Q*Pt" << endl;      
      cout <<  " 0 : "  << "   "  << b1jettag << " " << b1jetvq  << " : "  << qptsumb1<< endl;
      cout <<  " 1 : "  << "   "  << b2jettag << " " << b2jetvq  << " : "  << qptsumb2<< endl;*/
      if(b1jettag!=b2jettag){
	if(b1jettag>b2jettag){
	  b1jetvq=b1jetvq;   
	  b2jetvq=b2jetvq;   
	}else if(b2jettag>b1jettag){
	  float temp1q =b1jetvq;   
	  float temp2q =b2jetvq; 
	  float temp0 = b1jettag;
	  float temp1 = b2jettag;
	  float temp1qp =qptsumb1;   
	  float temp2qp =qptsumb2; 
	  b1jettag= temp1;
	  b2jettag= temp0;
	  b1jetvq=temp2q;   
	  b2jetvq=temp1q;  
	  qptsumb1=temp2qp; 
	  qptsumb2=temp1qp; 
	}
	/*	cout << " Sorted Charge & Q*Pt " << endl;
      	cout <<  " 0 : "  << b1jettag << " : " << b1jetvq  << " : "  << qptsumb1<< endl;
	cout <<  " 1 : "  << b2jettag << " : " << b2jetvq  << " : "  << qptsumb2<< endl;*/
      }

      int sortPositions[nJETS];
      float tempBtag[nJETS];
      for(int i = 0 ; i < nJETS ; i++) 
	{
	  tempBtag[i] = bTag[i];
	  sortPositions[i] = i;
	}
      for(int i = 0 ; i < nJETS ; i++)
	{
	  int highestBtagPosition = i;
	  for(int j = i ; j < nJETS ; j++)
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
	  _jetsCharge[i]= jetCharge[sortPositions[i]];
	  _jetsMass[i]= jetMass[sortPositions[i]];
	}

      h_btagz->Fill(_jetsBtag[0]);
      h_btagz1->Fill(_jetsBtag[1]);

      //      if(b2jettag>0)

      //      cout <<  " Tag 0 : "  << "   "  << _jetsBtag[0] << ", False Jet Charge   " << _jetsCharge[0] << endl;
      //      cout <<  " Tag 1 : "  << "   "  << _jetsBtag[1] << ", False Jet Charge   " << _jetsCharge[1] << endl;
      TLorentzVector xtemp1 = *_jets[0];
      cout << endl;
      //      cout << " Jets b1 b2 pointers " << & xtemp1 << " : " <<  & *_jets[1] << endl;
      h_btag->Fill(_jetsBtag[0]);
      h_btag1->Fill(_jetsBtag[1]);
      h_btag2->Fill(_jetsBtag[2]);
      TLorentzVector b0test = *_jets[0];
      TVector3 b1 = TVector3(b0test.X(),b0test.Y(),b0test.Z());
      TLorentzVector b1test = *_jets[1];
      TVector3 b2 = TVector3(b1test.X(),b1test.Y(),b1test.Z());
      TLorentzVector b2test = *_jets[2];
      TVector3 b3 = TVector3(b2test.X(),b2test.Y(),b2test.Z());

      //      cout << "Tag value 0 : " <<  _jetsBtag[0] << " Cos theta 0 = " << b1.CosTheta() << endl; 
      //      cout << "Tag value 1 : " <<  _jetsBtag[1] << " Cos theta 1 = " << b2.CosTheta() << endl; 

      btagvstheta1->Fill(b1.CosTheta(), _jetsBtag[0]);
      btagvstheta2->Fill(b2.CosTheta(), _jetsBtag[1]);
      btagvstheta3->Fill(b3.CosTheta(), _jetsBtag[2]);
      btagvsen1->Fill(b0test.E(), _jetsBtag[0]);
      btagvsen2->Fill(b1test.E(), _jetsBtag[1]);
     
      float  sumcharge=-15;
      //   if( _jetsBtag[0]>0.3 && _jetsBtag[1]>0.3 && fabs(_jetsCharge[0])<5 && fabs(_jetsCharge[1])<5  )
      sumcharge=_jetsCharge[0]-_jetsCharge[1];
      if(fabs(sumcharge)<15 )   h_Chargesum->Fill(sumcharge);
      //    cout << " Charge Summ = " << sumcharge << endl;
      jets34=*_jets[2]+*_jets[3];
      jets35=*_jets[2]+*_jets[4];
      jets36=*_jets[2]+*_jets[5];
      
      jets45=*_jets[3]+*_jets[4];
      jets46=*_jets[3]+*_jets[5];
      jets56=*_jets[4]+*_jets[5];
      
      _jetsMass34=jets34.M();
      _jetsMass35=jets35.M();
      _jetsMass36=jets36.M();
      _jetsMass45=jets45.M();
      _jetsMass46=jets46.M();
      _jetsMass56=jets56.M();

      h_mass34->Fill(_jetsMass34);
      h_mass35->Fill(_jetsMass35);
      h_mass36->Fill(_jetsMass36);
      /* 
	 float diff34 = fabs(_jetsMass34 - WMass);
	 float diff35 = fabs(_jetsMass35 - WMass);
	 float diff36 = fabs(_jetsMass36 - WMass);
	 float diff45 = fabs(_jetsMass45 - WMass);
	 float diff46 = fabs(_jetsMass46 - WMass);
	 float diff56 = fabs(_jetsMass56 - WMass);	  
	 if( diff34 < diff35 &&  diff34 < diff36 &&  diff34 < diff45 &&  diff34 < diff46 &&  diff34 < diff56 ){
	 W1=jets34;
	 W2=jets56;
	 }
	 if( diff35 < diff34 &&  diff35 < diff36 &&  diff35 < diff45 &&  diff35 < diff46 &&  diff35 < diff56 ){
	 W1=jets35;
	 W2=jets46;
	 }
	 if( diff36 < diff34 &&  diff36 < diff35 &&  diff36 < diff45 &&  diff36 < diff46 &&  diff36 < diff56 ){
	 W1=jets36;
	 W2=jets45;
	 }
	 if( diff45 < diff34 &&  diff45 < diff34 &&  diff45 < diff36 &&  diff45 < diff46 &&  diff45 < diff56 ){
	 W1=jets45;
	 W2=jets36;
	 }
	 if( diff46 < diff34 &&  diff46 < diff35 &&  diff46 < diff36 &&  diff46 < diff45 &&  diff46 < diff56 ){
	 W1=jets46;
	 W2=jets35;
	 }
	 if( diff56 < diff34 &&  diff56 < diff35 &&  diff56 < diff36 &&  diff56 < diff45 &&  diff56 < diff46 ){
	 W1=jets56;
	 W2=jets34;
	 }
      
      */
	    
      difff3456 = fabs(_jetsMass34-WMass) + fabs(_jetsMass56-WMass);
      difff3546 = fabs(_jetsMass35-WMass) + fabs(_jetsMass46-WMass);
      difff3645 = fabs(_jetsMass36-WMass) + fabs(_jetsMass45-WMass);
      
      h_difff3456->Fill(difff3456) ; 
      h_difff3546->Fill(difff3546) ; 
      h_difff3645->Fill(difff3645) ; 

      if( difff3456 < difff3546 &&  difff3456 < difff3645 &&  difff3456 < 50 ){
	W1=jets34;
	W2=jets56;
      }
      else if( difff3546 < difff3456 &&  difff3546 < difff3645 &&  difff3546 < 50   ){
	W1=jets35;
	W2=jets46;
      }
      else  if( difff3645 < difff3546 &&  difff3645 < difff3456 &&  difff3645 < 50  ){
	W1=jets36;
	W2=jets45;
      }
      

         
      if(W1.M()>0)	_W1mass=W1.M();
      if(W2.M()>0)	_W2mass=W2.M();
      h_w1mass->Fill(_W1mass);
      h_w2mass->Fill(_W2mass);
      recW12fenergy->Fill(W1.E(),W2.E()); 
      //      if( _jetsBtag[0]>0.3 &&  _jetsBtag[1]>0.3 ){
      testtop1=W1+*_jets[0];
      testtop2=W2+*_jets[0];
      testtop3=W1+*_jets[1];
      testtop4=W2+*_jets[1];
      //      }
      b1used=*_jets[0];
      b2used=*_jets[1];
 
      b1en->Fill(b1used.E());
      b2en->Fill(b2used.E());

      if(W1.E()>0) w1en->Fill(W1.E());
      if(W1.E()>0) w2en->Fill(W2.E());

      bb=b1used+b2used;
      bbmass=bb.M();

      h_bbmass->Fill(bbmass);   
      if(isTop) h_bbmass_tag->Fill(bbmass);       
      else h_bbmass_notag->Fill(bbmass);       
//      float bplustheta = cos(b1used.Angle(mcbVec.Vect()));
      //      float bminustheta = cos(b1used.Angle(mcabVec.Vect()));


      float bplustheta = cos(b1used.Angle(BpMesonVec.Vect()));
      float bminustheta = cos(b1used.Angle(BmMesonVec.Vect()));

      if( bplustheta>bminustheta && bplustheta>0.8 && fabs(_jetsCharge[0])<3)	h_bplusartonQ2->Fill(_jetsCharge[0]); 
      else if(bminustheta>bplustheta && bminustheta>0.8  && fabs(_jetsCharge[0])<3 )	h_bminuspartonQ2->Fill(_jetsCharge[0]); 

      TVector3 b1usedmomx(b1used.X(),b1used.Y(),b1used.Z());
      //     b1usedmom((b1used.X(),b1used.Y(),b1used.Z()));
      b1usedmom=b1usedmomx;      
      //      float b2usedmass = b2used.M();
      //      float b1usedmass = b1used.M();

      float bmcbrec1 = -999;
      float cmcbrec1 = -999 ;


      float bmcbrec2 = -999;
      float cmcbrec2 = -999;
      
      float varc=10;


	if(fabs(b1jetvq)<5 && fabs(b2jetvq)<5) varc = b1jetvq  -  b2jetvq ;


	if(varc>0) {
	  cmcbrec1 = cos(mcacVec.Angle(b1used.Vect()));
	  bmcbrec1 = cos(mcabVec.Angle(b1used.Vect()));
	  bmcbrec2 = cos(mcbVec.Angle(b2used.Vect()));
	  cmcbrec2 = cos(mccVec.Angle(b2used.Vect()));
	  
	}


	if(varc<0) {
	  cmcbrec1 = cos(mccVec.Angle(b1used.Vect()));
	  bmcbrec1 = cos(mcbVec.Angle(b1used.Vect()));

	  bmcbrec2 = cos(mcabVec.Angle(b2used.Vect()));
	  cmcbrec2 = cos(mcacVec.Angle(b2used.Vect()));
	}
	
	if(cmcbrec1!=-999 )	coscmcbrec->Fill(cmcbrec1);
	if(cmcbrec2!=-999 )	coscmcbrec->Fill(cmcbrec2);
	cosb1mcb2mc->Fill(cos(mcbVec.Angle(mcabVec.Vect()))); 
	cosb1mccmc->Fill(cos(mcbVec.Angle(mccVec.Vect()))); 
	cosb2mccmc->Fill(cos(mcabVec.Angle(mccVec.Vect()))); 

	if(   cos(mcbVec.Angle(mcabVec.Vect())) >=0  ){
	cosb1mccmccut->Fill( cos(mcbVec.Angle(mccVec.Vect()))); 
	cosb2mccmccut->Fill( cos(mcabVec.Angle(mccVec.Vect()))); 
	}

	if(  cos(mcbVec.Angle(mccVec.Vect())) >= 0.7 || cos(mcabVec.Angle(mccVec.Vect())) >= 0.7  )	cosb1mcb2mccut->Fill(cos(mcbVec.Angle(mcabVec.Vect())));

	cosb1mcb2mcb1mccmc->Fill( cos(mcbVec.Angle(mcabVec.Vect())) , cos(mcbVec.Angle(mccVec.Vect()))  );
	cosb1mcb2mcb2mccmc->Fill( cos(mcbVec.Angle(mcabVec.Vect())) , cos(mcabVec.Angle(mccVec.Vect())) );

       	if(cmcbrec1!=-999 && bmcbrec1!=-999 ) {
	  coscmcbrecvsbmc1->Fill(cmcbrec1, bmcbrec1);
	  cout << " ---------------------------------------------------------------------- " << endl;
	  cout << " c MC and b rec " << cmcbrec1 << " b MC and b rec " << bmcbrec1 << endl;
	  cout << " ---------------------------------------------------------------------- " << endl;
	if( cmcbrec1 > bmcbrec1 ) 	coscmcbrecvsbmc21->Fill(cmcbrec1, bmcbrec1);
	if( cmcbrec1 > bmcbrec1  && cmcbrec1 > 0 )	coscmcbrecvsbmc31->Fill(cmcbrec1, bmcbrec1);
	if( cmcbrec1 > bmcbrec1  && cmcbrec1 > 0  && bmcbrec1 <= 0 ){
	coscmcbrecvsbmc41->Fill(cmcbrec1, bmcbrec1);
	btagsusp1->Fill( _jetsBtag[0]);
	}	
	cosbreccmcbtag1-> Fill( _jetsBtag[1], cmcbrec1);
	cosbreccmcctag1-> Fill( _jetsCtag[1], cmcbrec1);
	if( _jetsBtag[0] > 0.3) {
	  coscmcbrecvsbmcc1->Fill(cmcbrec1, bmcbrec1);
	  if( cmcbrec1 > bmcbrec1 ) 	coscmcbrecvsbmc2c1->Fill(cmcbrec1, bmcbrec1);
	  if( cmcbrec1 > bmcbrec1  && cmcbrec1 > 0 )	coscmcbrecvsbmc3c1->Fill(cmcbrec1, bmcbrec1);
	  if( cmcbrec1 > bmcbrec1  && cmcbrec1 > 0  && bmcbrec1 <= 0 )  coscmcbrecvsbmc4c1->Fill(cmcbrec1, bmcbrec1);
	  cosbreccmcbtagc1->Fill( _jetsBtag[0],cmcbrec1)  ; 
	  cosbreccmcctagc1->Fill( _jetsCtag[0],cmcbrec1)  ;
	}      
	
	cosbreccmcbtag1->Fill( _jetsBtag[0],cmcbrec1)  ;
	cosbreccmcctag1 ->Fill( _jetsCtag[0],cmcbrec1)  ;
	}
	
       	if(cmcbrec2!=-999 && bmcbrec2!=-999 ) {
	  coscmcbrecvsbmc2->Fill(cmcbrec2, bmcbrec2);
	  if( cmcbrec2 > bmcbrec2 ) 	coscmcbrecvsbmc22->Fill(cmcbrec2, bmcbrec2);
	  if( cmcbrec2 > bmcbrec2  && cmcbrec2 > 0 )	coscmcbrecvsbmc32->Fill(cmcbrec2, bmcbrec2);
	if( cmcbrec2 > bmcbrec2  && cmcbrec2 > 0  && bmcbrec2 <= 0 ){
	coscmcbrecvsbmc42->Fill(cmcbrec2, bmcbrec2);
	btagsusp2->Fill( _jetsBtag[1]);
	}	
	cosbreccmcbtag2-> Fill( _jetsBtag[1], cmcbrec2);
	cosbreccmcctag2-> Fill( _jetsCtag[1], cmcbrec2);
	if( _jetsBtag[1] > 0.3) {
	  coscmcbrecvsbmcc2->Fill(cmcbrec2, bmcbrec2);
	  if( cmcbrec2 > bmcbrec2 ) 	coscmcbrecvsbmc2c2->Fill(cmcbrec2, bmcbrec2);
	  if( cmcbrec2 > bmcbrec2  && cmcbrec2 > 0 )	coscmcbrecvsbmc3c2->Fill(cmcbrec2, bmcbrec2);
	  if( cmcbrec2 > bmcbrec2  && cmcbrec2 > 0  && bmcbrec2 <= 0 )  coscmcbrecvsbmc4c2->Fill(cmcbrec2, bmcbrec2);
	  cosbreccmcbtagc2->Fill( _jetsBtag[1],cmcbrec2)  ; 
	  cosbreccmcctagc2->Fill( _jetsCtag[1],cmcbrec2)  ;
	}      
	
	cosbreccmcbtag2->Fill( _jetsBtag[1],cmcbrec2)  ;
	cosbreccmcctag2 ->Fill( _jetsCtag[1],cmcbrec2)  ;
	}

	//	coscmcbrec->Fill(cos(mcacVec.Angle(b2used.Vect())));


      testop1Mass = testtop1.M();   
      testop2Mass = testtop2.M();   
      testop3Mass = testtop3.M();   
      testop4Mass = testtop4.M();   
     
      Double_t bpstar1 = testtop1.Gamma()*_jets[1]->Energy()*(1 - testtop1.Beta()*TMath::Cos(testtop1.Angle(-_jets[1]->Vect())));
      Double_t bpstar2 = testtop2.Gamma()*_jets[1]->Energy()*(1 - testtop2.Beta()*TMath::Cos(testtop2.Angle(-_jets[1]->Vect())));
      Double_t bpstar3 = testtop3.Gamma()*_jets[0]->Energy()*(1 - testtop3.Beta()*TMath::Cos(testtop3.Angle(-_jets[0]->Vect())));
      Double_t bpstar4 = testtop4.Gamma()*_jets[0]->Energy()*(1 - testtop4.Beta()*TMath::Cos(testtop4.Angle(-_jets[0]->Vect())));      
      
      float difft1=pow((testop1Mass-TopMass)/6.3,2) + pow((testtop1.E()-250.0)/8.0,2)+pow((bpstar1-68)/10.0,2);
      float difft2=pow((testop2Mass-TopMass)/6.3,2) + pow((testtop2.E()-250.0)/8.0,2)+pow((bpstar2-68)/10.0,2);;
      float difft3=pow((testop3Mass-TopMass)/6.3,2) + pow((testtop3.E()-250.0)/8.0,2)+pow((bpstar3-68)/10.0,2);;
      float difft4=pow((testop4Mass-TopMass)/6.3,2) + pow((testtop4.E()-250.0)/8.0,2)+pow((bpstar4-68)/10.0,2);;
           
      cout << difft1 << ":" <<  difft4 << ":" << difft3 << ":" << difft4  << endl;

      cout << " ------------------  Tag Test ------------------------- " << endl;
      cout <<  _jetsBtag[0] << " = " << b1jettag << " : " <<  _jetsBtag[1] << " = " << b2jettag << endl;
      Double_t bp1=0,bp2=0;
      if(difft1<difft2){
	top1=testtop1;
	difftx1=difft1;
	bp1=bpstar1;
      } else if(difft2<difft1){
	top1=testtop2;
	difftx1=difft2;
	bp1=bpstar2;
      }
      if(difft3<difft4){
	top2=testtop3;
	difftx2=difft3;
	bp2=bpstar3;
      }else if(difft4<difft3){
	top2=testtop4;
	difftx2=difft4;
	bp2=bpstar4;
      }

      if(top1.M()>5 )      mtforchi->Fill(top1.M()) ;
      if(top1.E()>5 )      etforchi->Fill(top1.E()) ;
      if(bp1>5 )      pbforchi->Fill(bp1) ;




      float	testqsumx = -20;
      if(fabs(b1jetvq)<5 && fabs(b2jetvq)<5)  testqsumx= b1jetvq  -  b2jetvq ;
      if(fabs(testqsumx)<20)QSum->Fill(testqsumx);
      if( fabs(testqsumx)<3 ){
	QSum3->Fill(testqsumx);	
      }
      bool good1=false;
      bool good2=false;
      //      bool rightq1=true;
      //      bool rightq2=true;
      cout << " for setting good1 and good2 " << endl;
      cout << " top1.M() 140-210  " << top1.M() << ", Chi1 " <<  difftx1  << endl;
      cout << " top2.M() 140-210  " << top2.M() << ", Chi2 " <<  difftx2  << endl;

      if( top1.M()>140 && top1.M()<210 && difftx1 <20 ) good1=true;
      if( top2.M()>140 && top2.M()<210 && difftx2 <40 ) good2=true;

      b1mom->Fill(b1usedmom.Mag());

      float	testqptsum = -20;
      testqptsum= qptsumb1  - qptsumb2; 
      _t1mass=top1.M();
      _t2mass=top2.M();

      NEvts->Fill(1);
	Topmasseff->Fill(_t1mass);
	Topmasseff->Fill(_t2mass);
	if(_t1mass>10 && _t2mass>10) {
	Topmasseff1->Fill(_t1mass);
	Topmasseff1->Fill(_t2mass);
	if(cos(mcbVec.Angle(b1used.Vect()))>=0.75){
	if(MCTopf1.CosTheta()>0 && top1.CosTheta()>0)	Topmasseff1t->Fill(_t1mass);
	if(MCTopf2.CosTheta()>0 && top2.CosTheta()>0)	Topmasseff1t->Fill(_t2mass);	  
	if(MCTopf1.CosTheta()<0 && top1.CosTheta()<0)	Topmasseff1t->Fill(_t1mass);
	if(MCTopf2.CosTheta()<0 && top2.CosTheta()<0)	Topmasseff1t->Fill(_t2mass);	  
	} else if (cos(mcbVec.Angle(b1used.Vect()))>=0.9){
	  if(MCTopf2.CosTheta()>0 && top1.CosTheta()>0)	Topmasseff1t->Fill(_t1mass);
	  if(MCTopf1.CosTheta()>0 && top2.CosTheta()>0)	Topmasseff1t->Fill(_t2mass);	  
	  if(MCTopf2.CosTheta()<0 && top1.CosTheta()<0)	Topmasseff1t->Fill(_t1mass);
	  if(MCTopf1.CosTheta()<0 && top2.CosTheta()<0)	Topmasseff1t->Fill(_t2mass);	  
	  
	}

	//	if( _jetsBtag[0]>0.3 &&  _jetsBtag[1]>0.3){
	if( _jetsBtag[0]>0.8 &&  _jetsBtag[1]>0.3){
	  Topmasseff2->Fill(_t1mass);
	  Topmasseff2->Fill(_t2mass);
	  if(MCTopf1.CosTheta()>0 && top1.CosTheta()>0)	Topmasseff2t->Fill(_t1mass);
	  if(MCTopf2.CosTheta()>0 && top2.CosTheta()>0)	Topmasseff2t->Fill(_t2mass);	  
	  if(MCTopf1.CosTheta()<0 && top1.CosTheta()<0)	Topmasseff2t->Fill(_t1mass);
	  if(MCTopf2.CosTheta()<0 && top2.CosTheta()<0)	Topmasseff2t->Fill(_t2mass);	  
	  if (top1.M()> 140 && top1.M()<210 && top2.M()> 140 && top2.M()<210){
	    Topmasseff3->Fill(_t1mass);
	    Topmasseff3->Fill(_t2mass);
	    if(MCTopf1.CosTheta()>0 && top1.CosTheta()>0)	Topmasseff3t->Fill(_t1mass);
	    if(MCTopf2.CosTheta()>0 && top2.CosTheta()>0)	Topmasseff3t->Fill(_t2mass);	  
	    if(MCTopf1.CosTheta()<0 && top1.CosTheta()<0)	Topmasseff3t->Fill(_t1mass);
	    if(MCTopf2.CosTheta()<0 && top2.CosTheta()<0)	Topmasseff3t->Fill(_t2mass);	  
	    //	    if(difftx1 <20 && difftx2 <40 ){
	    if(difftx1 <30 && difftx2 <30 ){
	      Topmasseff4->Fill(_t1mass);
	      Topmasseff4->Fill(_t2mass);
	      if(MCTopf1.CosTheta()>0 && top1.CosTheta()>0)	Topmasseff4t->Fill(_t1mass);
	      if(MCTopf2.CosTheta()>0 && top2.CosTheta()>0)	Topmasseff4t->Fill(_t2mass);	  
	      if(MCTopf1.CosTheta()<0 && top1.CosTheta()<0)	Topmasseff4t->Fill(_t1mass);
	      if(MCTopf2.CosTheta()<0 && top2.CosTheta()<0)	Topmasseff4t->Fill(_t2mass);	  
	      if(testqsumx!=0){
		Topmasseff5->Fill(_t1mass);
		Topmasseff5->Fill(_t2mass);
		if(cos(mcbVec.Angle(b1used.Vect()))>=0.75){
		  if(MCTopf1.CosTheta()>0 && top1.CosTheta()>0)	Topmasseff5t->Fill(_t1mass);
		  if(MCTopf2.CosTheta()>0 && top2.CosTheta()>0)	Topmasseff5t->Fill(_t2mass);	  
		  if(MCTopf1.CosTheta()<0 && top1.CosTheta()<0)	Topmasseff5t->Fill(_t1mass);
		  if(MCTopf2.CosTheta()<0 && top2.CosTheta()<0)	Topmasseff5t->Fill(_t2mass);	  
		}else if (cos(mcbVec.Angle(b1used.Vect()))>=0.9){
		  if(MCTopf2.CosTheta()>0 && top1.CosTheta()>0)	Topmasseff1t->Fill(_t1mass);
		  if(MCTopf1.CosTheta()>0 && top2.CosTheta()>0)	Topmasseff1t->Fill(_t2mass);	  
		  if(MCTopf2.CosTheta()<0 && top1.CosTheta()<0)	Topmasseff1t->Fill(_t1mass);
		  if(MCTopf1.CosTheta()<0 && top2.CosTheta()<0)	Topmasseff1t->Fill(_t2mass);	  
		}	      
	      }
	    }
	  }
	}
	} //  if(_t1mass>10 && _t2mass>10)
      
      if(_t1mass>10) h_top1mass->Fill(_t1mass);
      if(_t2mass>10) h_top2mass->Fill(_t2mass);
      if(_t1mass>10 && _t2mass>10 ) h_top1top2mass->Fill(_t1mass, _t2mass);
      _t1t2mass = _t1mass+_t2mass;
      if(_t1mass>0)   topmasschia->Fill(_t1mass,difftx1);
      if(_t2mass>0)   topmasschib->Fill(_t2mass,difftx2);
      if(_t1mass>0)   topmasschiaE->Fill(top1.E(),difftx1);
      if(_t2mass>0)   topmasschibE->Fill(top2.E(),difftx2);
      if(_t1mass>0)   topmasschiaP->Fill(bp1,difftx1);
      if(_t2mass>0)   topmasschibP->Fill(bp2,difftx2);
      if(_t1mass>10) top1momvsmass->Fill(_t1mass, b1usedmom.Mag());
      cout << " Checking good1 =  " << good1 << ", and good2 = " << good2  << endl;
      if(good1 && good2 && _jetsBtag[0]>0.3 &&  _jetsBtag[1]>0.3  ){
	if(b1jettag>0.5 && b2jettag>0.5){
	  if(cos(mcbVec.Angle(b1used.Vect()))>0.999 ){
	    _costhetatq = top1.CosTheta();	
	    if( fabs(b1jetvq)<10 && fabs(b2jetvq)<10 && fabs(testqsumx)<20 ){
	      h_costhetatop->Fill(_costhetatq);
	      gud_bjbmcangle->Fill(cos(mcbVec.Angle(b1used.Vect()))); 
	      gudbjmom->Fill(b1usedmom.Mag());
	      gudbjtag->Fill(b1jettag);
	      if(fabs(testqsumx)<4)QSum_mcbp->Fill(testqsumx);
	      if(fabs(qptsumoptf1v)<10)QPtSum1V_S->Fill(qptsumoptf1v);
	      if(fabs(qptsumoptf2v)<10)QPtSum2V_S->Fill(qptsumoptf2v);
	      h_costhetatopbar->Fill(top2.CosTheta());	  
	      if(fabs(testqptsum)<4) QPtSumb_C->Fill(qptsumb1);  
	      if(fabs(b1jetvq)<10)  qmcb_b->Fill(qmcb,b1jetvq);
	      if(b1jetvq>0){
		//		Wrong_Q->Fill(top1.CosTheta());  
		//		rightq1=false;
	      }
	    }	   
	    if(fabs(b1jetvq)<10)  b1jetvq_b->Fill(b1jetvq);
	    if(fabs(b2jetvq)<10)   b2jetvq_b->Fill(b2jetvq);
	    QPtSum_mcbp->Fill(qptsum/ptf);
	  }
	  if(cos(mcbVec.Angle(b2used.Vect()))>0.999 && b2jetvq>0 && fabs(b2jetvq)<10 && fabs(testqsumx)<20){
	    //	    Wrong_Q->Fill(-1*top2.CosTheta());  
	    //	    rightq2=false;
	  }	  
	  if(cos(mcabVec.Angle(b1used.Vect()))>0.999 ){
	    _costhetatbarq = top1.CosTheta();	
	    if( fabs(b1jetvq)<10 && fabs(b2jetvq)<10 && fabs(testqsumx<20)){  
	      h_costhetatopbar->Fill(_costhetatbarq);	  
	      h_costhetatop->Fill(top2.CosTheta());	  
	      if(fabs(testqsumx)<4) QSum_mcabp->Fill(testqsumx);
	      QPtSum_mcabp->Fill(qptsum/ptf);
	      if(fabs(b1jetvq)<10)   b1jetvq_ab->Fill(b1jetvq);
	      if(fabs(b2jetvq)<10)  b2jetvq_ab->Fill(b2jetvq);
	      if(fabs(testqptsum)<4) QPtSumbbar_C->Fill(qptsumb1);  
	      if(fabs(b1jetvq)<10)  qmcab_ab->Fill(qmcb,b1jetvq);
	      if(b1jetvq<0){
		//      Wrong_Q->Fill(-1*top1.CosTheta());  
		//		rightq1=false;	   
	      }	    
	    }  
	  }
	  if(cos(mcabVec.Angle(b2used.Vect()))>0.999  && b2jetvq<0 && fabs(b2jetvq)<10 && fabs(testqsumx)<20){
	    //	    Wrong_Q->Fill(top2.CosTheta());  
	    //	    rightq2=false;	  
	  }

	  cout << " Test Q SUM X " << testqsumx << endl; 
	  if(testqsumx>0 ) {
	    //	  if(testqsumx>0) {
	    cout << " fabs(b1jetvq) 10 : fabs(b2jetvq) 10 : fabs(testqsumx) 7 " << endl;
	    cout << fabs(b1jetvq)<< " : " <<  fabs(b2jetvq) << " : " <<  fabs(testqsumx) << endl; 
	    if( fabs(b1jetvq)<10 && fabs(b2jetvq)<10 &&  fabs(testqsumx)<7 ){
	      topbarquarkangleq->Fill(top1.CosTheta());
	      topquarkangleq->Fill(top2.CosTheta());
	      bad_bjbmcangle->Fill(cos(mcbVec.Angle(b1used.Vect()))); 	  
	      badbjmom->Fill(b1usedmom.Mag());
	      badbjtag->Fill(b1jettag);
	      float cosinusthetaa =-1*(top1.CosTheta());
	      topangleq->Fill(top2.CosTheta());	
	      topangleq->Fill(cosinusthetaa);	
	      Topanglecomp1->Fill(top2.CosTheta());
	      Topanglecomp1->Fill(cosinusthetaa);
	     
	      if( b1jettag>0.5 && b2jettag>0.5) {
		Topanglecomp2->Fill(top2.CosTheta());
		Topanglecomp2->Fill(cosinusthetaa);
	      }
	      if( b1jettag>0.7 && b2jettag>0.7) {
		Topanglecomp3->Fill(top2.CosTheta());
		Topanglecomp3->Fill(cosinusthetaa);
	      }
	      if( b1jettag>0.9 && b2jettag>0.8) {
		Topanglecomp4->Fill(top2.CosTheta());
		Topanglecomp4->Fill(cosinusthetaa);
	      }	
	      //	      float diffang =  cos( b1used.CosTheta()) - cos (mcabVec.CosTheta());
	      cosbmcbrec->Fill(cos(mcabVec.Angle(b1used.Vect())));
	      if( cos(mcabVec.Angle(b1used.Vect()))>=0.8){
		b1momsel->Fill(b1usedmom.Mag()); // selected b jets momentum 
		topangleq_c->Fill(top2.CosTheta());
		topangleq_c->Fill(cosinusthetaa);
		Topanglecomp1c->Fill(top2.CosTheta());
		Topanglecomp1c->Fill(cosinusthetaa);	      
		if( b1jettag>0.5 && b2jettag>0.5) {
		  Topanglecomp2c->Fill(top2.CosTheta());
		  Topanglecomp2c->Fill(cosinusthetaa);
		}
		if( b1jettag>0.7 && b2jettag>0.7) {
		  Topanglecomp3c->Fill(top2.CosTheta());
		  Topanglecomp3c->Fill(cosinusthetaa);
		}
		if( b1jettag>0.9 && b2jettag>0.8) {
		  Topanglecomp4c->Fill(top2.CosTheta());
		  Topanglecomp4c->Fill(cosinusthetaa);
		}	    		

	      }	    
	    }		
	  }
	  
	  if(testqsumx<0){
	    if( fabs(b1jetvq)<10 && fabs(b2jetvq)<10  && fabs(testqsumx)<7 ){
	      topquarkangleq->Fill(top1.CosTheta());
	      topbarquarkangleq->Fill(top2.CosTheta());
	      float cosinusthetab =-1* (top2.CosTheta());
	      //	    if(cos(top1.Angle(top2.Vect()))>0.9)
	      topangleq->Fill(cosinusthetab);	
	      topangleq->Fill(top1.CosTheta());	
	      Topanglecomp1->Fill(top1.CosTheta());
	      Topanglecomp1->Fill(cosinusthetab);
	      if( b1jettag>0.5 && b2jettag>0.5) {
		Topanglecomp2->Fill(top1.CosTheta());
		Topanglecomp2->Fill(cosinusthetab);
	      }
	      if( b1jettag>0.7 && b2jettag>0.7) {
		Topanglecomp3->Fill(top1.CosTheta());
		Topanglecomp3->Fill(cosinusthetab);
	      }
	      if( b1jettag>0.9 && b2jettag>0.8) {
		Topanglecomp4->Fill(top1.CosTheta());
		Topanglecomp4->Fill(cosinusthetab);
	      }	    
	      cosbmcbrec->Fill(cos(mcbVec.Angle(b1used.Vect())));
	      //	      coscmcbrec->Fill(cos(mcacVec.Angle(b1used.Vect())));
	      if( cos(mcbVec.Angle(b1used.Vect()))>=0.8){
		topangleq_c->Fill(top1.CosTheta());
		topangleq_c->Fill(cosinusthetab);
		Topanglecomp1c->Fill(top1.CosTheta());
		Topanglecomp1c->Fill(cosinusthetab);
		b1momsel->Fill(b1usedmom.Mag());	 
		if( b1jettag>0.5 && b2jettag>0.5) {
		  Topanglecomp2c->Fill(top1.CosTheta());
		  Topanglecomp2c->Fill(cosinusthetab);
		}
		if( b1jettag>0.7 && b2jettag>0.7) {
		  Topanglecomp3c->Fill(top1.CosTheta());
		  Topanglecomp3c->Fill(cosinusthetab);
		}
		if( b1jettag>0.9 && b2jettag>0.8) {
		  Topanglecomp4c->Fill(top1.CosTheta());
		  Topanglecomp4c->Fill(cosinusthetab);
		}	    
	      }	    
	    }
	  }
	  
	  
	  //	  float	qbsum = b1jetvq + qmcb ;
	  //	  float	qabsum = b1jetvq + qmcab ;
	  //	  cout << "Charge vs MC Charge " << b1jetvq << " : " << qmcb  << " S1 : " << qbsum << " S2 : " << qabsum << endl;
	  if( b1jettag>0.5 && b2jettag>0.5 && b1jettag==_jetsBtag[0]  && b2jettag==_jetsBtag[1] ){
	    if(testqsumx<0 ){
	      float cosinusthetax =-1*(top2.CosTheta());
	      topangle_bc->Fill(top1.CosTheta());
	      topangle_bc->Fill(cosinusthetax);
	      if(testqsumx==-1 ){
		topangle_c1bc->Fill(top1.CosTheta());
		topangle_c1bc->Fill(cosinusthetax);
		b1jetvq_c1bc->Fill(b1jetvq);
		b2jetvq_c1bc->Fill(b2jetvq);
	      }	    
	      if(testqsumx<-1){
		topangle_c2bc->Fill(top1.CosTheta());
		topangle_c2bc->Fill(cosinusthetax);
		b1jetvq_c2bc->Fill(b1jetvq);
		b2jetvq_c2bc->Fill(b2jetvq);
	      }	    
	    }

	    if(testqsumx>0 ){
	      float cosinus = (-1)*(top1.CosTheta());		
	      topangle_bc->Fill(top2.CosTheta());
	      topangle_bc->Fill(cosinus);
	      if(testqsumx==1 ){
		topangle_c1bc->Fill(top2.CosTheta());
		topangle_c1bc->Fill(cosinus);
		b1jetvq_c1bc->Fill(b1jetvq);
		b2jetvq_c1bc->Fill(b2jetvq);
	      }	    
	      if(testqsumx>1){
		topangle_c2bc->Fill(top2.CosTheta());
		topangle_c2bc->Fill(cosinus);
		b1jetvq_c2bc->Fill(b1jetvq);
		b2jetvq_c2bc->Fill(b2jetvq);
	      }	    
	    }

	    float cosangleval=0.9;

	    if(testqsumx<0 ){
	      topquarkangle->Fill(top1.CosTheta());
	      topbarquarkangle->Fill(top2.CosTheta());
	      float cosinustheta =-1*(top2.CosTheta());
	      float difftmctrec1 ;
	      float difftmctrec2 ;
	      difftmctrec1  = top1.CosTheta() -  mctopVec.CosTheta();
	      difftmctrec2  = top2.CosTheta() -  mcatopVec.CosTheta();
	      costmctrecdiffbc->Fill(difftmctrec1);
	      costmctrecdiffbc->Fill(difftmctrec2);
	      if(cos(mcbVec.Angle(b1used.Vect()))>= cosangleval){
		costmctrecdiff->Fill(difftmctrec1);
		costmctrecdiff->Fill(difftmctrec2);
		topangle->Fill(top1.CosTheta());
		topangle->Fill(cosinustheta);

		h_top1top2mass_sel->Fill(_t1mass, _t2mass);
		h_bbmass_sel->Fill(bbmass);   
		float diff = mccostheta - (top1.CosTheta());	  
		float diff2 =  mccosinustheta - cosinustheta;
		topangle_diff->Fill(diff);
		topangle_diff->Fill(diff2);
		if(testqsumx==-1){
		  topangle_c1->Fill(top1.CosTheta());
		  topangle_c1->Fill(cosinustheta);
		  b1jetvq_c1->Fill(b1jetvq);
		  b2jetvq_c1->Fill(b2jetvq);	
		}	    
		if(testqsumx<-1){
		  topangle_c2->Fill(top1.CosTheta());
		  topangle_c2->Fill(cosinustheta);
		  b1jetvq_c2->Fill(b1jetvq);
		  b2jetvq_c2->Fill(b2jetvq);
		}	    
	      }
	      if(cos(mcbVec.Angle(b1used.Vect()))<cosangleval){
		Wrong_Q->Fill(top1.CosTheta());  
	        Wrong_Q->Fill(cosinustheta);  
	      }
	    }
	    if(testqsumx>0){
	      topbarquarkangle->Fill(top1.CosTheta());
	      topquarkangle->Fill(top2.CosTheta());
	      float cosinustheta1 =-1*(top1.CosTheta());
	      float difftmctrecx1 ;
	      float difftmctrecx2 ;
 
	      difftmctrecx1  = top2.CosTheta() -  mctopVec.CosTheta();
	      difftmctrecx2  = top1.CosTheta() -  mcatopVec.CosTheta();
	      costmctrecdiffbc->Fill(difftmctrecx1);
	      costmctrecdiffbc->Fill(difftmctrecx2);
	      if(cos(mcabVec.Angle(b1used.Vect()))>=cosangleval){
		costmctrecdiff->Fill(difftmctrecx1);
		costmctrecdiff->Fill(difftmctrecx2);
		topangle->Fill(top2.CosTheta());
		topangle->Fill(cosinustheta1);
		h_top1top2mass_sel->Fill(_t1mass, _t2mass);	
		h_bbmass_sel->Fill(bbmass);   
		float difff = mccostheta - (top2.CosTheta());	  
		float difff2 =  mccosinustheta - cosinustheta1;
		topangle_diff->Fill(difff);
		topangle_diff->Fill(difff2);
		if(testqsumx==1){
		  topangle_c1->Fill(top2.CosTheta());
		  topangle_c1->Fill(cosinustheta1);
		  b1jetvq_c1->Fill(b1jetvq);
		  b2jetvq_c1->Fill(b2jetvq);
		}	    
		if(testqsumx>1){
		  topangle_c2->Fill(top2.CosTheta());
		  topangle_c2->Fill(cosinustheta1);
		  b1jetvq_c2->Fill(b1jetvq);
		  b2jetvq_c2->Fill(b2jetvq);	
		}	    	      
	      }
	      if(cos(mcabVec.Angle(b1used.Vect()))<cosangleval){
		Wrong_Q->Fill(top2.CosTheta());  
		Wrong_Q->Fill(cosinustheta1);  
	      }
	    }
	  } //   if( b1jettag>0 && b2jettag>0){ before testqsumx<0 
	  
	  cost1t2->Fill(cos(top1.Angle(top2.Vect())));
	}
	
	if(testqsumx==0 && b1jettag>0.3 && fabs(testqptsum)<20 &&  fabs(b1jetvq)<10 &&  fabs(b2jetvq)<10 ) {
	  if( b1jetvq<0 ){
	    topquarkangleqpt->Fill(top1.CosTheta());
	    if(b2jetvq>0)topquarkangleqpt->Fill(-1*top2.CosTheta());
	  }
	  if( b1jetvq>0 ){
	    if(b2jetvq<0) topquarkangleqpt->Fill(top2.CosTheta());
	    topquarkangleqpt->Fill(-1*top1.CosTheta());
	  }
	}
	_tree->Fill();
      }
  }  catch(DataNotAvailableException &e) {
    cout << "WARNGING can not find collection with name " << _myJetsCollection << endl;
  }
}
        

void ttbar::check(LCEvent * evt){ 
  
}


void ttbar::end()
{ 
  ROOTfile->Write();
  ROOTfile->Close();
  delete ROOTfile;
  std::cout << "ttbarProcessor::end()  " << endl; 
  cout << "MC tops " << ntops <<  endl;
  cout <<" MC b quarks = " <<  nbqs << endl;
  cout <<" MC c quarks = " <<  ncqs << endl;
  cout <<" MC Ws = " <<  nws << endl;
  cout <<" MC (hadronic)Ws = " <<  nwhs << endl;
  cout <<" MC (leptonic)Ws = " <<  nwls << endl;
  cout <<" MC quarks = " <<  nqs << endl;
}
