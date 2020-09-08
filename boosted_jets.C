///////////////////////////////////////////////////////
// Example program how to access jets that were found 
// in different reference frames
//
// Author     : Roman Kogler (roman.kogler@desy.de)
// Created    : 10.1.2008
// last update: 
// by         : 
//
///////////////////////////////////////////////////////
#include <stdlib.h>
#include <iostream>
// ROOT-header:
#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TObjArray.h>
#include <TApplication.h>
// H1 header:
#include "H1Skeleton/H1Tree.h"
#include "H1Steering/H1StdCmdLine.h"
#include "H1Pointers/H1FloatPtr.h"
#include "H1PhysUtils/H1BoostedJets.h"
#include "H1Steering/H1ErrorHandler.h"

#include "H1Steering/H1SteerManager.h"


#include "H1HadronicCalibration/H1HadronicCalibration.h"

//added the H1 calculator ones:
#include "H1Calculator/H1Calculator.h"
#include "H1Calculator/H1CalcElec.h"
#include "H1Calculator/H1CalcFs.h"
#include "H1Calculator/H1CalcKine.h"
#include "H1Calculator/H1CalcWeight.h"
#include "H1Calculator/H1CalcVertex.h"
#include "H1Calculator/H1CalcGenericInterface.h"

#include "H1Cuts/H1CutEventSelector.h"

#include <limits.h>
using namespace std;
using namespace H1CalcGenericInterface;

int main(int argc, char* argv[])
{
    
    // Some general cosmetics:
    gROOT->SetStyle("Plain");
    gStyle->SetHistLineWidth(2);

    // parse the command line to give the option, e.g. nr of events to 
    // be processed with dstar_mods -n 1000 
    // the 2 lines are not needed if no cmd line options are given
    H1StdCmdLine opts;
    opts.Parse(&argc, argv);
    
    // needed for graphic, but MUST be AFTER Parse(...):
    TApplication theApp("boosted_jets", &argc, argv);

    H1Tree::Instance()->Open(); //must be there!

    // needed kinematics for the boost:
    H1FloatPtr Q2("Q2s");     // viruality
    H1FloatPtr Y("Ys");       // inelasticity
    H1FloatPtr Q2Gen("Q2sGen"); // generated viruality
    H1FloatPtr YGen("YsGen");   // generated inelasticity
    H1FloatPtr Ee("EBeamE");  // electron beam energy
    H1FloatPtr Ep("EBeamP");  // proton beam energy

    H1FloatPtr Pol("Pol"); //polarization of the event 
 
    // electron quantities needed for the boost:
    H1FloatPtr ElecE("ElecE");           // energy of scattered electron
    H1FloatPtr ElecTheta("ElecTheta");   // theta of scattered electron
    H1FloatPtr ElecPhi("ElecPhi");       // phi of scattered electron

                                          
    H1FloatPtr genEnElec("GenEnElec"); //   Electron energy, combined with photon for FSR
    H1FloatPtr genPhElec("GenPhElec"); //   Electron phi, combined with photon for FSR
    H1FloatPtr genThElec("GenThElec"); //   Electron theta, combined with photon for FSR
                                       
    H1ShortPtr runtype("RunType"); // 
    H1IntPtr runnumber("RunNumber");//

    H1BoostedJets* boostedjets = H1BoostedJets::Instance();
 
    H1ErrorHandler::Instance()->SetMaxCount(5); // mop up those error messages

    H1CutEventSelector* DISSelector = new H1CutEventSelector("DISSelector");
    DISSelector->AddRecCut(new H1CutFloat(Kine_Q2es,       150.0, 150000.0), "Q^2 Rec");
    DISSelector->AddRecCut(new H1CutFloat(Kine_Yes,         0.08,   0.94), "Y Rec");
    DISSelector->AddRecCut(new H1CutFloat(Fs_Empz,       45.0, 65.0), "E-pz");
  
    // Electron : Scattered Electron in LAr with Energy > 11 GeV
    DISSelector->AddRecCut(new H1CutDiscreteInt(Elec_Type, 1), "Elec_Type");     // 1 for LAr, 4 for SpaCal
    DISSelector->AddRecCut(new H1CutBool(Elec_IsScatteredElectron), "ScattElec");
    DISSelector->AddRecCut(new H1CutLorentz(Elec_FirstElectron, H1CutLorentz::E, 11.0, FLT_MAX), "Electron Energy");



    H1CutEventSelector* TRIGGERSelector = new H1CutEventSelector("TRIGGERSelector");   
    // Trigger Cuts: Select St67 and St77
    H1Cut* SubTrig67 = new H1CutBool(Trig_L1ac_idx, 67);

    TRIGGERSelector->AddRecCut(new H1CutOr(new H1CutBool(Event_IsMC), SubTrig67), "SubTriggers");
 
  
    ////////////////////////Background selector. 
    const Float_t DegToRad = TMath::Pi()/180.0;
    H1CutEventSelector* BKGSelector = new H1CutEventSelector("BKGSelector");


    // colliding bunch
    BKGSelector->AddRecCut(new H1CutInt(Event_BunchType,         3,       3), "Bunch Type");


    H1Cut* FwdTheta   = new H1CutLorentz(Elec_FirstElectron, H1CutLorentz::Theta, -FLT_MAX, 30*DegToRad);
    H1Cut* CJCVtx     = new H1CutBool(Vertex_IsVertexFound_idx, H1CalcVertex::vtCJC);
    H1Cut* FwdRegion  = new H1CutAnd(FwdTheta, CJCVtx);
    H1Cut* OptimalVtx = new H1CutBool(Vertex_IsVertexFound_idx, H1CalcVertex::vtOptimalNC);

    BKGSelector->AddRecCut(new H1CutOr(OptimalVtx, FwdRegion) , "NCOptimalVertex");
    BKGSelector->AddRecCut(new H1CutFloat(Vertex_Z, -35.0, 35.0), "Z-Vertex");

    // ---------------- LAr fiducial volume cuts  ---------------------

    // electron Zimpact Criteria
    BKGSelector->AddRecCut(new H1CutOr(new H1CutFloat(Elec_Zimpact, -190.0,  15.0),  new H1CutFloat(Elec_Zimpact,   25.0,  FLT_MAX)), "Electron Zimpact");

    // Remove Cracks
    BKGSelector->AddRecCut(new H1CutOr(new H1CutFloat(Elec_PhiOctant, 2.0*DegToRad, 43.0*DegToRad),  new H1CutFloat(Elec_Zimpact, 300.0, FLT_MAX)), "PhiCracks");

    // ---------------- electron cuts ---------------------

    // Electron : Scattered Electron in LAr with Energy > 8 GeV (tighter cut in DIS Selector)
    BKGSelector->AddRecCut(new H1CutDiscreteInt(Elec_Type, 1), "Elec_Type");     // 1 for LAr, 4 for SpaCal
    BKGSelector->AddRecCut(new H1CutBool(Elec_IsScatteredElectron), "ScattElec");
    BKGSelector->AddRecCut(new H1CutLorentz(Elec_FirstElectron, H1CutLorentz::E, 8.0, FLT_MAX), "Electron Energy");

    // this cut gets rid of most of the photoproduction background due to soft, wrongly identified electrons in the forward direction
    BKGSelector->AddRecCut(new H1CutFloat(Kine_Ye,        0.0,     0.94), "Ye Rec Background Reduction");

    // ---------------- Anti QED Compton ------------------
    H1Cut* nEmParts = new H1CutInt(Elec_NEmParts, 2, INT_MAX);
    H1Cut* Energy1  = new H1CutLorentz(Elec_FirstElectron, H1CutLorentz::E, 4, FLT_MAX);
    H1Cut* Energy2  = new H1CutLorentz(Elec_SecondElectron, H1CutLorentz::E, 4, FLT_MAX);

    H1Cut* EnergySum = new H1CutFunction("[0]+[1]", 18.0, FLT_MAX);                                                                                                                                              
    EnergySum->AddLorentzVariable(0, Elec_SecondElectron, H1CutLorentz::E);
    EnergySum->AddLorentzVariable(1, Elec_FirstElectron, H1CutLorentz::E);
    H1Cut* Acoplanarity = new H1CutFunction("-cos(TMath::Abs([0]-[1]))", 0.95, FLT_MAX);                                                                                                                         
    Acoplanarity->AddLorentzVariable(0, Elec_SecondElectron, H1CutLorentz::Phi);                                                                                                                                 
    Acoplanarity->AddLorentzVariable(1, Elec_FirstElectron, H1CutLorentz::Phi);
    BKGSelector->AddRecCut(new H1CutNot(new H1CutAnd(nEmParts, Energy1, Energy2, EnergySum, Acoplanarity)), "AntiCompton");

    // ---------------------- NON EP BACKGROUND FINDERS ---------------------------
    H1Cut* Finder0 =  new H1CutBool(BgTiming_Ibg_idx, 0);
    H1Cut* Finder1 =  new H1CutBool(BgTiming_Ibg_idx, 1);
    H1Cut* Finder5 =  new H1CutBool(BgTiming_Ibg_idx, 5);
    H1Cut* Finder6 =  new H1CutBool(BgTiming_Ibg_idx, 6);
    H1Cut* Finder7 =  new H1CutBool(BgTiming_Ibg_idx, 7);
    // Finder 7 and  Pth/Pte < 0.1
    BKGSelector->AddRecCut(new H1CutNot(new H1CutAnd(Finder7, new H1CutFloat(Fs_HadPtElecPtRatio, -FLT_MAX, 0.1))), "BkgFinder1");

    // Finder 0 and 1 or by two finders out of 5-7 and  Pth/Pte >= 0.1
    H1Cut* Pair1  = new H1CutAnd(Finder0, Finder1);
    H1Cut* Pair2  = new H1CutAnd(Finder5, Finder6);
    H1Cut* Pair3  = new H1CutAnd(Finder5, Finder7);
    H1Cut* Pair4  = new H1CutAnd(Finder6, Finder7);
    H1Cut* Pairs  = new H1CutOr(Pair1, Pair2, Pair3, Pair4);

    BKGSelector->AddRecCut(new H1CutNot(new H1CutAnd(Pairs, new H1CutFloat(Fs_HadPtElecPtRatio, 0.1, FLT_MAX))), "BkgFinder2");

    // Finder ( 5 || 6 ) && Pth/Pte < 0.5
    H1Cut* Finders = new H1CutOr(Finder5, Finder6);

    BKGSelector->AddRecCut(new H1CutNot(new H1CutAnd(Finders, new H1CutFloat(Fs_HadPtElecPtRatio, -FLT_MAX, 0.5))), "BkgFinder3");
    /////////////////////////////////////////////////////
    

  // some settings for the boosted jets - play with them!
    //boostedjets->SetCalibrationMethod(H1BoostedJets::eHighPtJetsOnly);
    //boostedjets->SetReferenceFrame(H1BoostedJets::eBreit);
    

    H1HadronicCalibration *hadronicCalibration=H1HadronicCalibration::Instance();
    // hadronicCalibration->SetCalibrationMethod(H1HadronicCalibration::eIterative);   
    hadronicCalibration->SetCalibrationMethod(H1HadronicCalibration::eHighPtJet);
    // hadronicCalibration->ApplyHadronicCalibration(kFALSE); 
    hadronicCalibration->ApplyHadronicCalibration(kTRUE);
    //    boostedjets->SysShiftHFSEnergy(0.02);
    //boostedjets->SysShiftHFSEnergy(-0.02);
    

    boostedjets->SetReferenceFrame(H1BoostedJets::eLab);
    boostedjets->UseFastJet();
    boostedjets->SetFastJetFinder(H1PartJet::eKt); ///eAntiKt
    boostedjets->SetJetFinderRecom(H1PartJet::ePt); //was H1PartJet::ePt //eE

    // if you would like to use an inclusive jet sample, where each jet 
    // should fulfill the requirements 
    //                   pt_min < pt_boostedjet < pt_max 
    //                               AND 
    //                eta_lab_min < eta_labjet < eta_lab_max 
    // you can use the CutFlags:
    boostedjets->SetPtCut(7.0, 50.0); 
    boostedjets->SetEtaLabCut(-1.0, 2.5);
    

    // open file to store histograms:
    TFile* file = new TFile(opts.GetOutput(), "RECREATE"); 
  

    // book some jet histograms
    TH1F* pthisto = new TH1F("pthisto", "P_{T} of inclusive jets in boosted frame", 50, 0, 50);
    TH1F* pthisto_lab = new TH1F("pthisto_lab", "P_{T} of inclusive jets in the lab frame", 50, 0, 50);

    TH1F* h_qt = new TH1F("h_qt", "qT of inclusive jets in the lab frame", 100, 0, 10);
    TH1F* h_z = new TH1F("h_z", "z spectrum of inclusive jets", 100,0,2.0);
    TH1F* h_zh = new TH1F("h_zh", "zh spectrum of hadrons-in-jet for inclusive jets", 100,0,2.0); 
    TH1F* h_jt = new TH1F("h_jt", "jt spectrum of hadrons-in-jet for inclusive jets, for 01<zh<0.5", 100,0.0,5.0);
    TH1F* dphi_lab = new TH1F("dphi_lab", "delta-phi of electrons and jets in the laboratory frame", 100, 0, TMath::Pi());
    TH1F* etahisto = new TH1F("etahisto", "#eta of inclusive jets in laboratory rest frame", 50, -1.5, 3);
   


    TTree *T = new TTree("T","jet Tree");

    //Tree variables:
    Float_t isBKG, passedDIScuts, passedTrigger ;
    Float_t event_x, event_y, event_Q2, event_genQ2, event_geny;
    Float_t event_x_es, event_y_es, event_Q2_es;
    Float_t polarization;
    
    Float_t e_pt, e_phi, e_rap, e_eta, e_theta, e_p;
    Float_t gene_pt, gene_phi, gene_rap, gene_eta, gene_theta, gene_p;
    Float_t vertex_z, ptmiss, Empz, Weight, WeightGen, ptratio,acoplanarity; 

    Float_t njets;  
    std::vector<UInt_t> nconstituents;
    std::vector<UInt_t> n_charged;
    std::vector<float> jet_pt;
    std::vector<float> jet_qt;
    std::vector<float> jet_phi;
    std::vector<float> jet_rap;
    std::vector<float> jet_eta;
    std::vector<float> jet_theta;
    std::vector<float> jet_p;
    std::vector<float> jet_z;
    std::vector<float> jet_dphi;

    std::vector<float> genjet_pt;
    std::vector<float> genjet_qt;
    std::vector<float> genjet_phi;
    std::vector<float> genjet_rap;
    std::vector<float> genjet_eta;
    std::vector<float> genjet_theta;
    std::vector<float> genjet_p;
    std::vector<float> genjet_z;
    std::vector<float> genjet_dphi;

    T->Branch("isBKG", &isBKG, "isBKG/F");
    T->Branch("passedDIScuts", &passedDIScuts, "passedDIScuts/F");
    T->Branch("passedTrigger", &passedTrigger, "passedTrigger/F");
    T->Branch("polarization", &polarization, "polarization/F");
    T->Branch("Weight", &Weight, "Weight/F");
    T->Branch("WeightGen", &WeightGen, "WeightGen/F");

    T->Branch("x", &event_x, "event_x/F");
    T->Branch("y", &event_y, "event_y/F");
    T->Branch("Q2", &event_Q2, "event_Q2/F");

    T->Branch("x_e", &event_x_es, "event_x_es/F");                                                                                                                                                                
    T->Branch("y_e", &event_y_es, "event_y_es/F");    
    T->Branch("Q2_e", &event_Q2_es, "event_Q2_es/F");


    T->Branch("genQ2", &event_genQ2, "event_genQ2/F");
    T->Branch("geny", &event_geny, "event_geny/F");

    T->Branch("vertex_z", &vertex_z, "vertex_z/F");
    T->Branch("ptmiss", &ptmiss, "ptmiss/F");
    T->Branch("ptratio", &ptratio, "ptratio/F");
    T->Branch("acoplanarity", &acoplanarity, "acoplanarity/F");
    T->Branch("Empz", &Empz, "Empz/F");
   
   
    T->Branch("e_pt", &e_pt, "e_pt/F");
    T->Branch("e_phi", &e_phi, "e_phi/F");
    T->Branch("e_rap",&e_rap, "e_rap/F");
    T->Branch("e_eta", &e_eta, "e_eta/F");
    T->Branch("e_p", &e_p, "e_p/F");
    T->Branch("e_theta", &e_theta, "e_theta/F");  

    T->Branch("gene_pt", &gene_pt, "gene_pt/F");                                                                                                                                                                  
    T->Branch("gene_phi", &gene_phi, "gene_phi/F");                                                                                                                                                              
    T->Branch("gene_rap",&gene_rap, "gene_rap/F");                                                                                                                                                                
    T->Branch("gene_eta", &gene_eta, "gene_eta/F");
    T->Branch("gene_p", &gene_p, "gene_p/F");
    T->Branch("gene_theta", &gene_theta, "gene_theta/F");

    T->Branch("njets", &njets, "njets/F");
    T->Branch("n_total",&nconstituents);
   
    T->Branch("jet_pt", &jet_pt);
    T->Branch("jet_qt", &jet_qt);
    T->Branch("jet_phi", &jet_phi);
    T->Branch("jet_rap",&jet_rap);
    T->Branch("jet_eta", &jet_eta);
    T->Branch("jet_theta", &jet_theta);
    T->Branch("jet_dphi", &jet_dphi);
    T->Branch("jet_p", &jet_p);
    T->Branch("jet_z", &jet_z);


    T->Branch("genjet_pt", &genjet_pt);
    T->Branch("genjet_qt", &genjet_qt);                                                                                                                                                                           
    T->Branch("genjet_phi", &genjet_phi);
    T->Branch("genjet_rap",&genjet_rap);
    T->Branch("genjet_eta", &genjet_eta);
    T->Branch("genjet_theta", &genjet_theta);
    T->Branch("genjet_dphi", &genjet_dphi);
    T->Branch("genjet_p", &genjet_p);
    T->Branch("genjet_z", &genjet_z);

    // open the canvas before the event loop
    TCanvas *canvas = new TCanvas("boosted_jets", "Plots", 10, 10, 600, 800);
    canvas->Divide(1, 3);

    // start event selection
    //cout << "Starting HAT selection: Q2 > 150" << endl;
    //std::cout << H1Tree::Instance()->SelectHat("Q2e > 150 && (Ye>0.2 && Ye<0.7)") << " events selected!" << std::endl;

    std::cout << H1Tree::Instance()->SelectHat("Q2e > 150") << " events selected!" << std::endl;

    if(hadronicCalibration->GetUseCalib()){                                                                                                                                                                      
        std::cout << "using calibration" << std::endl;                                                                                                                                                           
        }
    else{
      std::cout << "NOT USING CALIBRATION" << std::endl;                                                                                                                                                         
      }

    gH1Calc = H1Calculator::Instance();
    //Loop over events    
    Int_t events = 0;
    while (H1Tree::Instance()->Next() && (!opts.IsMaxEvent(events))) {
      //std::cout << *runnumber << " " << *runtype << std::endl;   
      //  H1HadronicCalibration::Instance()->Print();
      gH1Calc->Reset();
      boostedjets->Reset();
      H1Calculator::Instance()->Vertex()->SetPrimaryVertexType(H1CalcVertex::vtOptimalNC);
   
      //event selection:                                                                                                                                                                                                 //min Q2
      //if(gH1Calc->Kine()->GetQ2es()<150.0) continue;

      //std::cout << gH1Calc->Kine()->GetQ2es() << " " <<  gH1Calc->Kine()->GetQ2() << " "<< std::endl;
      if (!DISSelector->RecCuts()) continue;

      if (DISSelector->RecCuts()){
	passedDIScuts = 1.0;
      }
      else {
        passedDIScuts =0.0;
      }
      if(!BKGSelector->RecCuts()){
	isBKG=1.0;
      }
      else{
        isBKG=0.0;
      }
      if(TRIGGERSelector->RecCuts()){
        passedTrigger =1.0;
      }
      else{
        passedTrigger = 0.0;
      }
      //std::cout << isBKG << std::endl;

      //std::cout << gH1Calc->Kine()->GetQ2es() << " " << std::endl;

      //event selection:
      //min Q2   
      //      if(gH1Calc->Kine()->GetQ2es()<150.0) continue;
      //
      vertex_z = gH1Calc->Vertex()->GetHatPrimaryVertexZ();
      //if(fabs(vertex_z)>35 or vertex_z==0) continue;

      polarization =-999;      

      event_x = 0;
      event_y = 0;
      event_Q2 = 0;

      event_x_es = 0;
      event_y_es =0;
      event_Q2_es = 0;

      event_genQ2 = 0;
      event_geny = 0;

      vertex_z = 0;
      Empz = 0;
      ptmiss = 0;
      ptratio = 0;  
      acoplanarity = 0;

      Weight=0;
      WeightGen=0;

      jet_pt.clear();
      jet_qt.clear();
      jet_phi.clear();
      jet_rap.clear();
      jet_eta.clear();
      jet_theta.clear();
      jet_p.clear();
      jet_z.clear();
      jet_dphi.clear();

      genjet_pt.clear();
      genjet_qt.clear();
      genjet_phi.clear();
      genjet_rap.clear();
      genjet_eta.clear();
      genjet_theta.clear();
      genjet_p.clear();
      genjet_z.clear();
      genjet_dphi.clear();
          

      e_pt = 0;
      e_phi = 0;
      e_rap = 0; 
      e_eta = 0; 
      e_theta = 0;
      e_p = 0;

      gene_pt = 0;
      gene_phi = 0;
      gene_rap = 0;
      gene_eta = 0;
      gene_theta = 0;
      gene_p  = 0;

      njets =0;
      nconstituents.clear();
      n_charged.clear();

      // first thing to do in the eventloop:
      // reset the boosted jets 
   

      
      // and set reconstructed variables for the boost (must be there if a boost is used)
      // choose your preferred way of calculating these quantities:
      Float_t s = 4 * (*Ep) * (*Ee);
     

      Float_t ElecPx = *ElecE * TMath::Sin(*ElecTheta)*TMath::Cos(*ElecPhi);
      Float_t ElecPy = *ElecE * TMath::Sin(*ElecTheta)*TMath::Sin(*ElecPhi);
      Float_t ElecPz = *ElecE * TMath::Cos(*ElecTheta);

      //std::cout << "Reco electron energy " << *ElecE << " true " << *genEnElec << std::endl;
       
      Float_t genElecPx  = *genEnElec * TMath::Sin(*genThElec)*TMath::Cos(*genPhElec);
      Float_t genElecPy = *genEnElec * TMath::Sin(*genThElec)*TMath::Sin(*genPhElec);
      Float_t genElecPz = *genEnElec * TMath::Cos(*genThElec);

      TLorentzVector genScatteredElec(genElecPx, genElecPy, genElecPz, *genEnElec);  
      TLorentzVector ScatteredElec(ElecPx, ElecPy, ElecPz, *ElecE);
      TLorentzVector proton(0.0,0.0,*Ep, *Ep);
      TLorentzVector incoming_electron(0.0,0.0,-*Ee,*Ee);
      TLorentzVector virtual_photon = incoming_electron - ScatteredElec;
      TLorentzVector genvirtual_photon = incoming_electron - genScatteredElec;

      //event variables
      // Float_t Empz = gH1Calc->Fs()->GetEmpz();                           
      //Float_t PtMiss = gH1Calc->Fs()->GetPtMiss();         

      Empz = gH1Calc->Fs()->GetEmpz();  
      ptmiss = gH1Calc->Fs()->GetPtMiss();
      vertex_z = gH1Calc->Vertex()->GetHatPrimaryVertexZ();
      ptratio =   gH1Calc->Fs()->GetHadPtElecPtRatio();
      acoplanarity = gH1Calc->Fs()->GetElecAcoplanarity();
      Weight    = gH1Calc->Weight()->GetWeight();
      WeightGen = gH1Calc->Weight()->GetWeightGen();

      //std::cout << virtual_photon.M2() << " " << *Q2 << std::endl;
      //      boostedjets->SetXandScatElecVect(xbj,ScatteredElec);
      
      // if generator information is available: set generated variables for the boost
      //  Float_t xbjGen = gH1Calc->Kine()->GetXesGen();
      

      // level of jet finding (reconstructed, hadron or parton)
      H1BoostedJets::JetType reco = H1BoostedJets::eRecLev;
      H1BoostedJets::JetType truth = H1BoostedJets::eHadLev;
	
      // get the jets boosted back to the lab frame:
      TObjArray* jetarray = boostedjets->GetLabJets(reco);
      
      // boolean flags to see which jets fulfill the cuts:
      Bool_t* cutflags = boostedjets->GetCutFlags(reco);

      // this is the number of jets found without cuts:
      Int_t NumJets = boostedjets->GetNumJets(reco);

      TObjArray* truth_jetarray = boostedjets->GetLabJets(H1BoostedJets::eHadLev);
      Bool_t*   truth_cutflag = boostedjets->GetCutFlags(H1BoostedJets::eHadLev);
      Int_t truth_NumJets = boostedjets->GetNumJets(H1BoostedJets::eHadLev);

      event_x =  gH1Calc->Kine()->GetXes();
      event_Q2 = gH1Calc->Kine()->GetQ2es();
      event_y =  gH1Calc->Kine()->GetYes();

      event_x_es = gH1Calc->Kine()->GetXe();
      event_y_es = gH1Calc->Kine()->GetYe();          
      event_Q2_es = gH1Calc->Kine()->GetQ2e();
          
      event_genQ2 = gH1Calc->Kine()->GetQ2esGen();
      event_geny  = gH1Calc->Kine()->GetYesGen();
      //std::cout << genvirtual_photon.M2() << " " << event_genQ2 << " " << event_Q2<< std::endl;
      //  std::cout << " y " << event_y << " " << event_y_es <<std::endl;
      //std::cout << " Q2 " << event_Q2 << " " << event_Q2_es <<std::endl;   
       
      e_pt = ScatteredElec.Pt();
      e_phi = ScatteredElec.Phi();
      e_rap= ScatteredElec.Rapidity();
      e_eta = ScatteredElec.Eta();
      e_theta = ScatteredElec.Theta();
      e_p = ScatteredElec.P();

      gene_pt = genScatteredElec.Pt();
      gene_phi = genScatteredElec.Phi();                                                                                                                                                                          
      gene_rap= genScatteredElec.Rapidity();
      gene_eta = genScatteredElec.Eta();
      gene_theta = genScatteredElec.Theta();
      gene_p = genScatteredElec.P();

      njets = (Float_t)boostedjets->GetNumJetsAfterCuts(reco);

      polarization = *Pol;
      //  std::cout << "Polarization of the event " << polarization << std::endl;     
      // loop over the number of jets:
      for (Int_t ijet = 0; ijet<NumJets; ++ijet){
	//H1PartJet* jet = (H1PartJet*) boostedjetarray->At(ijet); 
	H1PartJet* jet = (H1PartJet*) jetarray->At(ijet);

        double jetphi = jet->GetPhi();
        double ephi = ScatteredElec.Phi();
	Float_t dphi = TMath::Abs(TVector2::Phi_mpi_pi(jetphi-ephi));

	TVector2 jet2(jet->GetPx(),jet->GetPy());
	TVector2 electron2(ElecPx,ElecPy);
	TVector2 qT = jet2 + electron2;

        // z variable (Lorentz invariant). 
        double z = proton.Dot(jet->GetFourVector())/proton.Dot(virtual_photon);
  
	//std::cout << ijet << std::endl;
	if (cutflags[ijet]){
	  // std::cout << ijet << std::endl;

	  //std::cout << " jet pt " << jet->GetPt() <<  " " << jet->GetEta() << std::endl;
	  pthisto->Fill(jet->GetPt());
          pthisto_lab->Fill(jet->GetPt());
          dphi_lab->Fill(dphi);
          h_qt->Fill(qT.Mod());
          h_z->Fill(z);
	  etahisto->Fill(jet->GetEta());
	   

          float deltaR = 999;
          int matched_index = -999;                                                       
          for(Int_t itrue = 0; itrue < truth_NumJets; ++itrue)
	    {																				
  	     H1PartJet* genjet = (H1PartJet*) truth_jetarray->At(itrue);                                                                                                                                     
             float distance = genjet->GetFourVector().DeltaR(jet->GetFourVector()) ;                                                                                                                       
             if( distance < deltaR)
	       {
		 deltaR = distance;
		 matched_index = itrue;
	       }
	    }	
          if(*runtype==0 or (*runtype!=0 and matched_index>-1)){

          jet_pt.push_back(jet->GetPt());
          jet_qt.push_back(qT.Mod());
	  jet_eta.push_back(jet->GetEta());
	  jet_p.push_back(jet->GetP());
	  jet_phi.push_back(jet->GetPhi());
	  jet_rap.push_back(jet->GetRapidity());
	  jet_theta.push_back(jet->GetTheta());
	  jet_dphi.push_back(TMath::Abs(TVector2::Phi_mpi_pi(jet->GetPhi()-ephi)));
	  jet_z.push_back(z);
	  nconstituents.push_back(jet->GetNumOfParticles());
	  
          if(*runtype>0 and matched_index>-1){
          H1PartJet* true_jet = (H1PartJet*) truth_jetarray->At(matched_index);

	  TVector2 genjet2(true_jet->GetPx(),true_jet->GetPy());
	  TVector2 genelectron2(genElecPx,genElecPy);
	  TVector2 genqT = genjet2 + genelectron2;

	  // z variable (Lorentz invariant).
	  double genz = proton.Dot(true_jet->GetFourVector())/proton.Dot(genvirtual_photon);
           
	  // std::cout << " truth jet pt eta" << true_jet->GetPt() << " " << true_jet->GetEta() << std::endl;
	  genjet_pt.push_back(true_jet->GetPt());
	  genjet_qt.push_back(genqT.Mod());
	  genjet_eta.push_back(true_jet->GetEta());
          genjet_p.push_back(true_jet->GetP());
          genjet_phi.push_back(true_jet->GetPhi());
	  genjet_rap.push_back(true_jet->GetRapidity());
	  genjet_theta.push_back(true_jet->GetTheta());
	  genjet_dphi.push_back(TMath::Abs(TVector2::Phi_mpi_pi(true_jet->GetPhi()- genScatteredElec.Phi())));
          genjet_z.push_back(genz);
	  }
	  //fill truth-level variables
	  }

	}	  //end if jet is selected loop

	//loop over constituents
	for (int n = 0; n < jet->GetNumOfParticles(); n++){
	  const H1Part* track = jet->GetParticle(n);
          if(abs(track->GetCharge())!=1) continue;
	  double zh = jet->GetFourVector().Vect().Dot( track->GetFourVector().Vect() )/(jet->GetFourVector().P()*jet->GetFourVector().P());
          h_zh->Fill(zh);
          double r = TMath::Sqrt( pow(jet->GetFourVector().Phi() - track->GetFourVector().Phi(),2.0) + pow(jet->GetFourVector().Eta() - track->GetFourVector().Eta(),2.0));
	  TVector3 zaxis(0,0,1);
          TVector3 N = zaxis.Cross(jet->GetFourVector().Vect());
	  TVector3 S = N.Cross(jet->GetFourVector().Vect());
	  N = N.Unit();
	  S = S.Unit();
	  TVector3 jt  = track->GetFourVector().Vect().Dot(N)*N + track->GetFourVector().Vect().Dot(S)*S;
	  if( z>0.1 and z<0.5){
	    h_jt->Fill(jt.Mag());
          }

	}

      }//end jet loop
	
      // this is the number of jets found after cuts:   
   
      
      
      if (events % 1000 == 0) {
	cout <<"event " << events <<" **********\n";
        canvas->cd(1);
        pthisto->Draw();
        //gPad->SetLogy();
        gPad->Update();
        canvas->cd(2);
        etahisto->Draw();
        //gPad->SetLogy();
        gPad->Update();
        canvas->cd(3);
       
        gPad->SetLogy();
        gPad->Update();
      }
      T->Fill();
      ++events;
      if (opts.IsMaxEvent(events)) break;
    }//end loop

    canvas->cd(1);
    pthisto->Draw();
    //gPad->SetLogy();
    gPad->Update();
    canvas->cd(2);
    etahisto->Draw();
    //gPad->SetLogy();
    gPad->Update();
    canvas->cd(3);
    
    gPad->SetLogy();
    gPad->Update();
    
    
    boostedjets->Print("Example Program");

    cout << "\nProcessed " << events
	 << " events and wrote jet histos to file "
	 << opts.GetOutput()  << "." << endl;

    pthisto->Write();
    pthisto_lab->Write();
    dphi_lab->Write();
    h_qt->Write();
    h_z->Write();
    h_zh->Write();
    h_jt->Write();
    etahisto->Write();
 
    T->Write();
    file->Write();
    file->Close();
}
