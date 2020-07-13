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

using namespace std;

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
    H1FloatPtr Q2("Q2e");     // viruality
    H1FloatPtr Y("Ye");       // inelasticity
    H1FloatPtr Q2Gen("Q2eGen"); // generated viruality
    H1FloatPtr YGen("YeGen");   // generated inelasticity
    H1FloatPtr Ee("EBeamE");  // electron beam energy
    H1FloatPtr Ep("EBeamP");  // proton beam energy
 
    // electron quantities needed for the boost:
    H1FloatPtr ElecE("ElecE");           // energy of scattered electron
    H1FloatPtr ElecTheta("ElecTheta");   // theta of scattered electron
    H1FloatPtr ElecPhi("ElecPhi");       // phi of scattered electron

    H1BoostedJets* boostedjets = H1BoostedJets::Instance();
 
    H1ErrorHandler::Instance()->SetMaxCount(5); // mop up those error messages

    // some settings for the boosted jets - play with them!
    //    boostedjets->SetCalibrationMethod(H1BoostedJets::eHighPtJetsOnly);
    //boostedjets->SetReferenceFrame(H1BoostedJets::eBreit);
    boostedjets->SetReferenceFrame(H1BoostedJets::eLab);
    boostedjets->UseFastJet();
    boostedjets->SetFastJetFinder(H1PartJet::eAntiKt);
    boostedjets->SetJetFinderRecom(H1PartJet::eE); //was H1PartJet::ePt

    // if you would like to use an inclusive jet sample, where each jet 
    // should fulfill the requirements 
    //                   pt_min < pt_boostedjet < pt_max 
    //                               AND 
    //                eta_lab_min < eta_labjet < eta_lab_max 
    // you can use the CutFlags:
    boostedjets->SetPtCut(5.0, 50.0); 
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
    TH1F* jetmultiplicity = new TH1F("jetmultiplicity", "Number of jets after cuts", 6, -0.5, 5.5);


    TTree *T = new TTree("T","jet Tree");

    //Tree variables:
    Float_t event_x, event_y, event_Q2;
    Float_t e_pt, e_phi, e_rap, e_eta, e_theta, e_p;
  
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

    T->Branch("x", &event_x, "event_x/F");
    T->Branch("y", &event_y, "event_y/F");
    T->Branch("Q2", &event_Q2, "event_Q2/F");
  
    //jet variables
    T->Branch("e_pt", &e_pt, "e_pt/F");
    T->Branch("e_phi", &e_phi, "e_phi/F");
    T->Branch("e_rap",&e_rap, "e_rap/F");
    T->Branch("e_eta", &e_eta, "e_eta/F");
    T->Branch("e_p", &e_p, "e_p/F");
    T->Branch("e_theta", &e_theta, "e_theta/F");  

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

    // open the canvas before the event loop
    TCanvas *canvas = new TCanvas("boosted_jets", "Plots", 10, 10, 600, 800);
    canvas->Divide(1, 3);

    // start event selection
    cout << "Starting HAT selection: Q2 > 100" << endl;
    cout << H1Tree::Instance()->SelectHat("Q2e > 100 && (Ye>0.2 && Ye<0.7)") << " events selected!" << endl;

    
    Int_t events = 0;
    while (H1Tree::Instance()->Next() && (!opts.IsMaxEvent(events))) {
      event_x = 0;
      event_y = 0;
      event_Q2 = 0;

      jet_pt.clear();
      jet_qt.clear();
      jet_phi.clear();
      jet_rap.clear();
      jet_eta.clear();
      jet_theta.clear();
      jet_p.clear();
      jet_z.clear();
      jet_dphi.clear();

      e_pt = 0;
      e_phi = 0;
      e_rap = 0; 
      e_eta = 0; 
      e_theta = 0;
      e_p = 0;

      nconstituents.clear();
      n_charged.clear();

      // first thing to do in the eventloop:
      // reset the boosted jets 
      boostedjets->Reset();
      
      // and set reconstructed variables for the boost (must be there if a boost is used)
      // choose your preferred way of calculating these quantities:
      Float_t s = 4 * (*Ep) * (*Ee);
      Float_t xbj = (*Q2)/((*Y) * s);

      Float_t ElecPx = *ElecE * TMath::Sin(*ElecTheta)*TMath::Cos(*ElecPhi);
      Float_t ElecPy = *ElecE * TMath::Sin(*ElecTheta)*TMath::Sin(*ElecPhi);
      Float_t ElecPz = *ElecE * TMath::Cos(*ElecTheta);
      TLorentzVector ScatteredElec(ElecPx, ElecPy, ElecPz, *ElecE);
      TLorentzVector proton(0.0,0.0,*Ep, *Ep);
      TLorentzVector incoming_electron(0.0,0.0,-*Ee,*Ee);
      TLorentzVector virtual_photon = incoming_electron - ScatteredElec;
      //std::cout << virtual_photon.M2() << " " << *Q2 << std::endl;
      boostedjets->SetXandScatElecVect(xbj,ScatteredElec);
      
      // if generator information is available: set generated variables for the boost
      Float_t xbjGen = (*Q2Gen)/((*YGen) * s);
      boostedjets->SetXGen(xbjGen);

      // level of jet finding (reconstructed, hadron or parton)
      H1BoostedJets::JetType level = H1BoostedJets::eRecLev;
	
      // get the jets in the boosted frame:
      TObjArray* boostedjetarray = boostedjets->GetBoostedJets(level);

      // get the jets boosted back to the lab frame:
      TObjArray* labjetarray = boostedjets->GetLabJets(level);

      // boolean flags to see which jets fulfill the cuts:
      Bool_t* cutflags = boostedjets->GetCutFlags(level);

      // this is the number of jets found without cuts:
      Int_t NumJets = boostedjets->GetNumJets(level);

      event_x = xbj;
      event_Q2 = *Q2;
      event_y = *Y;

      e_pt = ScatteredElec.Pt();
      e_phi = ScatteredElec.Phi();
      e_rap= ScatteredElec.Rapidity();
      e_eta = ScatteredElec.Eta();
      e_theta = ScatteredElec.Theta();
      e_p = ScatteredElec.P();
     
      // loop histograms:
      for (Int_t ijet = 0; ijet<NumJets; ++ijet){
	H1PartJet* jet = (H1PartJet*) boostedjetarray->At(ijet); 
	H1PartJet* labjet = (H1PartJet*) labjetarray->At(ijet);

        double jetphi = labjet->GetPhi();
        double ephi = ScatteredElec.Phi();
	Float_t dphi = TMath::Abs(TVector2::Phi_mpi_pi(jetphi-ephi));

	TVector2 jet2(labjet->GetPx(),labjet->GetPy());
	TVector2 electron2(ElecPx,ElecPy);
	TVector2 qT = jet2 + electron2;

        // z variable (Lorentz invariant). 
        double z = proton.Dot(labjet->GetFourVector())/proton.Dot(virtual_photon);
  
        //calculate delta-phi
	if (cutflags[ijet]){
	  pthisto->Fill(jet->GetPt());
          pthisto_lab->Fill(labjet->GetPt());
          dphi_lab->Fill(dphi);
          h_qt->Fill(qT.Mod());
          h_z->Fill(z);
	  etahisto->Fill(labjet->GetEta());

          jet_pt.push_back(labjet->GetPt());
          jet_qt.push_back(qT.Mod());
          jet_eta.push_back(labjet->GetEta());
          jet_p.push_back(labjet->GetP());
          jet_phi.push_back(labjet->GetPhi());
          jet_rap.push_back(labjet->GetRapidity());
          jet_theta.push_back(labjet->GetTheta());
          jet_dphi.push_back(TMath::Abs(TVector2::Phi_mpi_pi(labjet->GetPhi()-ephi)));
          jet_z.push_back(z);
	}	  

        nconstituents.push_back(labjet->GetNumOfParticles());
	//loop over constituents
	for (int n = 0; n < labjet->GetNumOfParticles(); n++){
	  const H1Part* track = labjet->GetParticle(n);
          if(abs(track->GetCharge())!=1) continue;
	  double zh = labjet->GetFourVector().Vect().Dot( track->GetFourVector().Vect() )/(labjet->GetFourVector().P()*labjet->GetFourVector().P());
          h_zh->Fill(zh);
          double r = TMath::Sqrt( pow(labjet->GetFourVector().Phi() - track->GetFourVector().Phi(),2.0) + pow(labjet->GetFourVector().Eta() - track->GetFourVector().Eta(),2.0));
	  TVector3 zaxis(0,0,1);
          TVector3 N = zaxis.Cross(labjet->GetFourVector().Vect());
	  TVector3 S = N.Cross(labjet->GetFourVector().Vect());
	  N = N.Unit();
	  S = S.Unit();
	  TVector3 jt  = track->GetFourVector().Vect().Dot(N)*N + track->GetFourVector().Vect().Dot(S)*S;
	  if( z>0.1 and z<0.5){
	    h_jt->Fill(jt.Mag());
          }

	}

      }//end jet loop
	
      // this is the number of jets found after cuts:
      Int_t NumJetsAfterCuts = boostedjets->GetNumJetsAfterCuts(level);
      jetmultiplicity->Fill(NumJetsAfterCuts);

      
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
        jetmultiplicity->Draw();
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
    jetmultiplicity->Draw();
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
    jetmultiplicity->Write();
    T->Write();
    file->Write();
    file->Close();
}
