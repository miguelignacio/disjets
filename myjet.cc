// main71.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace Pythia8;

int main() {
  // Settings
  int  nEvent = 1000000;
  const double etaMax  = 5.0;    // Pseudorapidity range of detector.

  Pythia pythia;
  Event& event = pythia.event;

  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ)=on");
  pythia.readString("PhaseSpace:Q2Min = 150");
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL=off");
  pythia.readString("SpaceShower:pTmaxMatch=2");
  pythia.readString("SpaceShower:dipoleRecoil=on");

    pythia.readString("HadronLevel:Hadronize = off");

  
  // Initialisation, p pbar @ 1.96 TeV
  pythia.readString("Beams:idB = 11");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:eA  = 920");
  pythia.readString("Beams:eB  = 27.6");
  pythia.readString("Beams:frameType = 2");
  pythia.init();

  // Histograms
  Hist dSigma1("Jet pt cross-section", 6, 10, 100.0,true);
  Hist dSigma_eta("Jet pseudorapidity cross-section", 5, -1.0,2.5);
  Hist dSigma_qt("Electron-jet momentum imbalance ", 20,1,6,true);

  // Fastjet analysis - select algorithm and parameters
  double Rparam = 1.0;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::kt_algorithm, Rparam,
                                      recombScheme, strategy);
  // Fastjet input
  std::vector <fastjet::PseudoJet> fjInputs;

  // Statistics for later

  bool firstEvent = true;
 
  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Reset Fastjet input
    fjInputs.resize(0);

   //general event info
    float evid = iEvent;
    float xsec = pythia.info.sigmaGen();
    float ntrials = pythia.info.nTried();

    // four-momenta of proton, electron, virtual photon/Z^0/W^+-.
    Vec4 pProton = event[1].p();
    Vec4 peIn = event[4].p();
    Vec4 peOut = event[6].p();
    Vec4 pPhoton = peIn - peOut;
    Vec4 ptot_had = pPhoton + pProton;
    
   // Q2, W2, Bjorken x, y, nu.
    float Q2 = -pPhoton.m2Calc();
    float W2 = (pProton + pPhoton).m2Calc();
    float x = Q2 / (2. * pProton * pPhoton);
    float y = (pProton * pPhoton) / (pProton * peIn);


    //DIS selection
    if(y>0.7 or y<0.2) continue;
    if(Q2<150.0) continue;
    
    Vec4   pTemp;
    double mTemp;
    int nAnalyze = 0;
    //loop over particles in the event and store them as input for FastJet
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {

      // Require visible/charged particles inside detector.
      if (!pythia.event[i].isFinal())        continue;
   // No neutrinos                                                                                                                                                                                       
      if (pythia.event[i].idAbs() == 12 || pythia.event[i].idAbs() == 14 ||
          pythia.event[i].idAbs() == 16)     continue;
      if(pythia.event[i].idAbs() ==11) continue;
      
      if ( !event[i].isVisible()  ) continue;
      if ( event[i].mother1()==6 ) continue; //remove scattered lepton 
      // if(event[i].e()==event[6].e()) continue; 
      if (abs(event[i].eta()) > etaMax) continue;
      
      // Create a PseudoJet from the complete Pythia particle.
      fastjet::PseudoJet particleTemp(event[i].px(), event[i].py(), event[i].pz(),event[i].e());
      if (particleTemp.pt()<0.100) continue;

      //particleTemp.set_user_info(new MyUserInfo(event[i].id(),i,event[i].charge()));
           
      fjInputs.push_back( particleTemp);
      ++nAnalyze;
    } //end loop over particles


    if (fjInputs.size() == 0) {
      cout << "Error: event with no final state particles" << endl;
      continue;
    }

    // Run Fastjet algorithm
    vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);

    // For the first event, print the FastJet details
    if (firstEvent) {
      cout << "Ran " << jetDef->description() << endl;
      cout << "Strategy adopted by FastJet was "
           << clustSeq.strategy_string() << endl << endl;
      firstEvent = false;
    }

    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV)
    inclusiveJets = clustSeq.inclusive_jets(10.0);
    sortedJets    = sorted_by_pt(inclusiveJets);

    /*    // Keep track of jets with pT > 20/25 GeV
    int  jetCount20 = 0, jetCount25 = 0;
    // For the deltaR calculation below
    bool vetoEvent = false;
    fastjet::PseudoJet fjElec(pythia.event[idxElec].px(),
                              pythia.event[idxElec].py(),
                              pythia.event[idxElec].pz(),
                              pythia.event[idxElec].e());
    */
    
    for (unsigned int i = 0; i < sortedJets.size(); i++) {
      // Only count jets that have |eta| < 2.0
      if (abs(sortedJets[i].eta()) > 2.5) continue;
      // Check distance between W decay electron and jets
      //if (fjElec.squared_distance(sortedJets[i]) < 0.52 * 0.52)
      //  { vetoEvent = true; break; }

      // Fill dSigma histograms and count jets with ET > 25.0
      if (sortedJets[i].perp() > 10.0){
	dSigma1.fill(sortedJets[i].perp());
	dSigma_eta.fill(sortedJets[i].eta());
	float qt = sqrt((peOut.px() + sortedJets[i].px())*(peOut.px() + sortedJets[i].px())+ (peOut.py() + sortedJets[i].py())* (peOut.py() + sortedJets[i].py()));
	dSigma_qt.fill(qt/sqrt(Q2)+1.0);
      }
    }//end of loop over jets
		    
    

  // End of event loop.
  } 

  // Statistics
  pythia.stat();

  // Output histograms
  //double sigmapb = pythia.info.sigmaGen() * 1.0E9;
  
  //dSigma1 = dSigma1*sigmapb/nEvent;
  
  //for (int i = 1; i <= 4; i++)
  //  (*dSigmaHist[i]) = ((*dSigmaHist[i]) * sigmapb) / nEvent / dSigmaBin[i];
  cout << dSigma1;
  cout << dSigma_eta;
  cout << dSigma_qt;
  int nentries = dSigma1.getEntries();
  dSigma1 = dSigma1/nentries;
  dSigma_eta = dSigma_eta/nentries;
  dSigma_qt = dSigma_qt/nentries;
  
  HistPlot hpl( "jetpt"); 
  hpl.frame( "pTdist", "Boson pT distributions", "pT (GeV)", "sigma"); 

  hpl.add( dSigma1, "o-");
  hpl.plot(true);

  HistPlot hpl_eta( "jeteta");
  hpl_eta.frame( "jeteta", "eta distributions", "jet eta", "sigma");

  
  hpl_eta.add(dSigma_eta,"-o");
  hpl_eta.plot();

  HistPlot hpl_qt( "jetqt");
  hpl_qt.frame( "jetqt", "qt distributions", "jet qt", "sigma");


  hpl_qt.add(dSigma_qt,"-o");
  hpl_qt.plot(true);


  


  // Done.
  delete jetDef;
  return 0;
}
