// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "HepMC/GenParticle.h"
#include <iostream>

using HepMC::GenParticle;
using HepMC::GenVertex;

namespace Rivet {

  const double leading_lep_min_pT = 25*GeV;
  const double trailing_lep_min_pT = 13*GeV;
  const double lep_max_eta = 2.5*GeV;
  const double dilep_min_mass = 12*GeV;
  const double dilep_min_pT = 30*GeV;
  const double trailing_lep_min_mT = 30*GeV;
  const double higgs_min_mT = 60*GeV;

  const double jet_min_pT = 30*GeV;
  const double jet_max_eta = 4.7;

  const double lep_cone_size = 0.1;

  const bool acceptance_on = 0;

  bool cuts_0 [1];
  bool cuts_1 [3];
  bool cuts_2 [4];

  bool Zcuts [2];
  bool Wcuts [12];

  double m_sumw;

  /// @brief Add a short analysis description here
  class VH_WW : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(VH_WW);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // the basic final-state projection: 
      // all final-state particles within 
      // the given eta acceptance
      const FinalState fs;
      declare(fs, "fs");

      m_sumw = 0.0;

      // the final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      //FastJets jets(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      //declare(jets, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      declare(bare_leps, "bare_leps"); 
      // dress the prompt bare leptons with prompt photons within dR < 0.1
      DressedLeptons dressed_leps(photons, bare_leps, lep_cone_size);
      declare(dressed_leps, "leptons");

      //missing momentum
      MissingMomentum Met(fs);
      declare(Met, "Met");

      //book histograms
      book(H_PT, "H_PT", {0.0, 20.0, 45.0, 80.0, 120.0, 200.0});
      book(N_jets, "N_jets", {0, 1, 2, 3, 4});
      book(acc_H_PT, "acc_H_PT", {0.0, 20.0, 45.0, 80.0, 120.0, 200.0});
      book(acc_N_jets, "acc_N_jets", {0, 1, 2, 3, 4});

      book(leading_lep_PT, "leading_lep_PT", 50,0,100);
      book(trailing_lep_PT, "trailing_lep_PT", 50,0,50);
      //book(leading_lep_eta, "leading_lep_eta", 25,0,5);
      //book(trailing_lep_eta, "trailing_lep_eta", 25,0,5);
      book(dilep_mass, "dilep_mass", 50,0,150);
      book(dilep_PT, "dilep_PT", 25,0,200);
      //book(trailing_lep_mT_hist, "trailing_lep_mT", 25,0,120);
      //book(higgs_mT_hist, "H_mT", 25,0,140);
      book(dphi_hist, "delta_phi", 25,0,3.14);
      //book(cos_dphi, "cos_dphi", 25,-1,1);  

      book(acceptance, "acceptance", {0,1,2,3,4,5,6,7,8,9});
      book(lepton_no, "lepton_no", {0,1,2,3,4,5,6,7,8,9,10});
      //book(V_mass, "V_mass", 50, 0, 200);
      //book(V_PT, "V_PT", {0, 75, 150, 250, 400});
      //book(acc_V_PT, "acc_V_PT", {0, 75, 150, 250, 400});

      book(total_lepton_charge, "total_lepton_charge", {-2, -1, 0, 1, 2});
      book(lep3_PT, "lep3_PT", 50, 0, 50);
      book(additional_lepton_PT, "additional_lepton_PT", 25, 0, 30);
      book(X_mass, "X_mass", 50, 0, 100);
      book(max_b_jet_PT, "max_b_jet_PT", 25, 0, 50);
      book(dilepton_Z_mass_diff, "dilepton_Z_mass_diff", 50, -30, 30);
      book(missing_momentum, "missing_momentum", 50, 0, 200);
      book(total_lepton_mass, "total_lepton_mass", 50, 0, 250); 
      book(min_dilep_mass, "min_dilep_mass", 25, 0, 200);
      book(total_lep_charge, "total_lep_charge", {-2, -1, 0, 1, 2});
      book(leading_jet_PT, "leading_jet_PT", 25, 0, 100); 
      book(cut_hist, "cut_hist", {0,1,2,3,4,5,6,7,8,9,10,11,12});
     
	 
    }

    bool originateFrom(const Particle& p, const Particles& ptcls ) {
      const GenVertex* prodVtx = p.genParticle()->production_vertex();
      if (prodVtx == nullptr) return false;
      for (auto ancestor : Rivet::HepMCUtils::particles(prodVtx, HepMC::ancestors)){ 
          for ( auto part:ptcls ) 
              if ( ancestor==part.genParticle() ) return true;
          }
          return false; 
      }
     
    bool originateFrom(const Particle& p, const Particle& p2){
        Particles ptcls = {p2}; return originateFrom(p,ptcls);
    }

    /*
    void setCuts(const Event& event, Particles leptons, FourMomentum LL, double trailing_lep_mT, double higgs_mT) {
	//set cuts	
	cuts[0] = (leptons[0].charge() == leptons[1].charge()) || (leptons[0].abspid() == leptons[1].abspid());
	cuts[1] = leptons[0].pT() < leading_lep_min_pT;
	cuts[2] = leptons[1].pT() < trailing_lep_min_pT;
	cuts[3] = (leptons[0].abseta() > lep_max_eta || leptons[1].abseta() > lep_max_eta);
	cuts[4] = LL.mass() < dilep_min_mass;
	cuts[5] = LL.pT() < dilep_min_pT;
	cuts[6] = trailing_lep_mT < trailing_lep_min_mT;
	cuts[7] = higgs_mT < higgs_min_mT;
    }
    */   
    /*
    void LepCuts_0(const Event& event, FourMomentum p_miss) {
	//find extra variables for cutting                                      

	//set cuts
	cuts_0[0] = p_miss.pT() < 170*GeV;	
    }


    void LepCuts_1(const Event& event, FourMomentum p_miss, Particles temp_leptons, double dphi) {
	//find extra variables for cutting                                      

	//set cuts
	FourMomentum V = p_miss + temp_leptons[0];
	cuts_1[0] = V.pT() < 100*GeV;

	for (const Particle p : temp_leptons) {
          if (((p.pid()) == 11) || ((p.pid()) == -11)) {
		cuts_1[1] = temp_leptons[0].pT() < 30*GeV;
          } else if (((p.pid()) == 13) || ((p.pid()) == -13)) {
		cuts_1[1] = temp_leptons[0].pT() < 25*GeV;
	  } else {
		cuts_1[1] = 1;
	  }
	cuts_1[2] = dphi > 2; 
        }
     }

    void LepCuts_2(const Event& event, Particles temp_leptons) {
	//find extra variables for cutting                                      

	//set cuts
	bool pair_found = 0;
      	Particle leading_lep;
      	Particle trailing_lep;
      	if (temp_leptons.size()>1) {
        	leading_lep = temp_leptons[0];
        	for (const Particle p: temp_leptons) {
          		if (!((p.pid())==(-leading_lep.pid()))) continue;
          		trailing_lep = p;
          		pair_found = 1;
          		break;
        	}	
      	}

	if (pair_found) {
	  cuts_2[0] = 0;
          FourMomentum LL = (leading_lep.momentum() + trailing_lep.momentum());
          V_mass->fill(LL.mass());
        
	  for (const Particle p : temp_leptons){
	    if (((p.pid()) == 13) || ((p.pid()) == -13)) {
	      cuts_2[1] = (LL.pT() < 50*GeV || LL.pT() > 150*GeV);
	    } else if (((p.pid()) == 11) || ((p.pid()) == -11)) {
	      cuts_2[1] = LL.pT() < 150*GeV;
	    } else {
	      cuts_2[1]= 1;
	    }
	  cuts_2[2]= (leading_lep.pT() < 20* GeV) || (trailing_lep.pT() < 20*GeV);
	  cuts_2[3] = (LL.mass() < 75*GeV) || (LL.mass() > 105*GeV);
	  }
	} else {
	  cuts_2[0] = 1;
	}
	 
    }
    */



    //checks all cuts with one exception
    bool checkZCutsException(int exception_index) {
    	for (int i = 0; i < 2; i++) {
		if (i != exception_index) {
			if (Zcuts[i] == true) return true;
		}
	}
	return false;	
    }

    bool checkZCutsExceptions(int exception_index1, int exception_index2) {
        for (int i = 0; i < 13; i++) {
                if (i != exception_index1 && i != exception_index2) {
                        if (Zcuts[i] == true) return true;
                }
        }
        return false;
    }


    bool checkWCutsException(int exception_index) {
        for (int i = 0; i < 12; i++) {
                if (i != exception_index) {
                        if (Wcuts[i] == true) return true;
                }
        }
        return false;
    }

    bool checkWCutsExceptions(int exception_index1, int exception_index2) {
        for (int i = 0; i < 12; i++) {
                if (i != exception_index1 && i != exception_index2) {
                        if (Wcuts[i] == true) return true;
                }
        }
        return false;
    }


 

    bool checkWcutsPreselection() {
	for (int i = 0; i < 6; i++) {
        	if (Wcuts[i] == true) return true;
        }
        return false;
    }


    bool checkWcut(int cut_index) {
	if (Wcuts[cut_index] == true) {
		return true;
	}else {
		return false;
	}
    }



    bool checkWcuts() {
        for (bool cut : Wcuts) {
                if (cut) return true;
        }
        return false;
    }


    bool checkZcuts() {
        for (bool cut : Zcuts) {
                if (cut) return true;
        }
        return false;
    }
    
    void setCutsZH(Particles leptons, Particle Z, FourMomentum p_miss, Jets eta_cut_jets) {
	Particle Z_leps [2];
	Particles X_leps;

	FourMomentum Z_mom = Z.momentum();


	double min_dm = 9999*GeV;

	int lep_index1;
	int lep_index2;

	bool Z_pair_found = 0;

	for (unsigned int i=1; i<leptons.size(); i++) {
                for (unsigned int j=0; j<i; j++) {
                        if (leptons[i].pid() == -leptons[j].pid()) {
                                FourMomentum LL = (leptons[i].momentum() + leptons[j].momentum());
                                if (abs(LL.mass() - 91.2*GeV) < min_dm) {
					Z_pair_found = 1;
                                        min_dm = abs(LL.mass() - 91.2*GeV);
                                        lep_index1 = i;
                                        lep_index2 = j;
                                }
                        }
                }
        }
	
	Zcuts[0] = 0;
	if (!Z_pair_found) {
		Zcuts[0] = 1;
		return;
	}
		


	Particle Z_lep1 = leptons[lep_index1];
	Particle Z_lep2 = leptons[lep_index2];

	Z_leps[0] = Z_lep1;
	Z_leps[1] = Z_lep2;

	for (unsigned int i=0; i<leptons.size(); i++) {
			if (i != lep_index1 && i != lep_index2) {
				Particle p = leptons[i];				 
				X_leps += p;
			}
	} 
                
	

	FourMomentum LL = (Z_leps[0].momentum() + Z_leps[1].momentum());
	FourMomentum X_mom = (X_leps[0].momentum() + X_leps[1].momentum());
	FourMomentum Total_mom = (Z_leps[0].momentum() + Z_leps[1].momentum() + X_leps[0].momentum() + X_leps[1].momentum());
		
	/*
	Zcuts[1] = (leptons[0].charge() + leptons[1].charge() + leptons[2].charge() + leptons[3].charge() != 0);
	Zcuts[2] = (leptons[0].pT()  < 25*GeV);
	Zcuts[3] = (leptons[1].pT() < 15*GeV);
	Zcuts[4] = (leptons[2].pT() < 10*GeV || leptons[3].pT() < 10*GeV);
	if (leptons.size() > 4) {
		Zcuts[5] = (leptons[4].pT() > 10*GeV);
	} else {
		Zcuts[5] = 0;
	}

	Zcuts[6] = (LL.mass() < 4*GeV);
	Zcuts[7] = (X_mom.mass() < 4*GeV);
	Zcuts[8] = 0;
        if (eta_cut_jets.size() > 0) {
                for(unsigned int i = 0; i < eta_cut_jets.size(); i++) {
                        if (eta_cut_jets[i].bTagged()) {
                                if (eta_cut_jets[i].pT() > 20*GeV) {
                                                Zcuts[8] = 1;
                                                break;
                                }
                        }
                }
	}
	*/
	if (((X_leps[0].pid())==(-X_leps[1].pid()))) {
		//Zcuts[9] = (abs(LL.mass()- 91.2*GeV) > 15*GeV);
		//Zcuts[10] = (X_mom.mass() < 10*GeV || X_mom.mass() > 50*GeV);
		//Zcuts[11] = (p_miss.pT() < 35*GeV || p_miss.pT() > 100*GeV);
		//Zcuts[12] = (Total_mom.mass() < 140*GeV);
		Zcuts[0] = (p_miss.pT() < 35*GeV || p_miss.pT() > 100*GeV);
		Zcuts[1] = (Total_mom.mass() < 140*GeV);
		
	} else {
                //Zcuts[9] = (abs(LL.mass()- 91.2*GeV)  > 15*GeV);
                //Zcuts[10] = (X_mom.mass() < 10*GeV || X_mom.mass() > 70*GeV);
                //Zcuts[11] = (p_miss.pT() < 20*GeV);
                //Zcuts[12] = 0;
                Zcuts[0] = (p_miss.pT() < 20*GeV);
                Zcuts[1] = 0;
	}
	/*

	if (!checkZCutsException(1)) total_lepton_charge->fill(leptons[0].charge() + leptons[1].charge() + leptons[2].charge() + leptons[3].charge());
        if (!checkZCutsException(2)) leading_lep_PT->fill(leptons[0].pT()/GeV);
        if (!checkZCutsException(3)) trailing_lep_PT->fill(leptons[1].pT()/GeV);
        if (!checkZCutsException(4)) lep3_PT->fill(leptons[2].pT()/GeV);
        if (!checkZCutsException(5)) {
		if (leptons.size() > 4) {
			additional_lepton_PT->fill(leptons[4].pT()/GeV);
		} else {
			additional_lepton_PT->fill(0);
		}
	}
        if (!checkZCutsException(6)) dilep_mass->fill(LL.mass()/GeV);
        if (!checkZCutsExceptions(7,10)) X_mass->fill(X_mom.mass()/GeV);
        if (!checkZCutsException(8)) {
		double b_jet_PT = 0; 
		if (eta_cut_jets.size() > 0) {
                	for(unsigned int i = 0; i < eta_cut_jets.size(); i++) {
                        	if (eta_cut_jets[i].bTagged()) {
                                	b_jet_PT=eta_cut_jets[i].pT()/GeV; 
                                        break;
                                }
                        }
                }
		max_b_jet_PT->fill(b_jet_PT);
	}
	if (!checkZCutsException(9)) dilepton_Z_mass_diff->fill((LL.mass()- 91.2*GeV)/GeV);
	//if (!checkCutsException(10)) X_mass->fill(X_mom.mass()/GeV);
	if (!checkZCutsException(11)) missing_momentum->fill(p_miss.pT()/GeV);
	if (!checkZCutsException(12)) total_lepton_mass->fill(Total_mom.mass()/GeV);
	*/

	if (!checkZCutsException(0)) missing_momentum->fill(p_miss.pT()/GeV);
        if (!checkZCutsException(1)) total_lepton_mass->fill(Total_mom.mass()/GeV);
	
        //}
    }
   
    void setCutsWH(Particles leptons, FourMomentum p_miss, Particle Z, Jets eta_cut_jets) {
	FourMomentum Z_mom = Z.momentum();
	FourMomentum LL01 = (leptons[0].momentum()+leptons[1].momentum());
	FourMomentum LL12 = (leptons[1].momentum()+leptons[2].momentum());
	FourMomentum LL02 = (leptons[0].momentum()+leptons[2].momentum());
	FourMomentum LLL = (leptons[0].momentum()+leptons[1].momentum()+leptons[2].momentum());

	double min_ll = 9999;
	
	double min_dm = 9999;

	double dphi = deltaPhi(LLL, p_miss);

	//bool Z_pair_found = 0;


	for (unsigned int i=1; i<leptons.size(); i++) {
                for (unsigned int j=0; j<i; j++) {
                        if (leptons[i].pid() == -leptons[j].pid()) {
                                FourMomentum LL = (leptons[i].momentum() + leptons[j].momentum());
                                if (abs(LL.mass() - 91.2*GeV) < min_dm) {
                                        //Z_pair_found = 1;
                                        min_dm = abs(LL.mass() - 91.2*GeV);
                                        //lep_index1 = i;
                                        //lep_index2 = j;
                                }
                        }
                }
        }



        Wcuts[0] = (leptons[0].pT()  < 25*GeV);
	Wcuts[1] = (leptons[1].pT() < 20*GeV);
	Wcuts[2] = (leptons[2].pT() < 15*GeV);
	if (leptons.size() > 3) {
                Wcuts[3] = (leptons[3].pT() > 10*GeV);
        } else {
                Wcuts[3] = 0;
        }
	
	if (leptons[0].charge() == -leptons[1].charge()) {
		if (LL01.mass() < min_ll) {
			min_ll = LL01.mass();
		}
	}
	if (leptons[1].charge() == -leptons[2].charge()) {
		if (LL12.mass() < min_ll) {
			min_ll = LL12.mass();
		}
	}
	if (leptons[0].charge() == -leptons[2].charge()) {
                if (LL02.mass() < min_ll) {
                        min_ll = LL02.mass();
                }
        }
	Wcuts[4] = (min_ll < 12*GeV);
	Wcuts[5] = (abs(leptons[0].charge() + leptons[1].charge() + leptons[2].charge()) != 1);
	bool pair_found = 0;
	if (abs(leptons[0].pid()) == abs(leptons[1].pid()) &&  abs(leptons[1].pid()) == abs(leptons[2].pid())) {
		Wcuts[6] = 1;
		Wcuts[7] = 1;
		Wcuts[8] = 1;
		Wcuts[9] = 1;
		Wcuts[10] = 1;
		Wcuts[11] = 1;

	} else {
		for (int i=1; i<3; i++) {
			for (int j=0; j<i; j++) {
                        	if ((leptons[i].pid())==(-leptons[j].pid())) {
					FourMomentum LL = (leptons[i].momentum() +leptons[j].momentum());
                        		pair_found = 1;
                        		break;
				}	
			}
		}
	}
	
		if (pair_found == 1) {
			Wcuts[6] = 0;
			Wcuts[7] = 0;
			if (eta_cut_jets.size() > 0) {
				Wcuts[6] = eta_cut_jets[0].pT()  > 30*GeV;
                		for(unsigned int i = 0; i < eta_cut_jets.size(); i++) {
                        		if (eta_cut_jets[i].bTagged()) {
                                		if (eta_cut_jets[i].pT() > 20*GeV) {
                                                	Wcuts[7] = 1;
                                                	break;
                                		}	
                        		}			
                		}	
       	 		}
	
			Wcuts[8] = p_miss.pT() < 50*GeV;
			Wcuts[9] = min_ll > 100*GeV;
			Wcuts[10] = min_dm < 25*GeV;
			Wcuts[11] = dphi < 2.2;
	
		} else{
			Wcuts[6] = 0;
                        Wcuts[7] = 0;
                        if (eta_cut_jets.size() > 0) {
                                Wcuts[6] = eta_cut_jets[0].pT()  > 30*GeV;
                                for(unsigned int i = 0; i < eta_cut_jets.size(); i++) {
                                        if (eta_cut_jets[i].bTagged()) {
                                                if (eta_cut_jets[i].pT() > 20*GeV) {
                                                        Wcuts[7] = 1;
                                                        break;
                                                }
                                        }
                                }
                        }

                	Wcuts[8] = 0;
                	Wcuts[9] = 0;
                	Wcuts[10] = 0;
                	Wcuts[11] = dphi < 2.5;
		}

		
        if (!checkWCutsException(0)) leading_lep_PT->fill(leptons[0].pT()/GeV);
        if (!checkWCutsException(1)) trailing_lep_PT->fill(leptons[1].pT()/GeV);
        if (!checkWCutsException(2)) lep3_PT->fill(leptons[2].pT()/GeV);
        if (!checkWCutsException(3)) {
                if (leptons.size() > 3) {
                        additional_lepton_PT->fill(leptons[3].pT()/GeV);
                } else {
                        additional_lepton_PT->fill(0);
                }
	}
	if (!checkWCutsExceptions(4, 9)) min_dilep_mass ->fill(min_ll/GeV);
	if (!checkWCutsException(5)) total_lep_charge ->fill(abs(leptons[0].charge() + leptons[1].charge() + leptons[2].charge()));
	
        if (!checkWCutsException(6)) {
		if (eta_cut_jets.size() > 0) {
                                leading_jet_PT->fill(eta_cut_jets[0].pT()/GeV);
		}
	}
        if (!checkWCutsException(7)) {
                double b_jet_PT = 0;
                if (eta_cut_jets.size() > 0) {
                        for(unsigned int i = 0; i < eta_cut_jets.size(); i++) {
                                if (eta_cut_jets[i].bTagged()) {
                                        b_jet_PT=eta_cut_jets[i].pT()/GeV;
                                        break;
                                }
                        }
                }
                max_b_jet_PT->fill(b_jet_PT);
        }
        if (!checkWCutsException(8)) missing_momentum->fill(p_miss.pT()/GeV);
	
	if (!checkWCutsException(10)) {
		if (pair_found == 1) {
			dilepton_Z_mass_diff->fill(min_dm/GeV);
		} else {dilepton_Z_mass_diff->fill(0);
		}
	}
	if (!checkWCutsException(11)) dphi_hist->fill(dphi);
	
    	//}
    }
	
    /*
    //goes through cuts, if event failed any cuts, vetoEvent
    void applyVeto() {
        for (bool cut : cuts) {
		if (cut) vetoEvent;
	}
    }
    */
    /*
    bool checkCuts_0() {
	for (bool cut : cuts_0) {
		if (cut) return true;
	}
	return false;
    }

    bool checkCuts_1() {
        for (bool cut : cuts_1) {
                if (cut) return true;
        }
        return false;
    }

    bool checkCuts_2() {
        for (bool cut : cuts_2) {
                if (cut) return true;
        }
        return false;
    }
    */

    
    



    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weights()[0];
      m_sumw += weight;

      //grab higgs
      Particle higgs;
      for (const GenParticle *ptcl : Rivet::HepMCUtils::particles(event.genEvent())){
          if (!PID::isHiggs(ptcl->pdg_id())) continue;
          higgs = Particle(ptcl);
      }

      //grab vector boson
      Particle V;
      for (const GenParticle *ptcl : Rivet::HepMCUtils::particles(event.genEvent())){
          if (((ptcl->pdg_id()) != 23) && ((ptcl->pdg_id()) != 24) && ((ptcl->pdg_id()) != -24)) continue;
          V = Particle(ptcl);
      }

      Particle Z;
      for (const GenParticle *ptcl : Rivet::HepMCUtils::particles(event.genEvent())){
          if ((ptcl->pdg_id()) != 23) continue;
          Z = Particle(ptcl);
      }
      /*
      Particle W;
      for (const GenParticle *ptcl : Rivet::HepMCUtils::particles(event.genEvent())){
          if (((ptcl->pdg_id()) != 24) && ((ptcl->pdg_id()) != -24)) continue;
          W = Particle(ptcl);
      }
      */

      // grab jets not including particles from decay  
      Particles fs = apply<FinalState>(event, "fs").particles();
      Particles hadrons;
      for (const Particle &p : fs) {
          if (originateFrom(p,higgs) || originateFrom(p,V)) continue;
          hadrons += p;
      }
      FinalState fps_temp;
      FastJets jets(fps_temp, FastJets::ANTIKT, 0.4);
      jets.calc(hadrons); 
      Jets jets30 = jets.jetsByPt(Cuts::pT>jet_min_pT && Cuts::abseta<jet_max_eta);
      Jets eta_cut_jets = jets.jetsByPt(Cuts::abseta<jet_max_eta);

      //fill histograms before acceptance cuts
      H_PT->fill(higgs.pT()/GeV);
      //V_PT->fill(V.pT()/GeV);
      N_jets->fill(jets30.size());

      //start acceptance
      //Particles temp_leptons = apply<PromptFinalState>(event, "bare_leps").particlesByPt(Cuts::pT>20*GeV && Cuts::abseta<2.5);  
      Particles leptons = apply<DressedLeptons>(event, "leptons").particlesByPt();
      

      lepton_no->fill(leptons.size());

      FourMomentum p_miss = apply<MissingMomentum>(event, "Met").missingMomentum();
      								  
      //fill histograms (if event not been vetoed)
      /*       
      if (temp_leptons.size() == 0) {
	FourMomentum p_miss = apply<MissingMomentum>(event, "Met").missingMomentum();
        LepCuts_0(event, p_miss);
	if (checkCuts_0()) vetoEvent;
        acceptance->fill(0);
      } else if (temp_leptons.size() == 1) {
	FourMomentum p_miss = apply<MissingMomentum>(event, "Met").missingMomentum();
	double dphi = deltaPhi(temp_leptons[0], p_miss);
	LepCuts_1(event, p_miss, temp_leptons, dphi);
	if (checkCuts_1()) vetoEvent;
        acceptance->fill(1);
      } else {
	LepCuts_2(event, temp_leptons);
	if (checkCuts_2()) vetoEvent;
        acceptance->fill(2);
      } 

      */


      FourMomentum Z_mom = Z.momentum();

      if (leptons.size() == 3) {
	setCutsWH(leptons, p_miss, Z, eta_cut_jets);
	for (unsigned int i=0; i<12; i++) {
		if (checkWcut(i)) {
			cut_hist->fill(i);
		}
	}
      }
 

      if (leptons.size() > 3) {
        setCutsWH(leptons, p_miss, Z, eta_cut_jets);
	//setCutsZH(leptons, Z, p_miss, eta_cut_jets);
	if (checkWcutsPreselection())  {
		acceptance->fill(0);
		vetoEvent;
	}
      } else {if (leptons.size() > 2) {
                setCutsWH(leptons, p_miss, Z, eta_cut_jets);
                if (checkWcutsPreselection()) {
                        acceptance->fill(1);
                        vetoEvent;
                }
              } else {
                        acceptance->fill(2);
                        vetoEvent;
                }
      }
     
      if (leptons.size() > 3) {
        setCutsWH(leptons, p_miss, Z, eta_cut_jets);
	setCutsZH(leptons, Z, p_miss, eta_cut_jets);
        if (!checkWcuts() && !checkZcuts()) {
		acceptance->fill(3);
                //vetoEvent;
	} if (checkWcuts() && !checkZcuts()) {
		acceptance->fill(4);
                //vetoEvent;
        }
	if (!checkWcuts() && checkZcuts()) {
                acceptance->fill(5);
                //vetoEvent;
        }
	if (checkWcuts() && checkZcuts())  {
		acceptance->fill(6);
		vetoEvent;
	}
      } else {if (leptons.size() > 2) {
                setCutsWH(leptons, p_miss, Z, eta_cut_jets);
                if (!checkWcuts()) {
                        acceptance->fill(7);
                        //vetoEvent;
                } else {
			acceptance->fill(8);
			vetoEvent;
		}
              }
      }
        
     


	
      /*
      if (leptons.size() > 3) {
	setCutsWH(leptons, p_miss, Z, eta_cut_jets);
	if (checkWcuts()){
		if (Z_pair_found == 1) {
			setCutsZH(leptons, Z, p_miss, eta_cut_jets, lep_index1, lep_index2);
			if (checkZcuts()) { 
				acceptance->fill(0);
				vetoEvent;
			}
		} else {
			acceptance->fill(1);
			vetoEvent;
		}
      	}
      } else {if (leptons.size() > 2) {
                setCutsWH(leptons, p_miss, Z, eta_cut_jets);
                if (checkWcuts()) {
			acceptance->fill(2);
			vetoEvent;
		}
              } else {
			acceptance->fill(3);
                        vetoEvent;
                }
      }
	
     */

	
      acc_H_PT->fill(higgs.pT()/GeV);
      //acc_V_PT->fill(V.pT()/GeV);                                                                   
      acc_N_jets->fill(jets30.size());
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double sf = m_sumw>0?1.0/m_sumw:1.0;
      for (auto hist:{H_PT, acc_H_PT, N_jets, acc_N_jets}) scale(hist, sf);
    }
    //@}

    
    /// @name Histograms
    //@{
    //Histo1DPtr leptons_PT;
    Histo1DPtr H_PT;
    Histo1DPtr N_jets;
    Histo1DPtr acc_H_PT;
    Histo1DPtr acc_N_jets;
    Histo1DPtr acceptance;
    Histo1DPtr lepton_no;

    Histo1DPtr leading_lep_PT;
    Histo1DPtr trailing_lep_PT;
    //Histo1DPtr leading_lep_eta;
    //Histo1DPtr trailing_lep_eta;
    Histo1DPtr dilep_mass;
    Histo1DPtr dilep_PT;
    //Histo1DPtr trailing_lep_mT_hist;
    //Histo1DPtr higgs_mT_hist;
    Histo1DPtr dphi_hist;
    //Histo1DPtr cos_dphi;
    //Histo1DPtr V_mass;
    //Histo1DPtr V_PT;
    //Histo1DPtr acc_V_PT;

    Histo1DPtr total_lepton_charge;
    Histo1DPtr lep3_PT;
    Histo1DPtr additional_lepton_PT;
    Histo1DPtr X_mass;
    Histo1DPtr max_b_jet_PT;
    Histo1DPtr dilepton_Z_mass_diff;
    Histo1DPtr missing_momentum;
    Histo1DPtr total_lepton_mass;
    Histo1DPtr min_dilep_mass;
    Histo1DPtr total_lep_charge;
    Histo1DPtr leading_jet_PT; 
    Histo1DPtr cut_hist;


    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(VH_WW);


}
