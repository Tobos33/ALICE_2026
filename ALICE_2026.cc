// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Analyses/AliceCommon.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/EventMixingFinalState.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include <cmath>
#include <string>

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2026 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2026);


    /// @name Analysis methods
    /// @{

    double disCPhi(double a1, double a2, double MinBin = -1.408300, double MaxBin = 4.874885) const {
      double diff = a1 - a2;
      while (diff < MinBin)               diff += 2*M_PI;
      while (diff > MaxBin)                diff -= 2*M_PI;
      return diff;
    }

     


    /// Book histograms and initialise projections before the run
    void init() {

    const double etamax = 0.8;
    const double pTmin = 0.5; // GeV
    const double pTmax = 8.0; //GeV
    
      // Primary particles.
    PrimaryParticles pp({Rivet::PID::PIPLUS, Rivet::PID::KPLUS, Rivet::PID::K0S, Rivet::PID::K0L,
                           Rivet::PID::PROTON, Rivet::PID::NEUTRON, Rivet::PID::LAMBDA,
                           Rivet::PID::SIGMAMINUS, Rivet::PID::SIGMAPLUS, Rivet::PID::XIMINUS,
                           Rivet::PID::XI0, Rivet::PID::OMEGAMINUS},
                          Cuts::abseta < etamax && Cuts::pT > pTmin * GeV && Cuts::pT < pTmax * GeV);
    declare(pp,"APRIM");

    ChargedFinalState cfsMult(Cuts::abseta < etamax);//rapidity cut

      // The event mixing projection
    declare(EventMixingFinalState(cfsMult, pp, 10, 0, 100, 10, defaultWeightIndex()),"EVM");

    // Trigger projection.
    declare(ALICE::V0AndTrigger(), "V0-AND");

    // The particle pairs.
    pid = { {2212, 3122}, {2212, -3122}, {2212, 3312}, {2212, -3312}, {3122, 3122}, {3122, -3122}, 
            {3312, 3312}, {3312, -3312}, {2212, 2212}, {2212, -2212}, {2212, 3334}, {2212, -3334},
            {3334, 3334}, {3334, -3334} };
      // The differing pT cuts per pair, in GeV.
    pTcuts = {{0.5, 0.5}, {0.5, 0.5}, {0.5, 0.9}, {0.5, 0.9}, {0.5, 0.5}, {0.5, 0.5},
              {0.9, 0.9}, {0.9, 0.9}, {0.5, 0.5}, {0.5, 0.5}, {0.5, 0.9}, {0.5, 0.9}, {0.9, 0.9}, {0.9, 0.9}}; 
    //The differing pT cuts fo lambdas histograms
    pTcuts_lamb = {{0.5, 1.25}, {1.25, 2.5}, {2.5, 4}, {4, 8}};

    vector<string> Names = { "p_lambda", "p_lambda_bar", "p_Xi", "p_Xi_bar", "lambda_lambda", "lambda_lambda_bar", 
                              "Xi_Xi", "Xi_Xi_bar", "p_p", "p_p_bar", "p-Omega", "p-Omega_bar", "Omega-Omega", "Omega-Omega_bar" }; 

    vector<pair<string, string>> Names_lamb = {{"ll_05_125pT", "l_lbar_05_125pT"}, {"ll_125_250pT", "l_lbar_125_125pT"}, {"ll_25_40pT", "l_lbar_25_40pT"}, {"ll_4_8pT", "l_lbar_4_8pT"}};

    //cout<<"test"<<endl;
    //cout<<"test"<<endl;
    signal.resize(Nhistos);
    background.resize(Nhistos);
    ratioCP.resize(Nhistos);
    nsp.resize(Nhistos);
    nmp.resize(Nhistos);
    signal_lamb.resize(pTcuts_lamb.size());
    background_lamb.resize(pTcuts_lamb.size());
    ratioCP_lamb.resize(pTcuts_lamb.size());
    nsp_lamb.resize(pTcuts_lamb.size());
    nmp_lamb.resize(pTcuts_lamb.size());

    //booking histograms for different partciles pairs
    for(int i = 0 ; i < Nhistos; i++){
      book(ratioCP[i], i+1, 1, 1);
      book(signal[i], "TMP/" + Names[i] + "_signal", refData(1+i, 1, 1));
      book(background[i], "TMP/" + Names[i] + "_background", refData(1+i, 1, 1));
      book(nsp[i], "TMP/" + Names[i] + "_signal_counter");
      book(nmp[i], "TMP/" + Names[i] + "_background_counter");
    }
    //booking histograms for lambda pairs for different pt regions
    for(int i =0; i < (int)Names_lamb.size(); i++){
      //lambda-lambda
      book(ratioCP_lamb[i].first, (int)Names.size()+i*2+1, 1, 1);
      book(signal_lamb[i].first, "TMP/" + Names_lamb[i].first + "_signal", refData(1, 1, 1));
      book(background_lamb[i].first, "TMP/" + Names_lamb[i].first + "_background", refData(1, 1, 1));
      book(nsp_lamb[i].first, "TMP/" + Names_lamb[i].first + "_signal_counter");
      book(nmp_lamb[i].first, "TMP/" + Names_lamb[i].first + "_background_counter");
      
      //lambda-antilambda
      book(ratioCP_lamb[i].second, (int)Names.size()+i*2+2, 1, 1);
      book(signal_lamb[i].second, "TMP/" + Names_lamb[i].second + "_signal", refData(1, 1, 1));
      book(background_lamb[i].second, "TMP/" + Names_lamb[i].second + "_background", refData(1, 1, 1));
      book(nsp_lamb[i].second, "TMP/" + Names_lamb[i].second + "_signal_counter");
      book(nmp_lamb[i].second, "TMP/" + Names_lamb[i].second + "_background_counter");
      
    }
    

    }

    void fillPair(const Particle& p1, const Particle& p2, vector<Histo1DPtr>& histosPhi, vector<CounterPtr>& sow) 
    {
	    if (isSame(p1,p2)) return;
      if (abs(p1.eta() - p2.eta()) > 1.3) return;
      // Figure out which pid pair we are looking at.
      int iPair = -1;
      for (int i = 0, N = pid.size(); i < N; ++i) {
        if (pid[i].first == p1.pid() && pid[i].second == p2.pid()) { //(+ +) and (+ -) pairs
          iPair = i;
          break;
        }
        else if((pid[i].first == pid[i].second || abs(pid[i].first) != abs(pid[i].second)) && -pid[i].first == p1.pid() && -pid[i].second == p2.pid()){ 
          //(- -), (- +) pairs <-- automatically added to (+ +) (+ -)
          iPair = i;
          break;
        }
      }
      // If the pair is not in the analysis, don't fill anything.
      if (iPair < 0) return;
      // Apply min pT cuts, varies for different species.
      if (p1.pT() < pTcuts[iPair].first || p2.pT() < pTcuts[iPair].second || p1.pT() > pTmax_not_lamb || p2.pT() > pTmax_not_lamb) return;

      //fill DeltaPhi histogram
      double d_phi = disCPhi(p1.phi(), p2.phi());
      histosPhi[iPair]->fill(d_phi);

      //fill counter
      sow[iPair]->fill();
    }


    void fillPair_lamb(const Particle& p1, const Particle& p2, vector<pair<Histo1DPtr, Histo1DPtr>>& histosPhi, vector<pair<CounterPtr, CounterPtr>>& sow) 
    {
      if (abs(p1.pid()) != 3122 || abs(p2.pid()) != 3122) return;
	    if (isSame(p1,p2)) return;
      if (abs(p1.eta() - p2.eta()) > 1.3) return;
      // Figure out which pT cut bin we are looking for.
      int iPair_pT = -1;
      
      for (int i = 0, N = pTcuts_lamb.size(); i < N; ++i) {
        //itering for pTcuts for lambda
        if (p1.pT() > pTcuts_lamb[i].first && p1.pT() <= pTcuts_lamb[i].second && p2.pT() > pTcuts_lamb[i].first && p2.pT() <= pTcuts_lamb[i].second) { 
          iPair_pT = i;
          break;
        }
      }
      // If the pair is not in the analysis, don't fill anything.
      if (iPair_pT < 0) return;
      // Check if is it particle-particle or particle-antiparticle pair
      if(p1.pid() == p2.pid()){
        //fill DeltaPhi histogram
        double d_phi = disCPhi(p1.phi(), p2.phi());
        histosPhi[iPair_pT].first->fill(d_phi);
        
        //fill counter
        sow[iPair_pT].first->fill();
      }
      else if(p1.pid() == -p2.pid()){
        //fill DeltaPhi histogram
        double d_phi = disCPhi(p1.phi(), p2.phi());
        histosPhi[iPair_pT].second->fill(d_phi);
        
        //fill counter
        sow[iPair_pT].second->fill();
      }
    
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (!apply<ALICE::V0AndTrigger>(event, "V0-AND")()) return;

      const PrimaryParticles& pp = apply<PrimaryParticles>(event,"APRIM");
      const EventMixingFinalState& evm = apply<EventMixingFinalState>(event, "EVM");

      // Test if we have enough mixing events available to continue.
      if (!evm.hasMixingEvents()) return;


      for (const Particle& p1 : pp.particles()) {
	      // First do the signal histograms.
        for (const Particle& p2 : pp.particles()){
          fillPair(p1, p2, signal, nsp);
          fillPair_lamb(p1, p2, signal_lamb, nsp_lamb);
        }
        // Then do the background
        for (const Particle& p2 : evm.particles()){
          fillPair(p1, p2, background, nmp);
          fillPair_lamb(p1, p2, background_lamb, nmp_lamb);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(int pid = 0; pid < Nhistos; pid++){
        divide(signal[pid], background[pid], ratioCP[pid]);
        scale(ratioCP[pid], nmp[pid]->sumW()/nsp[pid]->sumW());
      }
      for(int pid = 0; pid < (int)pTcuts_lamb.size(); pid++){
        divide(signal_lamb[pid].first, background_lamb[pid].first, ratioCP_lamb[pid].first);
        scale(ratioCP_lamb[pid].first, nmp_lamb[pid].first->sumW()/nsp_lamb[pid].first->sumW());   
        divide(signal_lamb[pid].second, background_lamb[pid].second, ratioCP_lamb[pid].second);
        scale(ratioCP_lamb[pid].second, nmp_lamb[pid].second->sumW()/nsp_lamb[pid].second->sumW());     
      }

    }

    /// @}


    vector<pair<int, int> > pid;
    vector<pair<double, double> > pTcuts;
    vector<pair<double, double> > pTcuts_lamb;
    vector<Histo1DPtr> signal;
    vector<pair<Histo1DPtr, Histo1DPtr>> signal_lamb;
    vector<Histo1DPtr> background;
    vector<pair<Histo1DPtr, Histo1DPtr>> background_lamb;
    vector<Estimate1DPtr> ratioCP;
    vector<pair<Estimate1DPtr, Estimate1DPtr>> ratioCP_lamb;
    vector<CounterPtr> nsp;
    vector<pair<CounterPtr, CounterPtr>> nsp_lamb;
    vector<CounterPtr> nmp;
    vector<pair<CounterPtr, CounterPtr>> nmp_lamb;
    const int Nhistos = 14;
    const double pTmax_not_lamb = 4.0; //GeV


  };


  RIVET_DECLARE_PLUGIN(ALICE_2026);

}
