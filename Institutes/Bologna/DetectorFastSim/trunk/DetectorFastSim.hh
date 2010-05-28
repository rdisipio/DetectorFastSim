#ifndef RIVET_FASTSIM_HH
#define RIVET_FASTSIM_HH


#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"

#include <boost/random.hpp>

#include <map>

#include "fastjet/SISConePlugin.hh"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

#define ECAL_FRAC 0.8
#define HCAL_FRAC 0.25

namespace Rivet {

  struct trigger_chain {
    std::string name;
    double threshold;
  };


  class Trigger {

  public:
    Trigger() {};
    virtual ~Trigger() { m_triggerChains.clear(); };

    void addChain( trigger_chain chain ) { m_triggerChains.push_back( chain ); };

  private:
    std::vector< trigger_chain > m_triggerChains;
    

  };


  class DetectorFastSim : public Analysis {

  public:
    //Default contructor
    DetectorFastSim();
    
    DetectorFastSim::DetectorFastSim(const double B_field, 
				     const double indet_eta_coverage,
				     const double calo_eta_coverage,
				     const double calo_eta_crack,
				     const unsigned int calo_cells_eta,
				     const unsigned int calo_cells_phi,
				     const double calo_em_res_const,
				     const double calo_em_res_sqrtE,
				     const double calo_had_res_sqrtE,
				     const double calo_etmiss_res,
				     const double jet_ptmin, 
				     const double jet_aperture,
				     const FastJets::JetAlgName jet_algorithm  
				       );

    /// Factory method
    static Analysis* create() { 
      return new DetectorFastSim(); 
    }
        
    //@}
        
    /// @name Publication metadata
    //@{
    /// Return the name of this analysis
    string name() const {
      return "DetectorFastSim";
    }
        
    /// Get the SPIRES ID code
    string spiresId() const {
      return "NONE";
    }

    /// Get a description of the analysis.
    string description() const {
      ostringstream os;
      os << "plug-in to simulate a detector";

      return os.str();
    }
        
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "Generic";
    }

    /// Collider
    string collider() const {
      return "LHC";
    }
    
        /// When published (preprint year according to SPIRES).
    string year() const {
      return "NONE";
    }

    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Riccardo Di Sipio <disipio@cern.ch>";
      return rtn;
    }

    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "All event types will be accepted.";
      return os.str();
    }
    string status() const {
      return "EXAMPLE";
    }
    /// No journal or preprint references: this is a demo.
    vector<string> references() const {
      vector<string> ret;
      return ret;
    }

    /// A short description of the analysis.
    string summary() const {
      return "plug-in to simulate a detector";
    }

    //@}
        
    /// @name Analysis methods
    //@{
    void init();
    void clear();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:   //functions

    /// Hide the assignment operator
    DetectorFastSim& operator=(const DetectorFastSim&);

    double rand01();
    double norm01();

    bool inCrack(const double& eta);

    unsigned int get_eta_index(const double& eta);
    unsigned int get_phi_index(const double& phi);
    double get_true_eta(const int& index);
    double get_true_phi(const int& index);

    inline double deltaR(const FourMomentum& p1, const FourMomentum& p2) {
      return sqrt( pow(p1.eta() - p2.eta(),2) + pow(p1.phi() - p2.phi(), 2)  );
    };

    inline double deltaR(const double& eta, const double& phi, const FourMomentum& trk) {
      return sqrt( pow( eta - trk.eta(),2) + pow( phi - trk.phi(), 2)  );
    }

    double getTowerIsolationSum(const double& eta, const double& phi);
    double getTrackPtSum(const double& eta, const double& phi, const double& cone);
    double getCaloSumEnergy(const double& eta, const double& phi, const double& cone);


    //Particle propagateTrack(const FourMomentum& track);
    bool doTruth(const Event& event);
    bool findTracks(const Event& event);
    bool fillCalo(const Event& event);
    bool findClusters();
    bool findMuons(const Event& event);
    bool findElectrons();
    bool doOverlapRemoval();
    bool doTagging();
    bool calculateMissingET();
    bool findJets(const vector<fastjet::PseudoJet>& pseudojets);
    bool applyTrigger();
    bool triggerChain_Leptons();
    bool triggerChain_ETmiss();
    bool triggerChain_Jets();
    
  private:   //member variables

    /// @name Detector parameters
    //@{
    ///Inner Detector
    double m_param_B_field;
    double m_param_indet_eta_coverage;

    ///Calorimeter
    double m_param_calo_eta_coverage; // |eta|<5
    double m_param_calo_eta_crack;
     
    unsigned int m_param_calo_cells_eta;
    unsigned int m_param_calo_cells_phi;
    double m_calo_cells_eta_width;
    double m_calo_cells_phi_width;

    double m_param_calo_em_res_const;
    double m_param_calo_em_res_sqrtE;
    double m_param_calo_had_res_sqrtE;

    double m_param_calo_etmiss_res;

    double m_param_jet_ptmin; //20
    double m_param_calo_jet_aperture; //0.4?
    FastJets::JetAlgName m_param_calo_jet_algorithm;

    ///Trigger
    bool m_param_trigger_jet_prescale;


    //@}


    /// @name Containers
    //@{    
    std::map< std::string, bool > m_trigger_chains;
    bool    m_trig_e10;
    bool    m_trig_e20;
    bool    m_trig_e30;
    bool    m_trig_mu10;
    bool    m_trig_mu20;
    bool    m_trig_mu30;
    bool    m_trig_xe20;
    bool    m_trig_xe30;
    bool    m_trig_xe50;
    bool    m_trig_xe100;
    bool    m_trig_4j20;
    bool    m_trig_2j40;
    bool    m_trig_3j40;
    bool    m_trig_2j100;



    std::vector< FourMomentum >  m_indet_tracks;

    //std::vector< std::vector< double > > m_calo_em;
    //std::vector< std::vector< double > > m_calo_had;
    TH2F       *   m_calo_em;
    TH2F       *   m_calo_had;
    //ParticleVector m_calo_towers;
    std::vector< FourMomentum >  m_calo_tracks;
  

    fastjet::JetDefinition                     m_jet_def;
    shared_ptr<fastjet::ClusterSequence>       m_cluster_seq;
    shared_ptr<fastjet::JetDefinition::Plugin> m_jet_plugin;

    ParticleVector m_mu_tracks;
    Trigger        m_trigger;
  
    //ParticleVector m_fs_ele;
    std::vector< FourMomentum > m_fs_ele;
    int              m_fs_ele_n;
    vector<double> * m_fs_ele_pt;
    vector<double> * m_fs_ele_eta;
    vector<double> * m_fs_ele_phi;
    vector<double> * m_fs_ele_iso;
    vector<int>   * m_fs_ele_matched;
    
    ParticleVector   m_truth_ele;
    int              m_truth_ele_n;
    vector<double> * m_truth_ele_pt;
    vector<double> * m_truth_ele_eta;
    vector<double> * m_truth_ele_phi;

    ParticleVector m_fs_gamma;

    /** muons */
    std::vector< FourMomentum > m_fs_mu;
    int              m_fs_mu_n;
    vector<double> * m_fs_mu_pt;
    vector<double> * m_fs_mu_eta;
    vector<double> * m_fs_mu_phi;
    vector<double> * m_fs_mu_iso;
    vector<int>   * m_fs_mu_matched;

    ParticleVector   m_truth_mu;
    int              m_truth_mu_n;
    vector<double> * m_truth_mu_pt;
    vector<double> * m_truth_mu_eta;
    vector<double> * m_truth_mu_phi;

    
    /** taus  */
    Jets            m_fs_tjets;
    int             m_fs_tjets_n;

    ParticleVector m_truth_taus;

    
    /** jets */
   
    ParticleVector m_truth_bquarks;


    //jets before overlap removal
    Jets               m_raw_jets;
    int                m_raw_jets_n;
    vector<double> *   m_raw_jets_pt;
    vector<double> *   m_raw_jets_eta;
    vector<double> *   m_raw_jets_phi;

    //jets after overlap removal
    Jets               m_fs_jets;
    int                m_fs_jets_n;
    vector<double> *   m_fs_jets_pt;
    vector<double> *   m_fs_jets_eta;
    vector<double> *   m_fs_jets_phi;
    vector<int>   *   m_fs_jets_matched;

    int                m_fs_bjets_n;
    vector<int>    *   m_fs_jets_tag;
    vector<double> *   m_fs_jets_tagweight;

    Jets               m_fs_bjets;

    Jets               m_truth_jets;
    int                m_truth_jets_n;
    vector<double> *   m_truth_jets_pt;
    vector<double> *   m_truth_jets_eta;
    vector<double> *   m_truth_jets_phi;



    double         m_fs_calo_met;
    double         m_fs_etmiss_px;
    double         m_fs_etmiss_py;
    double         m_fs_etmiss_pt;
    double         m_fs_etmiss_phi;
    double         m_fs_sumet;

    double         m_truth_etmiss_px;
    double         m_truth_etmiss_py;
    double         m_truth_etmiss_pt;
    double         m_truth_etmiss_phi;
    double         m_truth_sumet;

    
    

    //@}


    TTree * m_tree;
    TFile * m_outfile;
  };


  extern "C" {
    AnalysisBuilders getAnalysisBuilders() {
      AnalysisBuilders fns;
      fns["FASTSIM"] = Rivet::DetectorFastSim::create;
      return fns;
    }
  }
   

}; //end namespace Rivet



#endif /** RIVET_FASTSIM_HH */
