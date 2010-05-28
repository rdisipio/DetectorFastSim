// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Math/MathUtils.hh"

#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/InvMassFinalState.hh"

#include "HepMC/Units.h"

#include "DetectorFastSim.hh"


namespace Rivet {

  DetectorFastSim::DetectorFastSim() :
    m_param_B_field(4.0),
    m_param_indet_eta_coverage(2.5),
    m_param_calo_eta_coverage(5.0),
    m_param_calo_eta_crack(1.2),
    m_param_calo_cells_eta(500),
    m_param_calo_cells_phi(500),
    m_param_calo_em_res_const(0.0044),
    m_param_calo_em_res_sqrtE(0.024),
    m_param_calo_had_res_sqrtE(0.8),
    m_param_calo_etmiss_res(0.2),
    m_param_jet_ptmin(20.),
    m_param_calo_jet_aperture(0.4),
    m_param_calo_jet_algorithm(FastJets::SISCONE)
  {
    
    setBeams(PROTON, PROTON);
    setNeedsCrossSection(false);
  }


  DetectorFastSim::DetectorFastSim(const double B_field            = 4.0, 
				   const double indet_eta_coverage = 2.5,
				   const double calo_eta_coverage  = 5.0,
				   const double calo_eta_crack     = 1.2,
				   const unsigned int calo_cells_eta     = 500,
				   const unsigned int calo_cells_phi     = 500,
				   const double calo_em_res_const  = 0.0044,
				   const double calo_em_res_sqrtE  = 0.024,
				   const double calo_had_res_sqrtE = 0.8,
				   const double calo_etmiss_res    = 0.2,
				   const double jet_ptmin          = 20., 
				   const double jet_aperture       = 0.4,
				   const FastJets::JetAlgName jet_algorithm = FastJets::SISCONE
				   ) :
    m_param_B_field(B_field),
    m_param_indet_eta_coverage(indet_eta_coverage),
    m_param_calo_eta_coverage(calo_eta_coverage),
    m_param_calo_eta_crack(calo_eta_crack),
    m_param_calo_cells_eta(calo_cells_eta),
    m_param_calo_cells_phi(calo_cells_phi),
    m_param_calo_em_res_const(calo_em_res_const),
    m_param_calo_em_res_sqrtE(calo_em_res_sqrtE),
    m_param_calo_had_res_sqrtE(calo_had_res_sqrtE),
    m_param_calo_etmiss_res(calo_etmiss_res),
    m_param_jet_ptmin(jet_ptmin),
    m_param_calo_jet_aperture(jet_aperture),
    m_param_calo_jet_algorithm(jet_algorithm)
  {
    setBeams(PROTON, PROTON);
    setNeedsCrossSection(false);
  }
  

  //////////////////////////////////////////////////////


  void DetectorFastSim::init()
  {
    Log log = getLog();
    log << Log::DEBUG << "Initializing detector" << endl;


    ///Open output ROOT file
    string outname = tree().storeName();
    outname = outname.substr( 0, outname.size() - 5 ); //get rid of .aida
    outname = outname + ".root";
    m_outfile = new TFile(outname.c_str(), "RECREATE");

    log << Log::INFO << "Outfile " << outname << " opened" << endl;
    
    m_tree = new TTree("events", "Events");

    /** ele */
    m_fs_ele_pt  = new vector<double>;
    m_fs_ele_eta = new vector<double>;
    m_fs_ele_phi = new vector<double>;
    m_fs_ele_iso = new vector<double>;
    m_fs_ele_matched = new vector<int>;

    m_truth_ele_pt  = new vector<double>;
    m_truth_ele_eta = new vector<double>;
    m_truth_ele_phi = new vector<double>;

    /** muons */
    m_fs_mu_pt  = new vector<double>;
    m_fs_mu_eta = new vector<double>;
    m_fs_mu_phi = new vector<double>;
    m_fs_mu_iso = new vector<double>;
    m_fs_mu_matched = new vector<int>;

    m_truth_mu_pt  = new vector<double>;
    m_truth_mu_eta = new vector<double>;
    m_truth_mu_phi = new vector<double>;

    /** taus */

    /** jets */
    m_raw_jets_pt  = new vector<double>;
    m_raw_jets_eta = new vector<double>;
    m_raw_jets_phi = new vector<double>;

    m_fs_jets_pt  = new vector<double>;
    m_fs_jets_eta = new vector<double>;
    m_fs_jets_phi = new vector<double>;
    m_fs_jets_matched = new vector<int>;

    m_fs_jets_tag = new vector<int>;
    m_fs_jets_tagweight = new vector<double>;

    m_truth_jets_pt  = new vector<double>;
    m_truth_jets_eta = new vector<double>;
    m_truth_jets_phi = new vector<double>;

   
    // Add tree branches
    m_tree->Branch("reco_raw_jets_n", &m_raw_jets_n,    "reco_raw_jets_n/i");
    m_tree->Branch("reco_raw_jets_pt", m_raw_jets_pt);
    m_tree->Branch("reco_raw_jets_eta", m_raw_jets_eta);
    m_tree->Branch("reco_raw_jets_phi", m_raw_jets_phi);

    m_tree->Branch("reco_jets_n",  &m_fs_jets_n,    "reco_jets_n/i");
    m_tree->Branch("reco_jets_pt",  m_fs_jets_pt);
    m_tree->Branch("reco_jets_eta", m_fs_jets_eta);
    m_tree->Branch("reco_jets_phi", m_fs_jets_phi);
    m_tree->Branch("reco_jets_matched", m_fs_jets_matched);

    m_tree->Branch("reco_jets_tag", m_fs_jets_tag);
    m_tree->Branch("reco_jets_tagweight", m_fs_jets_tagweight);
    m_tree->Branch("reco_bjets_n", &m_fs_bjets_n,   "reco_bjets_n/i");

    m_tree->Branch("reco_tjets_n", &m_fs_tjets_n,   "reco_tjets_n/i");


    m_tree->Branch("truth_jets_n",  &m_truth_jets_n,    "truth_jets_n/i");
    m_tree->Branch("truth_jets_pt",  m_truth_jets_pt);
    m_tree->Branch("truth_jets_eta", m_truth_jets_eta);
    m_tree->Branch("truth_jets_phi", m_truth_jets_phi);


    /** ele */

    m_tree->Branch("reco_ele_n",      &m_fs_ele_n, "reco_ele_n/i");
    m_tree->Branch("reco_ele_pt",      m_fs_ele_pt);
    m_tree->Branch("reco_ele_eta",     m_fs_ele_eta);
    m_tree->Branch("reco_ele_phi",     m_fs_ele_phi);
    m_tree->Branch("reco_ele_iso",     m_fs_ele_iso);
    m_tree->Branch("reco_ele_matched", m_fs_ele_matched);

    m_tree->Branch("truth_ele_n",      &m_truth_ele_n, "truth_ele_n/i");
    m_tree->Branch("truth_ele_pt",      m_truth_ele_pt);
    m_tree->Branch("truth_ele_eta",     m_truth_ele_eta);
    m_tree->Branch("truth_ele_phi",     m_truth_ele_phi);

    /** muon */
    m_tree->Branch("reco_mu_n",      &m_fs_mu_n, "reco_mu_n/i");
    m_tree->Branch("reco_mu_pt",      m_fs_mu_pt);
    m_tree->Branch("reco_mu_eta",     m_fs_mu_eta);
    m_tree->Branch("reco_mu_phi",     m_fs_mu_phi);
    m_tree->Branch("reco_mu_iso",     m_fs_mu_iso);
    m_tree->Branch("reco_mu_matched", m_fs_mu_matched);

    m_tree->Branch("truth_mu_n",      &m_truth_mu_n, "truth_mu_n/i");
    m_tree->Branch("truth_mu_pt",      m_truth_mu_pt);
    m_tree->Branch("truth_mu_eta",     m_truth_mu_eta);
    m_tree->Branch("truth_mu_phi",     m_truth_mu_phi);

    /** MissingEt */
    m_tree->Branch("reco_etmiss_px",  &m_fs_etmiss_px,  "reco_etmiss_px/d");
    m_tree->Branch("reco_etmiss_py",  &m_fs_etmiss_py,  "reco_etmiss_py/d");
    m_tree->Branch("reco_calo_met",   &m_fs_calo_met,   "reco_calo_met/d");
    m_tree->Branch("reco_etmiss_pt",  &m_fs_etmiss_pt,  "reco_etmiss_pt/d");
    m_tree->Branch("reco_etmiss_phi", &m_fs_etmiss_phi, "reco_etmiss_phi/d");
    m_tree->Branch("reco_sumet",      &m_fs_sumet,      "reco_sumet/d");
    
    m_tree->Branch("truth_etmiss_px",  &m_truth_etmiss_px,  "truth_etmiss_px/d");
    m_tree->Branch("truth_etmiss_py",  &m_truth_etmiss_py,  "truth_etmiss_py/d");
    m_tree->Branch("truth_etmiss_pt",  &m_truth_etmiss_pt,  "truth_etmiss_pt/d");
    m_tree->Branch("truth_etmiss_phi", &m_truth_etmiss_phi, "truth_etmiss_phi/d");
    m_tree->Branch("truth_sumet",      &m_truth_sumet,      "truth_sumet/d");

    /** trigger */
    m_tree->Branch("trig_e10",    &m_trig_e10,      "trig_e10/b");
    m_tree->Branch("trig_e20",    &m_trig_e20,      "trig_e20/b");
    m_tree->Branch("trig_e30",    &m_trig_e30,      "trig_e30/b");
    m_tree->Branch("trig_mu10",   &m_trig_mu10,     "trig_e10/b");
    m_tree->Branch("trig_mu20",   &m_trig_mu20,     "trig_mu20/b");
    m_tree->Branch("trig_mu30",   &m_trig_mu30,     "trig_mu30/b");

    m_tree->Branch("trig_xe20",   &m_trig_xe20,     "trig_xe20/b");
    m_tree->Branch("trig_xe30",   &m_trig_xe30,     "trig_xe30/b");
    m_tree->Branch("trig_xe50",   &m_trig_xe50,     "trig_xe50/b");
    m_tree->Branch("trig_xe100",  &m_trig_xe100,    "trig_xe100/b");

    m_tree->Branch("trig_4j20",   &m_trig_4j20,     "trig_4j20/b");
    m_tree->Branch("trig_2j40",   &m_trig_2j40,     "trig_2j40/b");
    m_tree->Branch("trig_3j40",   &m_trig_3j40,     "trig_3j40/b");
    m_tree->Branch("trig_2j100",  &m_trig_2j100,    "trig_2j100/b");


    //Init calorimeter cells
    //Use eta X phi space
    m_calo_cells_eta_width = m_param_calo_eta_coverage / m_param_calo_cells_eta;
    m_calo_cells_phi_width = 2 * 3.1415 / m_param_calo_cells_phi;
    
    log << Log::INFO << "Calorimeter cells (" << m_calo_cells_eta_width 
	<< " X " << m_calo_cells_phi_width << ")" << endl;
    /*
      for( unsigned int eta = 0 ; eta < 2 * m_param_calo_cells_eta ; ++eta ) {
      std::vector< double > eta_row( m_param_calo_cells_phi, 0.0 );
      m_calo_em.push_back(  eta_row );
      m_calo_had.push_back( eta_row );
      }
    */
    m_calo_em  = new TH2F("calo_em", "EM Calorimeter",  
			  2*m_param_calo_cells_eta, -m_param_calo_eta_coverage,  m_param_calo_eta_coverage,
			  m_param_calo_cells_phi, 0.0, 2*3.1415 );

    m_calo_had = new TH2F("calo_had", "Had Calorimeter",  
			  2*m_param_calo_cells_eta, -m_param_calo_eta_coverage,  m_param_calo_eta_coverage,
			  m_param_calo_cells_phi, 0.0, 2*3.1415 );
    
    /// build final states

    FinalState fs(-5, 5, 0.5);
    addProjection(fs, "FS");
    
    FinalState fs_barrel(-2.5, 2.5);
    addProjection(fs_barrel, "FSBARREL");

    
    //charged f.s. particles
    const ChargedFinalState cfs(fs);
    addProjection(cfs, "ChFS");

    const ChargedLeptons chlep(cfs);
    addProjection(chlep, "ChLEP");

    //no neutrinos
    //VetoedFinalState vfs_visible( FinalState(-5, 5, 0.5) );
    VetoedFinalState vfs_visible( fs );
    vfs_visible
      .addVetoPairId(NU_E)
      .addVetoPairId(NU_MU)
      .addVetoPairId(NU_TAU)
      .addVetoId(1000022); // LSP
    addProjection(vfs_visible, "VFSVISIBLE");

    //make a fs with no W decay products - for jet clustering
    std::vector<std::pair<long,long> > vids;
    vids.push_back(make_pair(ELECTRON, NU_EBAR) );
    vids.push_back(make_pair(POSITRON, NU_E)    );
    vids.push_back(make_pair(MUON, NU_MUBAR)    );
    vids.push_back(make_pair(ANTIMUON, NU_MU)   );
    
    InvMassFinalState invfs(fs_barrel, vids, 65*GeV, 95*GeV);
    addProjection(invfs, "INVFS");
    
    //VetoedFinalState vfs_noW(fs_barrel);
    VetoedFinalState vfs_noW(fs);
    vfs_noW.addVetoOnThisFinalState(invfs);
    addProjection(vfs_noW, "VFSNOW");

    /** set jet definition */
    if ( m_param_calo_jet_algorithm == FastJets::KT) {
      m_jet_def = fastjet::JetDefinition(fastjet::kt_algorithm, m_param_calo_jet_aperture, fastjet::E_scheme);
    }
    else if ( m_param_calo_jet_algorithm == FastJets::SISCONE) {
      const double OVERLAP_THRESHOLD = 0.75;
      m_jet_plugin.reset(new fastjet::SISConePlugin(m_param_calo_jet_aperture, OVERLAP_THRESHOLD));
      m_jet_def = fastjet::JetDefinition( m_jet_plugin.get() );
    } 


    // Trigger stuff
    m_param_trigger_jet_prescale = 10; // fire 1 event every N

  }


  //////////////////////////////////////////////////////


  double DetectorFastSim::rand01()
  {
    boost::mt19937 rng(43);
    static boost::uniform_01<boost::mt19937> zeroone(rng);
    return zeroone();
  }


  double DetectorFastSim::norm01()
  {
    boost::lagged_fibonacci19937 engine;
    static boost::normal_distribution<double> norm_dist(0.0, 1.0);
    return norm_dist.operator() <boost::lagged_fibonacci19937>((engine));

  }
  
  bool DetectorFastSim::inCrack(const double& eta)
  { 
    return ( fabs( fabs(eta) - m_param_calo_eta_crack ) < 0.1 ); 
  }


  unsigned int DetectorFastSim::get_eta_index(const double& eta)
  {
    //double eta_shifted = eta + m_param_calo_eta_coverage ;
    //return (int)( 1 + eta_shifted /  m_calo_cells_eta_width );
    //static double eta_min = -(m_param_calo_cells_eta / 2) *  m_calo_cells_eta_width;
    //static double eta_max = -eta_min;

    unsigned int eta_index = 0;
    //if( eta > eta_min && eta < eta_max ) {
    eta_index = (unsigned int)( (eta + m_param_calo_eta_coverage )/ m_calo_cells_eta_width + 0.5 );

    if( eta_index > 2*m_param_calo_cells_eta ) eta_index = 1;
    else if( eta_index < 0 ) eta_index = 1;
    //}

    return eta_index;
  }
  

  unsigned int DetectorFastSim::get_phi_index(const double& phi)
  {
    // return (int)( 1 + phi / m_calo_cells_phi_width );
    unsigned int phi_index = (unsigned int)( phi / m_calo_cells_phi_width + 0.5 );
    if( phi_index > m_param_calo_cells_phi ) phi_index -= m_param_calo_cells_phi;
    else if( phi_index < 0 ) phi_index += m_param_calo_cells_phi;

    return phi_index;
  }

  double DetectorFastSim::get_true_eta(const int& index)
  {
    //static const double WEIRD_DIFF = 0.000000000000001;
    // static double eta_min = - ( m_param_calo_cells_eta / 2 ) * m_calo_cells_eta_width;

    //return eta_min + m_calo_cells_eta_width/2. + (double)(index-1)*m_calo_cells_eta_width;
    double eta =  ( (index * m_calo_cells_eta_width) - m_param_calo_eta_coverage);
    // if( fabs(eta) < WEIRD_DIFF || fabs(eta - m_param_calo_eta_coverage) < WEIRD_DIFF ) eta = 0.0;

    return eta;
  }

  double DetectorFastSim::get_true_phi(const int& index)
  {
    return m_calo_cells_phi_width/2. + (double)( (index-1)*m_calo_cells_phi_width);
    //return ( index * m_calo_cells_phi_width );
  }


  double DetectorFastSim::getTowerIsolationSum(const double& eta, const double& phi)
  {
    Log log = getLog();

    //return energy in 3x3 cells around given (eta,phi), except seed energy
    const unsigned int seed_eta_bin = get_eta_index(eta);   
    const unsigned int seed_phi_bin = get_phi_index(phi);

    //log << Log::DEBUG << "seed: (eta,phi)=(" << eta << ", " << phi << ")=(" 
    //    << seed_eta_bin << ", " << seed_phi_bin << ") " << endl;

    unsigned int phi_left = seed_phi_bin > 1 ? seed_phi_bin - 1 : m_param_calo_cells_phi;
    //if( phi_left < 1 ) phi_left = m_param_calo_cells_phi;

    unsigned int phi_right = seed_phi_bin < m_param_calo_cells_phi ? seed_phi_bin + 1 : 1;
    //if( phi_right > m_param_calo_cells_phi ) phi_right = 1;

    //look up true eta
    const double true_eta_up   = get_true_eta( seed_eta_bin + 1 );
    const double true_eta_down = get_true_eta( seed_eta_bin - 1 );

    const double theta_up   = 2 * atan( exp(-1*true_eta_up) );
    const double theta_seed = 2 * atan( exp(-1*eta) );
    const double theta_down = 2 * atan( exp(-1*true_eta_down) );
    
    const double s_up   = sin(theta_up);
    const double s_seed = sin(theta_seed);
    const double s_down = sin(theta_down);
    
    //upper row
    const double Eul = m_calo_em->GetBinContent(seed_eta_bin + 1, phi_left) 
      + m_calo_had->GetBinContent(seed_eta_bin + 1, phi_left);

    const double Euc =  m_calo_em->GetBinContent(seed_eta_bin + 1, seed_phi_bin) 
      + m_calo_had->GetBinContent(seed_eta_bin + 1, seed_phi_bin);

    const double Eur = m_calo_em->GetBinContent(seed_eta_bin + 1, phi_right) 
      + m_calo_had->GetBinContent(seed_eta_bin + 1, phi_right);

    //center row - SKIP SEED BIN!
    const double Ecl =  m_calo_em->GetBinContent(seed_eta_bin, phi_left) 
      + m_calo_had->GetBinContent(seed_eta_bin, phi_left);

    const double Ecc = 0.0;
    //m_calo_em->GetBinContent(seed_eta_bin, seed_phi_bin) 
    //  + m_calo_had->GetBinContent(seed_eta_bin, seed_phi_bin);
    //log << Log::DEBUG << " * seed E=" << Ecc << " ET=" << Ecc * s_seed << endl;

    const double Ecr = m_calo_em->GetBinContent(seed_eta_bin , phi_right) 
      + m_calo_had->GetBinContent(seed_eta_bin, phi_right);

    //lower row
    const double Edl = m_calo_em->GetBinContent(seed_eta_bin - 1, phi_left) 
      + m_calo_had->GetBinContent(seed_eta_bin - 1, phi_left);

    const double Edc =  m_calo_em->GetBinContent(seed_eta_bin - 1, seed_phi_bin) 
      + m_calo_had->GetBinContent(seed_eta_bin - 1, seed_phi_bin);

    const double Edr = m_calo_em->GetBinContent(seed_eta_bin - 1, phi_right) 
      + m_calo_had->GetBinContent(seed_eta_bin - 1, phi_right);

    const double sum = s_up * (  Eul + Euc + Eur )
      + s_seed * ( Ecl + Ecc + Ecr )
      + s_down * ( Edl + Edc + Edr );
	       
    return sum;
  }



  double DetectorFastSim::getTrackPtSum(const double& eta, const double& phi, const double& cone)
  {
    double ptsum     = 0.0;
    int    trk_found = 0;

    foreach(const FourMomentum& track, m_indet_tracks) {
      if( track.pT() < 0.5 ) continue;

      const double dR = sqrt( pow( eta - track.eta(),2 ) + pow( phi - track.phi(),2 ) );
      if( dR > cone ) continue;
      
      ++trk_found;
      ptsum += track.pT();
    }
   
    //if( trk_found == 0 ) ptsum = 100000000;
    return ptsum;
  }

  
  double DetectorFastSim::getCaloSumEnergy(const double& eta, const double& phi, const double& cone)
  {
    double SumE = 0.0;

    foreach( const FourMomentum& tower, m_calo_tracks ) {
      const double dR = deltaR( eta, phi, tower );
      
      if( dR > cone ) continue;

      SumE += tower.E();
    }

    return SumE;
  }

  //////////////////////////////////////////////////////
  

  void DetectorFastSim::clear() {
    //reset tracker
    m_indet_tracks.clear();
    m_calo_tracks.clear();

    //reset calorimeters
    /*for( unsigned int eta = 0 ; eta < 2 * m_param_calo_cells_eta ; ++eta ) {
      m_calo_em[eta].clear();
      m_calo_had[eta].clear();
      }*/
    // m_calo_em->Clear();
    //m_calo_had->Clear();
    m_calo_em->Reset();
    m_calo_had->Reset();

    //electrons
    
    m_fs_ele.clear();
    m_fs_ele_n = 0;
    m_fs_ele_pt->clear();
    m_fs_ele_eta->clear();
    m_fs_ele_phi->clear();
    m_fs_ele_iso->clear();
    m_fs_ele_matched->clear();

    m_truth_ele.clear();
    m_truth_ele_n = 0;
    m_truth_ele_pt->clear();
    m_truth_ele_eta->clear();
    m_truth_ele_phi->clear();

    //gamma
    m_fs_gamma.clear();

    //muons
    m_fs_mu.clear();
    m_fs_mu_n = 0;
    m_fs_mu_pt->clear();
    m_fs_mu_eta->clear();
    m_fs_mu_phi->clear();
    m_fs_mu_iso->clear();
    m_fs_mu_matched->clear();

    m_truth_mu.clear();
    m_truth_mu_n = 0;
    m_truth_mu_pt->clear();
    m_truth_mu_eta->clear();
    m_truth_mu_phi->clear();

    //jets
    m_truth_bquarks.clear();
    m_truth_taus.clear();

    m_raw_jets.clear();
    m_raw_jets_n = 0;
    m_raw_jets_pt->clear();    
    m_raw_jets_eta->clear();
    m_raw_jets_phi->clear();

    m_fs_jets.clear();
    m_fs_jets_n = 0;
    m_fs_jets_pt->clear();
    m_fs_jets_eta->clear();
    m_fs_jets_phi->clear();
    m_fs_jets_matched->clear();

    m_fs_jets_tag->clear();
    m_fs_jets_tagweight->clear();
    m_fs_bjets_n = 0;
    m_fs_tjets_n = 0;

    m_fs_bjets.clear();
    m_fs_tjets.clear();

    m_truth_jets.clear();
    m_truth_jets_n = 0;
    m_truth_jets_pt->clear();    
    m_truth_jets_eta->clear();
    m_truth_jets_phi->clear();


    m_fs_etmiss_px  = 0.0;
    m_fs_etmiss_py  = 0.0;
    m_fs_calo_met   = 0.0;
    m_fs_etmiss_pt  = 0.0;
    m_fs_etmiss_phi = 0.0;
    m_fs_sumet      = 0.0;

    m_truth_etmiss_px  = 0.0;
    m_truth_etmiss_py  = 0.0;
    m_truth_etmiss_pt  = 0.0;
    m_truth_etmiss_phi = 0.0;
    m_truth_sumet      = 0.0;

    m_trigger_chains.clear();
    m_trig_e10   = false;
    m_trig_e20   = false;
    m_trig_e30   = false;
    m_trig_mu10  = false;
    m_trig_mu20  = false;
    m_trig_mu30  = false;
    m_trig_xe20  = false;
    m_trig_xe30  = false;
    m_trig_xe50  = false;
    m_trig_xe100 = false;
    m_trig_4j20  = false;
    m_trig_2j40  = false;
    m_trig_3j40  = false;
    m_trig_2j100 = false;
    
  }


  //////////////////////////////////////////////////////


  void DetectorFastSim::analyze(const Event& event) {
    Log log = getLog();

    clear();

    bool status = false;

    //const FinalState& fs  = applyProjection<FinalState>(event, "FS");
    const VetoedFinalState& fs = applyProjection<VetoedFinalState>(event, "VFSVISIBLE");
    if (fs.isEmpty()) {
      log << Log::INFO << "Empty event!" << endl;
      return; // skip if empty!
    }

    status = doTruth(event);
    if(!status) {
      log << Log::INFO << "Impossible to find truth particles" << endl;
      return;
    }

    status = findTracks(event);
    if(!status) {
      log << Log::INFO << "Impossible to find inner detector tracks" << endl;
      return;
    }

    status = fillCalo(event);
    if(!status) {
      log << Log::INFO << "Impossible to fill calorimeter cells" << endl;
      return;
    }
    
    status = findClusters();
    if(!status) {
      log << Log::INFO << "Impossible to find clusters" << endl;
      return;
    }

    status = findElectrons();
    if(!status) {
      log << Log::INFO << "Impossible to find electrons" << endl;
      return;
    }

    status = findMuons(event);
    if(!status) {
      log << Log::INFO << "Impossible to find muons" << endl;
      return;
    }

    status = doOverlapRemoval();
    if(!status) {
      log << Log::INFO << "Impossible to apply overlap removal" << endl;
      return;
    }

    status = doTagging();
    if(!status) {
      log << Log::INFO << "Impossible to tag jets" << endl;
      return;
    }

    status = calculateMissingET();
    if(!status) {
      log << Log::INFO << "Impossible to calculate ETmiss" << endl;
      return;
    }

    status = applyTrigger();
    if(!status) {
      log << Log::INFO << "Impossible to apply trigger" << endl;
      return;
    }

    m_tree->Fill();
    return;
  }

  //////////////////////////////////////

  void DetectorFastSim::finalize() {
    Log log = getLog();

    log << Log::INFO << "Finalizing simulation" << endl;

    /////////////////////////////////////
    // Print Stats                     //
    /////////////////////////////////////

    //foreach(FourMomentum& track, m_indet_tracks) {
    //  log << Log::INFO << "Q = " << track.t() << " pT = " << track.perp() << endl;
    // }

    m_outfile->Write();
    m_outfile->Close();

    //delete m_outfile;
    //delete m_tree;
  }


  //////////////////////////////////////


  bool DetectorFastSim::doTruth(const Event& event)
  {
    Log log = getLog();
    log << Log::DEBUG << "Finding Truth Objets" << endl; 

    bool status = true;

    //b quarks
    for (GenEvent::particle_const_iterator p = event.genEvent().particles_begin();
         p != event.genEvent().particles_end(); ++p) {
      //++np;
      //if( np > 100 ) break;

      const GenVertex* pv = (*p)->production_vertex();
      //const GenVertex* dv = (*p)->end_vertex();
      
      const int status = (*p)->status();
      const int apid   = abs( (*p)->pdg_id() );

      if( apid == 5 || apid == 15 ) {

	//avoid double counting due to self-decay
	bool selfdecay = false;
	if(pv) {
	  for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin() ;
	       pp != pv->particles_in_const_end() ; ++pp) {
	    // Avoid double counting if particle == parent
	    if ( apid == (*pp)->pdg_id() )
	      selfdecay = true;
	  }	 
	}
	if(!selfdecay) {
	  if( apid == 5 )
	    m_truth_bquarks.push_back( Particle(*(*p)) );
	  else if( apid == 15 ) 
	    m_truth_taus.push_back( Particle(*(*p)) );
	}
      }
      else if(apid == 12 || apid == 14 || apid == 16) {
	// Neutrinos = Missing ET source

	if( status != 1 ) continue;
	const double theta = (*p)->momentum().theta();
	m_truth_sumet     += (*p)->momentum().e() * cos( theta );
	m_truth_etmiss_px += (*p)->momentum().px();
	m_truth_etmiss_py += (*p)->momentum().py();

      }
    }
    log << Log::DEBUG << " * Truth: b quarks N = "    << m_truth_bquarks.size() << endl;
    log << Log::DEBUG << " * Truth: tau leptons N = " << m_truth_taus.size() << endl;

    m_truth_etmiss_pt  = sqrt( pow(m_truth_etmiss_px,2) + pow(m_truth_etmiss_py,2) );
    m_truth_etmiss_phi = atan2( m_truth_etmiss_px, m_truth_etmiss_py );

    log << Log::DEBUG << " * Truth: ETmiss = "       << m_truth_etmiss_pt << endl;

    //Charged Leptons
    const ChargedLeptons chlep =  applyProjection<ChargedLeptons>(event, "ChLEP");
    const ParticleVector allLeptons = chlep.chargedLeptons();
    std::string l_name = "";

    foreach(const Particle& lept, allLeptons ) {
      if( lept.momentum().pT() > 10.*GeV ) {

	const double eta = lept.momentum().eta();
	const double phi = lept.momentum().phi();
	const double pt  = lept.momentum().pT();

	const int apid = abs( lept.pdgId() );

	if( apid == 11 )  {
	  m_truth_ele.push_back( lept );
	  m_truth_ele_pt->push_back( pt );
	  m_truth_ele_eta->push_back( eta );
	  m_truth_ele_phi->push_back( phi );

	  l_name = "ele";
	}
	else if( apid == 13 ) {
	  m_truth_mu.push_back( lept );
	  m_truth_mu_pt->push_back( pt );
	  m_truth_mu_eta->push_back( eta );
	  m_truth_mu_phi->push_back( phi );
	  l_name = "mu";
	}
	log << Log::DEBUG << " * Truth: "<< l_name << ": (" << eta << ", " << phi << ") pT=" << pt << endl;    
      }
    }
    
    //Jets
    const VetoedFinalState& vfs = applyProjection<VetoedFinalState>(event, "VFSNOW");

    addProjection( FastJets(vfs, m_param_calo_jet_algorithm, m_param_calo_jet_aperture), "Jets" );

    const FastJets& fastjets = applyProjection<FastJets>(event, "Jets");
    const Jets all_jets      = fastjets.jetsByPt( m_param_jet_ptmin ); //20.GeV
    
    foreach(const Jet& j, all_jets) {
      const double eta = j.momentum().eta();
      const double phi = j.momentum().phi();
      const double pt  = j.momentum().pT();

      bool overlap = false;
      foreach(const Particle& ele, m_truth_ele) {
	const double ele_eta = ele.momentum().eta();
	const double ele_phi = ele.momentum().phi();

	const double dR = sqrt( pow( eta - ele_eta, 2) + pow( phi - ele_phi, 2 ) );
	if( dR < 0.3 ) {
	  log << Log::DEBUG << "e/j overlap" << endl;
	  overlap = true;
	  break;
	}
      }

      if(!overlap) {
	log << Log::DEBUG << " * Truth: jet: (" << eta << ", " << phi << ") pT=" << pt << endl; 
	m_truth_jets.push_back( j );
	m_truth_jets_pt->push_back( pt );
	m_truth_jets_eta->push_back( eta );
	m_truth_jets_phi->push_back( phi );
      }
    }
    log << Log::DEBUG << " * Truth : found " << m_truth_jets.size() << " jets" << endl;
    
    return status;
  }

  //////////////////////////////////////


  bool DetectorFastSim::findTracks(const Event& event)
  {
    Log log = getLog();
    log << Log::DEBUG << "Finding Tracks" << endl; 

    bool status = true;

    const ChargedFinalState& fs = applyProjection<ChargedFinalState>(event, "ChFS");
    if (fs.isEmpty()) {
      return false; // skip if empty!
    }

    unsigned int all_tracks    = 0;
    unsigned int missed_tracks = 0;

    ParticleVector particles = fs.entities();
    foreach(const Particle& particle, particles) {
      double trk_pt =  particle.momentum().pT();
     
      if( trk_pt < 1.0 ) continue;

      ++all_tracks;
      //apply efficiency. 
      if( rand01() > 0.98 ) {
	++missed_tracks;
	continue;
      }

      const double tantheta = tan( particle.momentum().theta() );
      
      double charge = PID::threeCharge( particle.pdgId() ) / 3.0;

      //smear momentum;

      //calculate radius of curvature
      double radius = trk_pt / (0.3 * m_param_B_field );
      
      //remove tracks that curl up
      if( radius <= (m_param_indet_eta_coverage / 2.0 ) + 0.01 ) continue;

      //calculate sagitta
      const double sinphi = (m_param_indet_eta_coverage / 2.0) / radius;
      double sagitta = radius - radius * sqrt( 1.0 - sinphi*sinphi );
	
      //smear sagitta
      double sagitta_new = sagitta + 0.000013 * norm01();

      //check if smearing changed charge, too
      if( sagitta * sagitta_new < 0.0 ) charge = -1.0 * charge;
      sagitta = sagitta_new;

      //update back radius
      radius = sagitta/2.0 + pow(m_param_indet_eta_coverage, 2) / (8.0 * sagitta);
      
      //update back pT
      trk_pt = radius * ( 0.3 *  m_param_B_field );

      //change px, py but keep pz
      double trk_px = trk_pt * cos( particle.momentum().phi() );
      double trk_py = trk_pt * sin( particle.momentum().phi() );
      double trk_pz = trk_pt / tantheta;

      FourMomentum smeared_p( charge, trk_px, trk_py, trk_pz );
      m_indet_tracks.push_back( smeared_p );
    }

    log << Log::DEBUG << "InDet: found " << m_indet_tracks.size() << " charged tracks" << endl;
    log << Log::DEBUG << "Missed tracks / All Tracks = " 
	<< missed_tracks << " / " << all_tracks << endl;
    return status;
  }


  //////////////////////////////////////


  bool DetectorFastSim::fillCalo(const Event& event)
  {
    Log log = getLog();
    log << Log::DEBUG << "Filling Calorimeter Cells" << endl; 

    bool status = true;

    m_fs_etmiss_pt  = 0.0;
    m_fs_etmiss_phi = 0.0;
    m_fs_sumet      = 0.0;


    const VetoedFinalState& fs = applyProjection<VetoedFinalState>(event, "VFSVISIBLE");
    if (fs.isEmpty()) {
      return false; // skip if empty!
    }
     
    ParticleVector particles = fs.entities();
    //log << Log::DEBUG << "Looping over " << particles.size() << endl;
    foreach( Particle& particle, particles ) {
      const double pT   = particle.momentum().pT();
      if( pT < 0.5 ) continue; //particle doesn't reach calo surface

      const double eta = particle.momentum().eta();
      if( fabs(eta) > m_param_calo_eta_coverage ) continue;

      const double phi = particle.momentum().phi(); 
      const int apid = abs( particle.pdgId() );

      if( apid == 11 && inCrack(eta) ) continue;
      if( apid == 13 && pT < 1.0 )   continue;

      double E = particle.momentum().t();

      double e_em  = 0.0;
      double e_had = 0.0;

      //Is it EM?  e, gamma, pi0, eta
      if( apid == 11 || apid == 22 || apid == 111 || apid == 221 ) {
        e_em  = ECAL_FRAC * E;
	e_had = ( 1 - ECAL_FRAC ) * E;
	
	//calculate resolution
	e_em += e_em * m_param_calo_em_res_const * norm01() + 
	  sqrt(e_em) * m_param_calo_em_res_sqrtE * norm01();
	if(e_em < 0.0) e_em = 0.0;

	e_had += sqrt(e_had) * m_param_calo_had_res_sqrtE * norm01();
	if(e_had < 0.0) e_had = 0.0;
      } // end EM particle
      else if( apid == 13 ) {
	double mu_e_em = 0.5 + norm01() * 0.1;
	if( mu_e_em < 0.0 ) mu_e_em = 0.0;
	
	e_em = 0.0;
	if( E > mu_e_em ) {
	  e_em = mu_e_em;
	  E -= mu_e_em;
	}
	else {
	  e_em = E;
	  E = 0.0;
	}

	double mu_e_had = 2.0 + norm01() * 0.4;
	if( mu_e_had < 0.0 ) mu_e_had = 0.0;

	e_had = 0.0;
	if( E > mu_e_had ) {
	  e_had = mu_e_had;
	  E -= mu_e_had;
	}
	else {
	  e_had = E;
	  E = 0.0;
	}
	
      } // end muon
      else {
	//must be hadronic
	double ecal_frac = 0.25 + 0.05 * norm01();
	if( ecal_frac < 0.0 ) ecal_frac = 0.0;
	else if( ecal_frac > 1.0 ) ecal_frac = 1.0;

	//apply resolution
	if( rand01() < 0.5 ) {
	  e_em = 0.5 + 0.1 * norm01();
	  if( e_em > E ) ecal_frac = 1.0;
	  else if( e_em < 0.0 )  ecal_frac = 0.0;
	  else ecal_frac = e_em / E;
	}
	
	e_em  = ecal_frac * E;
	e_had = (1.0 -  ecal_frac ) * E;
	
	//apply resolution
	//~80%/sqrt(E) for hadrons in EM calorimeter
	e_em += sqrt(e_em) * m_param_calo_em_res_sqrtE * norm01();
	if( e_em < 0.0 ) e_em = 0.0;

	e_had +=  sqrt(e_had) * m_param_calo_had_res_sqrtE * norm01();
	if( e_had < 0.0 ) e_had = 0.0;

      } // end hadronic particle
				

      // now fill arrays
      const int eta_index = get_eta_index( eta );
      const int phi_index = get_phi_index( phi );

      //log << Log::DEBUG << "em=" << e_em << " e_had=" << e_had << endl; 
      //log << Log::DEBUG << "eta=" << eta << " (" << eta_index << ")  " 
      //                << "phi=" << phi << " (" << phi_index << ")" << endl;

      //m_calo_em[eta_index][phi_index]  += e_em;
      //m_calo_had[eta_index][phi_index] += e_had;
      double em_value  = e_em  + m_calo_em->GetBinContent( eta_index, phi_index );
      double had_value = e_had + m_calo_had->GetBinContent( eta_index, phi_index );

      m_calo_em->SetBinContent( eta_index, phi_index, em_value );
      m_calo_had->SetBinContent( eta_index, phi_index, had_value );

      // sum transverse energy
      const double theta = particle.momentum().theta();
      
      m_fs_etmiss_px += (e_em + e_had) * sin(theta)*cos(phi);
      m_fs_etmiss_py += (e_em + e_had) * sin(theta)*sin(phi);
      m_fs_sumet     += (e_em + e_had) * sin(theta);
      m_fs_calo_met   = sqrt( pow(m_fs_etmiss_px,2) + pow(m_fs_etmiss_py,2) );

    } // end loop over particles

 
    return status;
  }


  //////////////////////////////////////


  bool DetectorFastSim::findClusters()
  {
    Log log = getLog();
    log << Log::DEBUG << "Finding clusters" << endl; 

    vector<fastjet::PseudoJet> pseudojets;
    int counter = 1;
    for( unsigned int eta_bin = 0 ; eta_bin < 2*m_param_calo_cells_eta ; ++eta_bin ) {
      const double eta  = get_true_eta( eta_bin );
      const double theta = 2.0 * atan( exp( -eta ) );

      for( unsigned int phi_bin = 0 ; phi_bin < m_param_calo_cells_phi ; ++phi_bin ) {
	const double phi = phi_bin * m_calo_cells_phi_width;

	const double E_tower =  m_calo_em->GetBinContent(eta_bin, phi_bin) 
	  + m_calo_had->GetBinContent(eta_bin, phi_bin);
			
	const double E_T_tower = E_tower * sin(theta);

	//if( E_T_tower > 0.0 ) {
	//  log << Log::DEBUG << "Cell (" << eta << ", " << phi << ") E_T=" <<  E_T_tower << endl;
	//}

	if( E_T_tower < 1.0 ) continue;
	//log << Log::DEBUG << "Cell (" << eta << ", " << phi << ") E_T=" <<  E_T_tower << endl;

	const double px = E_T_tower * cos(phi);
	const double py = E_T_tower * sin(phi);
	const double pz = E_tower * cos(theta);

	//log << Log::DEBUG << "Tower 4mom=(" << px << ", " << py << ", " << pz << "; " << E_tower << ")" << endl;
	FourMomentum tower(E_tower, px, py, pz);
	m_calo_tracks.push_back( tower );
	//m_calo_tracks[counter] = tower;
	
	//log << Log::DEBUG << "Calo trk: true eta=" << eta << " calc eta=" << tower.eta() << endl;
	//log << Log::DEBUG << "Calo trk: (eta,phi)=(" << eta << ", " << phi << ")=(" 
	//    << eta_bin << ", " << phi_bin << ") ET=" <<  E_T_tower << endl;


	fastjet::PseudoJet pJet(px, py, pz, E_tower);
	pJet.set_user_index(counter);

	pseudojets.push_back( pJet );

	++counter;
      } //loop phi
    } //loop eta

    log << Log::DEBUG << "Running FastJet ClusterSequence construction" << endl;
    m_cluster_seq.reset(new fastjet::ClusterSequence( pseudojets, m_jet_def ) );

    bool found_jets = findJets( m_cluster_seq->inclusive_jets( m_param_jet_ptmin ) );
    if(!found_jets) {
      log << Log::ERROR << "Unable to perform jet finding!" << endl;
    }

    return true;
  }


  //////////////////////////////////////


  bool DetectorFastSim::findJets(const vector<fastjet::PseudoJet>& pseudojets )
  {
    Log log = getLog();

    
    log << Log::DEBUG << "Finding jets" <<  endl; 

    bool status = true;

    

    foreach (const fastjet::PseudoJet& pj, pseudojets) {
      Jet j;
      
      assert( m_cluster_seq.get() );
      const PseudoJets parts = (m_cluster_seq.get())->constituents(pj);

      //log << Log::DEBUG << "PJ has const = " << parts.size() << endl;

      foreach(const fastjet::PseudoJet& p, parts) {
	//int index = p.user_index();
	//j.addParticle( m_calo_tracks[index] );
	FourMomentum track( p.E(), p.px(), p.py(), p.pz() );
	j.addParticle( track );
      }
      
      if( j.momentum().pT() < 20. & fabs( j.momentum().eta() ) > 2.5 ) continue;
      m_raw_jets.push_back( j );
    }

    m_raw_jets_n = m_raw_jets.size();
    //log << Log::DEBUG << "Found " << m_raw_jets_n << " raw jets" << endl;
    foreach(const Jet& j, m_raw_jets) {
      const double jpt  = j.momentum().pT();
      const double jeta = j.momentum().eta();
      const double jphi = j.momentum().phi();

      log << Log::DEBUG << " * Reco: Raw jet: (" <<jeta << ",  " << jphi << ") pT = " << jpt  << endl;

      m_raw_jets_pt->push_back( jpt );      
      m_raw_jets_phi->push_back( jphi );      
      m_raw_jets_eta->push_back( jeta );          
    }

    
    return status;
  }


  /////////////////////////////////////


  bool DetectorFastSim::findElectrons()
  {
    Log log = getLog();
    log << Log::DEBUG << "Finding e/gamma" << endl; 

    bool status = true;

    const double eta_max_barrel = 2.5;

    //loop over calo clusters
    foreach(const FourMomentum& tower, m_calo_tracks) {
      const double eta = tower.eta();
      const double phi = tower.phi();
      const double theta = 2.0 * atan( exp( -eta ) );

      const unsigned int eta_bin = get_eta_index( eta );
      const unsigned int phi_bin = get_phi_index( phi );

      const double E_em      =  m_calo_em->GetBinContent(eta_bin, phi_bin);
      const double E_had     =  m_calo_had->GetBinContent(eta_bin, phi_bin);
      const double E_tower   =  E_em + E_had;
      const double E_T_tower =  E_tower * sin(theta);

      //log << Log::DEBUG << "Calo trk: (eta,phi)=(" << eta << ", " << phi << ")=(" 
      //	    << eta_bin << ", " << phi_bin << ") ET=" <<  E_T_tower << endl;

      //if( E_T_tower != tower.Et() ) 
      //log << Log::DEBUG << "E_T_tower = " << E_T_tower << " true ET = " << tower.Et() << endl;

      if( E_T_tower < 5. || fabs(eta) > eta_max_barrel  ) continue;

      const double E_ratio   = E_em > 0 ? E_had / E_em : 0.0;
      if( E_ratio > 0.5 ) continue; //had > 0.5 * em

      //const double et_isol   = getCaloSumEnergy(eta, phi, 0.1); 
      const double et_isol   = getTowerIsolationSum(eta, phi);
      const double pt_cone01 = getTrackPtSum(eta, phi, 0.01);
      const double pt_cone02 = getTrackPtSum(eta, phi, 0.2);
      const double P         = tower.p().rho();
      const double E_over_p  = P > 0.0 ? E_tower / P : -1.0;

      /*
	log << Log::DEBUG << "Candidate (" << eta << ", " << phi << ") ET = " 
	<< E_T_tower << " et_isol=" << et_isol << " pt_cone02=" << pt_cone02 
	<< " ptcone01=" << pt_cone01 << " E/P=" << E_over_p << endl;
      */
      //gamma candidate
      if( et_isol / E_T_tower < 0.5 
	  && pt_cone02 < 50.0 
	  && pt_cone01 > 0. && pt_cone01 < 3 ) {
	log << Log::DEBUG << " * Reco: gamma: (" << eta << ", " << phi << ") pT=" << pt_cone01 
	    << endl;
      }

      //electron candidate
      if( et_isol / E_T_tower < 0.1
	  && pt_cone01 > 0.
	  && ( pt_cone02 - pt_cone01 ) < 50.
	  && E_over_p > 0.5 && E_over_p < 3.) {
	
	
	m_fs_ele.push_back( tower );
	m_fs_ele_pt->push_back( tower.pT() );
	m_fs_ele_eta->push_back( eta );
	m_fs_ele_phi->push_back( phi );
	m_fs_ele_iso->push_back( pt_cone02 );
	++m_fs_ele_n;

	//matched?
	bool matched = false;
	foreach( const Particle& true_ele, m_truth_ele ) {
	  const double dR = deltaR( eta, phi, true_ele.momentum() );
	  if( dR < 0.2 ) {
	    matched = true;
	    break;
	  }
	}
	m_fs_ele_matched->push_back( matched );

	log << Log::DEBUG << " * Reco: ele: (" << eta << ", " << phi << ") pT=" 
	    << pt_cone01 << " matched=" << matched << endl;
      }
	    
    }

  
    return status;
  }


  //////////////////////////////////////


  bool DetectorFastSim::doOverlapRemoval()
  {
    Log log = getLog();
    log << Log::DEBUG << "Applying Overlap removal" << endl; 

    bool status = true;

    foreach(const Jet& j, m_raw_jets) {
      const double jeta = j.momentum().eta();
      const double jphi = j.momentum().phi();
      const double jpt  = j.momentum().pT();

      bool overlap = false;
      //foreach(const Particle& ele, m_truth_ele) {
      foreach(const FourMomentum& ele, m_fs_ele) {
	const double ele_eta = ele.eta(); //ele.momentum().eta();
	const double ele_phi = ele.phi(); //ele.momentum().phi();

	const double dR = sqrt( pow( jeta - ele_eta, 2) + pow( jphi - ele_phi, 2 ) );
	if( dR < 0.3 ) {
	  log << Log::DEBUG << "e/j overlap" << endl;
	  overlap = true;
	  break;
	}
      }

      if(!overlap) {
	m_fs_jets.push_back( j );
	m_fs_jets_eta->push_back( jeta );
	m_fs_jets_phi->push_back( jphi );
	m_fs_jets_pt->push_back( jpt );
	++m_fs_jets_n;

	bool matched = false;
	foreach( const Jet& true_jet, m_truth_jets ) {
	  const double dR = deltaR( jeta, jphi, true_jet.momentum() );
	  if( dR < 0.2 ) {
	    matched = true;
	    break;
	  }
	}
	m_fs_jets_matched->push_back( matched );

	log << Log::DEBUG << " * Reco: jet: (" << jeta << ", " << jphi << ") pT=" 
	    << jpt << " matched=" << matched <<endl; 

      }
    }
    log << Log::DEBUG << " * Reco: found " << m_fs_jets_n << " jets" << endl;

    return status;

  }


  //////////////////////////////////////


  bool  DetectorFastSim::doTagging()
  {
    Log log = getLog();
    log << Log::DEBUG << "Tagging jets" << endl; 
    
    bool status = true;
    static const double D_ctau = 439; //for D mesons, um

    foreach( const Jet& jet, m_fs_jets) {
      int tag = 0;
      double tag_weight = 0;

      // truth-based tag (reliable)
      foreach( const Particle& quark, m_truth_bquarks ) {
	const double dR = deltaR(jet.momentum(), quark.momentum());

	if( dR < 0.4 ) {
	  tag = 5;
	  ++m_fs_bjets_n;

	  break;
	}
      } // loop over bquarks

      if( tag == 0 ) {
	// no b quark? try taus
	foreach( const Particle& tau, m_truth_taus ) {
	  const double dR = deltaR(jet.momentum(), tau.momentum());
	  
	  if( dR < 0.4 ) {
	    tag = 15;
	    ++m_fs_tjets_n;
	    
	    break;
	  }
	} // loop over taus, if any
      }
      m_fs_jets_tag->push_back( tag );

      //vtx-based tag (test)
      /*
      foreach( const FourMomentum& track, jet.momenta() ){
	const double pt     = track.pT();
	const double theta  = track.theta();
	const double P      = track.mod();

	 

	//values taken from atlas detector
	const double angle = pt / P;
	const double d0 = D_ctau * sin(angle);
	const double z0 = rand01() * D_ctau;

	const double sigma_z0 = 11. * norm01() + (60. / pt * sqrt( sin(theta))) * norm01();
	const double sigma_d0 = 70. * norm01() + (100./ pt * sqrt( pow( sin(theta), 3 ) ) ) * norm01();

	const double Sd0 = d0 / sigma_d0;
	const double Sz0 = z0 / sigma_z0;

	tag_weight += Sz0 + Sd0;
      }
      */
      m_fs_jets_tagweight->push_back( tag_weight );
     
    }

    log << Log::DEBUG << " * Reco: btag: N=" << m_fs_bjets_n << endl;
    log << Log::DEBUG << " * Reco: taujet: N=" << m_fs_tjets_n << endl;


    return status;
  }

  
  //////////////////////////////////////


  bool DetectorFastSim::findMuons(const Event& event)
  {
    Log log = getLog();
    log << Log::DEBUG << "Finding muons" << endl; 

    bool status = true;

    foreach(const Particle& lept, m_truth_mu ) {
      
      //const double pt  = lept.momentum().pT();
      //if( pt < 6.*GeV ) continue;
     
      const double eta = lept.momentum().eta();
      const double phi = lept.momentum().phi();

     
      //is there a corresponding track in the inner detector?
      bool matched = false;
      foreach(const FourMomentum& track, m_indet_tracks) {
	const double dR = deltaR( lept.momentum(), track );
	if( dR < 0.3 ) {
	  if( lept.momentum().pT() - track.pT() > 5. ) continue;

	  matched  = true;

	  //appy detector geo efficiency. 
	  if(fabs(eta) > 2.7 ) continue;
	  if( fabs(eta) < 0.1) {
	    if( fabs(eta) < 0.05 ) {
	      //hole in central region
	      continue;
	    }
	    else {
	      //give it a chance
	      if( rand01() > 0.5 ) continue; 
	    }
	  }

	  const double etcone30 = getTrackPtSum( eta, phi, 0.3 );
	  //const double eta_bin  = get_eta_index( track.eta() );
	  //const double phi_bin  = get_phi_index( track.phi() );

	  const double E        = getCaloSumEnergy( eta, phi, 0.3 );
	
	  FourMomentum muon( E, track.px(), track.py(), track.pz() );
	  m_fs_mu.push_back( muon );
	  m_fs_mu_pt->push_back( muon.pT() );
	  m_fs_mu_eta->push_back( muon.eta() );
	  m_fs_mu_phi->push_back( muon.phi() );
	  m_fs_mu_iso->push_back( etcone30 );
	  ++m_fs_mu_n;

	  //matched?
	  bool matched = false;
	  foreach( const Particle& true_mu, m_truth_mu ) {
	    const double dR = deltaR( eta, phi, true_mu.momentum() );
	    if( dR < 0.2 ) {
	      matched = true;
	      break;
	    }
	  }
	  m_fs_mu_matched->push_back( matched );

	  log << Log::DEBUG << " * Reco: muon: (" << eta << ", " << phi 
	      << ") pT=" << muon.pT() << " matched=" << matched << endl;

	  break;
	}//matched track
      }//loop over indet tracks
    }//loop over true muons

    return status;
  }


  //////////////////////////////////////


  
  bool DetectorFastSim::calculateMissingET()
  {
    Log log = getLog();
    log << Log::DEBUG << "Calculating MissingET" << endl; 

    bool status = true;

    foreach( const FourMomentum& muon, m_fs_mu ) {
      m_fs_etmiss_px += muon.px();
      m_fs_etmiss_py += muon.py();
      m_fs_sumet    += muon.Et();
    }

    m_fs_etmiss_px  = -1 * m_fs_etmiss_px;
    m_fs_etmiss_py  = -1 * m_fs_etmiss_py;
    m_fs_etmiss_pt  = sqrt( pow(m_fs_etmiss_px,2) + pow(m_fs_etmiss_py,2) );
    m_fs_etmiss_phi = atan2( m_fs_etmiss_px, m_fs_etmiss_py );

    log << Log::DEBUG << "Calo MET=" << m_fs_calo_met << endl;
    log << Log::DEBUG << "Corr MET=" << m_fs_etmiss_pt << endl;
    log << Log::DEBUG << "SumET=" << m_fs_sumet << endl;

    return status;
  }


  //////////////////////////////////////



  bool DetectorFastSim::applyTrigger()
  {
    Log log = getLog();
    log << Log::DEBUG << "Applying Trigger" << endl; 

    bool status = true;

    status = triggerChain_Leptons();
    if(!status) {
      log << Log::INFO << "Impossible apply trigger on leptons" << endl;
      return status;
    }

    status = triggerChain_ETmiss();
    if(!status) {
      log << Log::INFO << "Impossible apply trigger on ETmiss" << endl;
      return status;
    }

    status = triggerChain_Jets();
    if(!status) {
      log << Log::INFO << "Impossible apply trigger on jets" << endl;
      return status;
    }


    for( std::map< std::string, bool >::const_iterator trigger = m_trigger_chains.begin() ;
	 trigger != m_trigger_chains.end() ; ++trigger ) {
      if( trigger->second == 1 ) 
	log << Log::DEBUG << " * Trigger: " << trigger->first << endl; 
    }
   
    
    return status;
  }


  //////////////////////////////////////

#define TRIGGER_TURNON_FUNCTION(x)  tanh(0.07*( (x) - 9.0 )) 

  bool DetectorFastSim::triggerChain_Leptons()
  {
    Log log = getLog();

    bool status = true;

    // consider 3 threshold: 10, 20, 30
    m_trigger_chains["e10"] = false;
    m_trigger_chains["e20"] = false;
    m_trigger_chains["e30"] = false;

    m_trigger_chains["mu10"] = false;
    m_trigger_chains["mu20"] = false;
    m_trigger_chains["mu30"] = false;

    // electrons

    foreach( const FourMomentum& ele, m_fs_ele ) {
      if( inCrack( ele.eta() ) ) continue; // missed completely

      const double pt = ele.pT() * ( 1. +  norm01() *  m_param_calo_etmiss_res );

      if( pt < 10. ) continue;

      if( rand01() > TRIGGER_TURNON_FUNCTION( pt ) ) {
	log << Log::DEBUG << " * Trigger: missed ele due to resolution" << endl;
	continue;
      }
      //const double r = rand01();
      

      if( pt > 10. ) {
	m_trigger_chains["e10"] = true;
      
	if( pt > 20. ) {
	  m_trigger_chains["e20"] = true;
      
	  if( pt > 30. ) {
	    m_trigger_chains["e30"]  = true;
	  }
	}
      }
    }

    // muons
    foreach( const FourMomentum& muon, m_fs_mu) {
      const double eta = muon.eta();

      //appy detector geo efficiency. 
      if(fabs(eta) > 2.7 ) continue;
      if( fabs(eta) < 0.1) {
	if( fabs(eta) < 0.05 ) {
	  //hole in central region
	  continue;
	}
	else {
	  //give it a chance
	  if( rand01() > 0.5 ) continue; 
	}
      }

      const double pt = muon.pT() * ( 1. +  norm01() *  m_param_calo_etmiss_res );

      if( pt < 10. ) continue;

      if( rand01() > TRIGGER_TURNON_FUNCTION( pt ) ) {
	log << Log::DEBUG << " * Trigger: missed mu due to resolution" << endl;    
	continue;
      }

      if( pt > 10. ) {
	m_trigger_chains["mu10"] = true;
      
	if( pt > 20. ) {
	  m_trigger_chains["mu20"] = true;
      
	  if( pt > 30. ) {
	    m_trigger_chains["mu30"] = true;
	  }
	}
      }
    }

    m_trig_e10 = m_trigger_chains["e10"]; 
    m_trig_e20 = m_trigger_chains["e20"]; 
    m_trig_e30 = m_trigger_chains["e30"]; 

    m_trig_mu10 = m_trigger_chains["mu10"]; 
    m_trig_mu20 = m_trigger_chains["mu20"]; 
    m_trig_mu30 = m_trigger_chains["mu30"]; 


    return status;
  }

  
  //////////////////////////////////////


  bool DetectorFastSim::triggerChain_ETmiss()
  {
    bool status = true;

    m_trigger_chains["xe20"] = false;
    m_trigger_chains["xe30"] = false;
    m_trigger_chains["xe50"] = false;
    m_trigger_chains["xe100"] = false;

    const double met_smeared = m_fs_etmiss_pt  +  m_fs_etmiss_pt * norm01() *  m_param_calo_etmiss_res;
  
    if( met_smeared >=  20. ) m_trigger_chains["xe20"]   = true;
    if( met_smeared >=  30. ) m_trigger_chains["xe30"]   = true;
    if( met_smeared >=  50. ) m_trigger_chains["xe50"]   = true;
    if( met_smeared >= 100. ) m_trigger_chains["xe100"]  = true;


    m_trig_xe20  = m_trigger_chains["xe20"];
    m_trig_xe30  = m_trigger_chains["xe30"];
    m_trig_xe50  = m_trigger_chains["xe50"];
    m_trig_xe100 = m_trigger_chains["xe100"];

    
    
    return status;

  } 


  //////////////////////////////////////

    
  bool DetectorFastSim::triggerChain_Jets() 
  {
    bool status = true;

    //static unsigned int prescale_counter;

    m_trigger_chains["4j20"] = false;
    m_trigger_chains["2j40"] = false;
    m_trigger_chains["3j40"] = false;
    m_trigger_chains["2j100"] = false;

    //if( prescale_counter < m_param_trigger_jet_prescale ) return true;

    unsigned int n_jets_20  = 0;
    unsigned int n_jets_40  = 0;
    unsigned int n_jets_100 = 0;

    foreach( const Jet& jet, m_fs_jets ) {
      const double pt = jet.momentum().pT();

      const double pt_smeared = pt + ( pt * norm01() * m_param_calo_etmiss_res );

      if( pt_smeared >= 20. ) {
	++n_jets_20;

	if( pt_smeared >= 40. ) {
	  ++n_jets_40;
	  
	  if( pt_smeared >= 100. ) {
	    ++n_jets_100;
	  }
	}
      } // pt thresholds
      
    }

    if( n_jets_20 >= 4 )  m_trigger_chains["4j20"] = true;

    if( n_jets_40 >= 2 ) {
      m_trigger_chains["2j40"] = true;

       if( n_jets_40 >= 3 ) {
	 m_trigger_chains["3j40"] = true;
       }
    }

    if( n_jets_100 >= 2 ) {
       m_trigger_chains["2j100"] = true;
    }


    m_trig_4j20  = m_trigger_chains["4j20"];
    m_trig_2j40  = m_trigger_chains["2j40"];
    m_trig_3j40  = m_trigger_chains["3j40"];
    m_trig_2j100 = m_trigger_chains["2j100"];


    //++prescale_counter;
    return status;

  }


}; //end namespace Rivet
