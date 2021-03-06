
// $Id: $
// Include files

// from Gaudi
#include "GaudiKernel/AlgFactory.h"

#include "GaudiKernel/ToolFactory.h" 
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "GaudiKernel/KeyedObject.h"

#include "LoKi/Particles36.h"
#include "LoKi/Constants.h"
#include "LoKi/BasicFunctors.h"

// local
#include "DoubleJpsiSel.h"

//#include <Kernel/DVAlgorithm.h>
#include <Kernel/GetIDVAlgorithm.h>
#include "LHCbMath/ParticleParams.h"

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp> 
#include <boost/regex.hpp> 

#include "Kernel/ParticleFilters.h"

#include "Kernel/ILHCbMagnetSvc.h"
//#include "Kernel/IOnOffline.h"

#include "Event/ODIN.h"
#include "Event/L0DUReport.h"

#include "Event/MCHeader.h"
#include "Event/GenCollision.h"
#include "Event/GenHeader.h"
#include "Event/MCVertex.h"
#include "HepMC/GenEvent.h"
#include "Event/Node.h"
#include "Event/FitNode.h"
#include "Kernel/IParticle2MCWeightedAssociator.h"
#include "Kernel/IParticle2MCAssociator.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

// ===========================================================================
// for Kullback divergence
// ===========================================================================
#include "GaudiAlg/TupleObj.h"
#include "Event/Particle.h"

// ===========================================================================
// for MCreconstructible
// ===========================================================================

#include "MCInterfaces/IMCReconstructible.h"
#include "MCInterfaces/IMCReconstructed.h"

/*
 * DecayTreeFitter
 */
#include "DecayTreeFitter/Fitter.h"



// ===========================================================================
// GlobalEventCounters
// ===========================================================================
#include "Kernel/Counters.h"
// ===========================================================================
// Primary vertex
// ===========================================================================
#include "Event/VertexBase.h"
#include "Event/RecVertex.h"
// ===========================================================================
// Filter particles
// ===========================================================================
#include "Kernel/ParticleFilters.h"
// ===========================================================================
// Magnetic field
// ===========================================================================
#include "Kernel/ILHCbMagnetSvc.h"
//#include "Kernel/IOnOffline.h"
// ===========================================================================
// ODIN and Trigger info
// ===========================================================================
#include "Event/ODIN.h"
#include "Event/L0DUReport.h"
//#include "Event/HltDecReports.h"
#include "Kernel/IANNSvc.h"
// ===========================================================================
// Reconstruction info
// ===========================================================================
#include "Event/Track.h"
#include "Event/RecSummary.h"


// ===========================================================================
// from Gaudi
// ===========================================================================
#include "GaudiKernel/AlgFactory.h"
 
#include "GaudiAlg/TupleObj.h"
#include "Event/Particle.h"
 


// ===========================================================================
// DecayTreeFitter
// ===========================================================================
#include "DecayTreeFitter/Fitter.h"
// ===========================================================================
// Lifetime fitter
// ===========================================================================
#include "Kernel/ILifetimeFitter.h"
// ===========================================================================
// HepMC and MonteCarlo stuff
// ===========================================================================
#include "HepMC/GenEvent.h"
#include "Event/HepMCEvent.h"
#include "Event/GenHeader.h"
#include "Event/MCHeader.h"
#include "Event/MCParticle.h"
#include "Event/MCVertex.h"
// ===========================================================================
// BackgroundCategory and MC associators
// ===========================================================================
#include "Kernel/IParticle2MCAssociator.h"
// ===========================================================================
// LoKi stuff
// ===========================================================================
//#include "LoKi/AParticles.h"
//#include "LoKi/Particles26.h"
//#include "LoKi/CoreCuts.h"
//#include "LoKi/ParticleCuts.h"
//#include "LoKi/PhysTypes.h"
///////#include "LoKi/AParticleCuts.h"
/////

////////#include "LoKi/ParticleCuts.h"
//#include "LoKi/Particles.h"
//#include "LoKi/ATypes.h"
// ===========================================================================
// Standard Libraries
// ===========================================================================
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
// ===========================================================================
// Boost Libraries
// ===========================================================================
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/regex.hpp>
 
using namespace boost::lambda;



using namespace boost::lambda;

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( DoubleJpsiSel );


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
DoubleJpsiSel::DoubleJpsiSel(const std::string &name,
			     ISvcLocator*       pSvcLocator)
  : DaVinciTupleAlgorithm (name, pSvcLocator)
  , m_inputType("DST")
  , m_jPsiID(0)
  , m_jPsiMass(0.)
  , m_recible(0)
  , m_rected(0)
  , m_dva(NULL)
  , m_dist(NULL)
  , m_TriggerTisTosName("TriggerTisTos")
  , m_L0TriggerTisTosName("L0TriggerTisTos")
  , m_stateProvider("TrackStateProvider")
  , m_HLT1TISList(0)
  , m_HLT1TOSList(0)
  , m_HLT2TISList(0)
  , m_HLT2TOSList(0)
  , m_L0TISList(0)
  , m_L0TOSList(0)
  , m_TriggerTisTosTool(NULL)
  , m_L0TriggerTisTosTool(NULL)
  , m_associator(NULL)
  , m_associatorType("MCMatchObjP2MCRelator")
{

  declareProperty("MaxJpsiChi2", m_maxJpsiChi2 = 20.0);
  declareProperty("MaxJpsiMass", m_maxJpsiMass = 3.2 * Gaudi::Units::GeV);
  declareProperty("MinJpsiMass", m_minJpsiMass = 3.0 * Gaudi::Units::GeV);
  declareProperty("MaxDTFChi2",        m_maxDTFChi2Ndof    = 5.0);
  declareProperty("MaxMuChi2",   m_maxMuChi2NDof = 5.0);
  declareProperty("MinMuDll",    m_minMuDll      = 0.0);
  declareProperty("MinMuPt",     m_minMuPt       = 0.65 * Gaudi::Units::GeV);
  declareProperty("MinMomentum", m_minMomentum    = 5.0 * Gaudi::Units::GeV);
  declareProperty("MassWindow",  m_jPsiMassWin = 10.0 * Gaudi::Units::GeV);

  declareProperty("HLT1TISList", m_HLT1TISList = std::vector<std::string>(0));
  declareProperty("HLT1TOSList", m_HLT1TOSList = std::vector<std::string>(0));
  declareProperty("HLT2TISList", m_HLT2TISList = std::vector<std::string>(0));
  declareProperty("HLT2TOSList", m_HLT2TOSList = std::vector<std::string>(0));
  declareProperty("L0TISList",   m_L0TISList   = std::vector<std::string>(0));
  declareProperty("L0TOSList",   m_L0TOSList   = std::vector<std::string>(0));

  declareProperty("TriggerTisTosName",   m_TriggerTisTosName);
  declareProperty("L0TriggerTisTosName", m_L0TriggerTisTosName);

  declareProperty("FillVerbose",            m_verbose                = true);
  declareProperty("FillAllJpsi",            m_fillAllJpsi            = true);
  declareProperty("FillSelectedJpsi",       m_fillSelectedJpsi       = true);
  declareProperty("FillAllDoubleJpsi",      m_fillAllDoubleJpsi      = true);
  declareProperty("FillSelectedDoubleJpsi", m_fillSelectedDoubleJpsi = true);
  declareProperty("MCTen", m_MCTen = false);
  declareProperty("MCTruth",m_MCTruth = false);
  declareProperty("fillMCTruthTuple",m_fillMCTruthTuple = false);
  declareProperty("NoAccCutMC",m_NoAccCutMC = false);
  declareProperty("InputType",              m_inputType);
  declareProperty( "ChargedThetaMin" , m_chargedThetaMin = 10 * Gaudi::Units::mrad ) ;
  declareProperty( "ChargedThetaMax" , m_chargedThetaMax = 400 * Gaudi::Units::mrad ) ;
  declareProperty( "NeutralThetaMin" , m_neutralThetaMin = 5 * Gaudi::Units::mrad ) ;
  declareProperty( "NeutralThetaMax" , m_neutralThetaMax = 400 * Gaudi::Units::mrad ) ;

  declareProperty("SelectFillDLL",          m_selFillDll             = 0x1);

  m_backCatTypes.push_back("BackgroundCategoryViaRelations");
  m_backCatTypes.push_back("BackgroundCategory");
 
  declareProperty("IBackgroundCategoryTypes", m_backCatTypes);
  declareProperty("IP2MCPAssociatorType",     m_associatorType);
}

//=============================================================================
// Destructor
//=============================================================================

DoubleJpsiSel::~DoubleJpsiSel() {} 


/******************************************************************************
 *                                                                            *
 * DoubleJpsiSel::initialize                                                  *
 *                                                                            *
 ******************************************************************************/
StatusCode DoubleJpsiSel::initialize() {
  if (msgLevel(MSG::DEBUG)) debug() << "==> Initialize" << endmsg;

  const LHCb::ParticleProperty* property = ppSvc()->find("p+");
  if (!property){
    err() << "Cannot find particle property for proton" << endmsg ;
    return StatusCode::FAILURE;
  }
  
  property = ppSvc()->find("mu+");
  if (!property){
    err() << "Cannot find particle property for muon" << endmsg ;
    return StatusCode::FAILURE;
  }
  m_pidMuon = property->particleID().pid();

  StatusCode sc = DaVinciTupleAlgorithm::initialize(); 
  if (sc.isFailure()) return sc;
  
  //
  // Initialize the magnetic field service.
  //

  m_magnetSvc = svc<ILHCbMagnetSvc>("MagneticFieldSvc", true);

  //
  // Initialize the distance calculator
  //
 
  m_dva = Gaudi::Utils::getIDVAlgorithm (contextSvc()) ;
  if (0 == m_dva) return Error("Couldn't get parent IDVAlgorithm", 
			       StatusCode::FAILURE);  
  m_dist = m_dva->distanceCalculator();

  //
  // Initialize the TISTOS tools, one for L0 and one for HLT(1,2)
  //
  m_tck = -1;


  m_associator = tool<IParticle2MCAssociator>(m_associatorType, this);
  if (!m_associator){
    err() << "Cannot find Particle2MC associator" << endmsg ;
      return StatusCode::FAILURE;
  }


  // m_TriggerTisTosTool   = tool<ITriggerTisTos>(m_TriggerTisTosName,
	//				       "TriggerTisTos", this);
  // m_L0TriggerTisTosTool = tool<ITriggerTisTos>(m_L0TriggerTisTosName,
  //			       "L0TriggerTisTos", this);

  //
  // Load the J/psi PDG id and its mass central value.
  //


  //
  // Lifetime fitter tool
  //
 
  m_fitLifetime = tool<ILifetimeFitter>("PropertimeFitter", this);
  if (!m_fitLifetime) return Error("Couldn't get ILifetimeFitter",
                                   StatusCode::FAILURE);


  //
  // Track state provider for the DecayTreeFitter
  //
 
   if (!m_stateProvider.empty()){
     sc = m_stateProvider.retrieve();
     if (sc.isFailure()) return Error("Couldn't retrieve TrackStateProvider",
                                     StatusCode::FAILURE);
   }



  info() << "MC Particle location  "
         << LHCb::MCParticleLocation::Default
         << endmsg;
 
  std::vector<std::string>::const_iterator ibkg;
  for (ibkg = m_backCatTypes.begin(); ibkg != m_backCatTypes.end(); ++ibkg){
    m_bkgs.push_back(tool<IBackgroundCategory>(*ibkg, *ibkg, this));
  }

  
  m_associator = tool<IParticle2MCAssociator>(m_associatorType, this);
  if (!m_associator){
    err() << "Cannot find Particle2MC associator" << endmsg ;
      return StatusCode::FAILURE;
  }
  
  

  m_reconstructible = tool<IMCReconstructible>("MCReconstructible");
  m_reconstructed   = tool<IMCReconstructed>("MCReconstructed");

 
  const LHCb::ParticleProperty* JpsiProperty = ppSvc()->find("J/psi(1S)");
  if (!JpsiProperty){ 
    err() << "Cannot find particle property for J/psi(1S)" << endmsg ;
    return StatusCode::FAILURE;
  }

  m_jPsiID   = LHCb::ParticleID(JpsiProperty->pdgID());
  m_jPsiMass = JpsiProperty->mass();
  m_pidJpsi = 443;
  // m_pidJpsi  = property->particleID().pid();

  info() << "J/psi mass = " << m_jPsiMass
	 << " pid = "       << m_jPsiID << endmsg;

  m_recible = tool<IMCReconstructible>("MCReconstructible");
  m_reconstructible = tool<IMCReconstructible>("MCReconstructible");
  m_reconstructed   = tool<IMCReconstructed>("MCReconstructed");


  //  m_Hep2MC          = tool<IHepMC2MC>("LoKi::HepMC2MC");

  //
  // Psi(2S), psi(3770) and chi_b2(1P) properties
  //
  
  property = ppSvc()->find("psi(2S)");
  if (!property){
    err() << "Cannot find particle property for psi(2S)" << endmsg ;
    return StatusCode::FAILURE;
  }
  //m_massPsi2S = property->mass();
    m_pidPsi2S  = property->particleID().pid();
 
   property = ppSvc()->find("chi_b2(1P)");
    if (!property){
      err() << "Cannot find particle property for chi_b2(1P)" << endmsg ;
      return StatusCode::FAILURE;
    }
    // m_massChiB21P = property->mass();
    m_pidChiB21P  = property->particleID().pid();


  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode DoubleJpsiSel::execute() {
  
  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;
  StatusCode sc = StatusCode::SUCCESS ;

  const LHCb::Particle::Range particles = this->particles();

  if(m_MCTruth)sc = MCTruth(); //single Jpsi
  //  if(m_MCTruth)sc = MCTruthTest();
  // if(m_MCTruth)sc = MCTruthDoubleJpsi(); // for double jpsi 
  else {
    if (m_MCTen)sc = MCTenLoop(particles);
      else sc = MakeJpsi(particles); //for normal reco
    //else sc = MakeDoubleJpsi(particles); //for double jpsi MC 
    // else sc = FromBDTFTest(particles);
  }
      
  setFilterPassed(true); 

  return StatusCode::SUCCESS;
}




/* ***********************************************************************(*** *
 *                                                                             *
 * DoubleJpsiCheckMuons:                                                       *
 *                                                                             *
 * Make sure that the muons in two J/psi are not the same particle.            *
 *                                                                             *
 * *************************************************************************** */
StatusCode
DoubleJpsiSel::DoubleJpsiCheckMuons(const LHCb::Particle *muPlus1,
                                   const LHCb::Particle *muMinus1,
                                   const LHCb::Particle *muPlus2,
                                   const LHCb::Particle *muMinus2){
 
  counter("Muon 1")++;
  
  //
  // Are the muons the same?
  //
 
  if (muPlus1->key()  == muPlus2->key())  return StatusCode::FAILURE;
  if (muMinus1->key() == muMinus2->key()) return StatusCode::FAILURE;
  counter("Muon 2")++;

  //
  // Are the muons the same protoparticle?
  //
 
  const LHCb::ProtoParticle *protoPlus1  = muPlus1->proto();
  const LHCb::ProtoParticle *protoPlus2  = muPlus2->proto();
 
  const LHCb::ProtoParticle *protoMinus1 = muMinus1->proto();
  const LHCb::ProtoParticle *protoMinus2 = muMinus2->proto();
 
  if (protoPlus1->key()  == protoPlus2->key())  return StatusCode::FAILURE;
  counter("Muon 3")++;
  if (protoMinus1->key() == protoMinus2->key()) return StatusCode::FAILURE;
  counter("Muon 4")++;

  //
  // Are the muons the same track?
  // 
  const LHCb::Track *trackPlus1  = protoPlus1->track();
  const LHCb::Track *trackPlus2  = protoPlus2->track();
 
  const LHCb::Track *trackMinus1 = protoMinus1->track();
  const LHCb::Track *trackMinus2 = protoMinus2->track();
 
  if (trackPlus1->key()  == trackMinus1->key())  return StatusCode::FAILURE;
  counter("Muon 5")++;
  if (trackPlus1->key()  == trackPlus2->key())   return StatusCode::FAILURE;
  counter("Muon 6")++;
  if (trackPlus1->key()  == trackMinus2->key())  return StatusCode::FAILURE;
  counter("Muon 7")++;
  if (trackMinus1->key() == trackPlus2->key())   return StatusCode::FAILURE;
  counter("Muon 8")++;
  if (trackMinus1->key() == trackMinus2->key())  return StatusCode::FAILURE;
  counter("Muon 9")++;
  if (trackPlus2->key()  == trackMinus2->key())  return StatusCode::FAILURE;
 
  return StatusCode::SUCCESS;
}

//=============================================================================
// Fill Tuples
//=============================================================================

/**************************************************************************** *
 *                                                                            *
 * FillDataPropertime                                                         *
 *                                                                            *
 * Information on the proper time of a particle.                              #
 *                                                                            *
 **************************************************************************** */
void
DoubleJpsiSel::FillDataPropertime(Tuples::Tuple        &ntuple,
                                 const LHCb::Particle *p,
                                 const std::string     prefix){
  double tau     = -100.0;
  double tauErr  = -100.0;
  double tauChi2 = -100.0;
 
  StatusCode sc1 = ProperTime(p, &tau, &tauErr, &tauChi2);
 
  ntuple->column(prefix + "_TAU",     tau);
  ntuple->column(prefix + "_TAUERR",  tau);
  ntuple->column(prefix + "_TAUCHI2", tau);
 
  return;
}

/**************************************************************************** *
 *                                                                            *
 * ProperTime                                                                 *
 *                                                                            *
 **************************************************************************** */
StatusCode
DoubleJpsiSel::ProperTime(const LHCb::Particle *mother,
                         double               *tau,
                         double               *tauErr,
                         double               *tauChi2){
 
  *tau     = -100.0;
  *tauErr  = -100.0;
  *tauChi2 = -100.0;
 
  const LHCb::VertexBase *pv = m_dva->bestVertex(mother);
 
  StatusCode sc = (!pv) ? StatusCode::FAILURE :
     m_fitLifetime->fit(*pv, *mother, *tau, *tauErr, *tauChi2);
 
  return sc;
}  
 


/**************************************************************************** *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * Information of composite particles in simulation.                          *
 *                                                                            *
 **************************************************************************** */
void
DoubleJpsiSel::FillMCComposite(Tuples::Tuple          &ntuple,
                              const LHCb::MCParticle *p,
                              const std::string       prefix){
  if (!p) return;
 
  ntuple->column(prefix + "_P",              double(p->p()));
  ntuple->column(prefix + "_PT",             double(p->pt()));
  ntuple->column(prefix + "_RAPIDITY",       double(p->momentum().Rapidity()));
  ntuple->column(prefix + "_PSEUDORAPIDITY", double(p->pseudoRapidity()));
  ntuple->column(prefix + "_MASS",           double(p->virtualMass()));
 
  return;
}




void 
DoubleJpsiSel::FillMCTrueMuNtuple(Tuples::Tuple          &ntuple,
				  const LHCb::MCParticle *muon,
				  const std::string       prefix){

  ntuple->column(prefix + "_RECONSTRUCTED",   m_reconstructed->reconstructed(muon));
  ntuple->column(prefix + "_RECONSTRUCTIBLE", m_reconstructible->reconstructible(muon));
  ntuple->column(prefix + "_P",               muon->p());
  ntuple->column(prefix + "_PT",              muon->pt());
  ntuple->column(prefix + "_PSEUDORAPIDITY",  muon->pseudoRapidity());
  ntuple->column(prefix + "_RAPIDITY",        muon->momentum().Rapidity());

  return;
}

void
DoubleJpsiSel::FillMCTrueJpsiNtuple(Tuples::Tuple          &ntuple,
				    const LHCb::MCParticle *jpsi,
				    const std::string       prefix){
  ntuple->column(prefix + "_P",              jpsi->p());
  ntuple->column(prefix + "_PT",             jpsi->pt());
  ntuple->column(prefix + "_RAPIDITY",       jpsi->momentum().Rapidity());
  ntuple->column(prefix + "_PSEUDORAPIDITY", jpsi->pseudoRapidity());
  ntuple->column(prefix + "_MASS",           jpsi->virtualMass());

  return;
}



void
DoubleJpsiSel::FillMCTruePolarizationNtuple(Tuples::Tuple          &ntuple,
					    const LHCb::MCParticle *jpsi,
					    const LHCb::MCParticle *muon,
					    const std::string       prefix){
  Gaudi::LorentzVector jpsi4Mom = jpsi->momentum();
  Gaudi::LorentzVector mu4Mom   = muon->momentum();

  ROOT::Math::Boost boost(jpsi4Mom.BoostToCM());

  const Gaudi::XYZVector boostedMu4Mom   = (boost(mu4Mom)).Vect().unit();
  const Gaudi::XYZVector boostedJpsi4Mom = jpsi4Mom.Vect().unit();

  double cosTheta = fabs(boostedMu4Mom.Dot(boostedJpsi4Mom));

  ntuple->column(prefix + "_COSTHETA", cosTheta);

  return;
}




void DoubleJpsiSel::FillEventNtuple(Tuples::Tuple &ntuple){
  LHCb::ODIN *odin(NULL);
  if (exist<LHCb::ODIN>(LHCb::ODINLocation::Default)){
    odin = get<LHCb::ODIN>(LHCb::ODINLocation::Default);
  } else if (exist<LHCb::ODIN>(LHCb::ODINLocation::Default, false)){
    odin = get<LHCb::ODIN>(LHCb::ODINLocation::Default, false);
  } else {
    std::cout << "Cannot load the ODIN data object" << std::endl;
  }

  const std::string prefix = "Event";

  ntuple->column(prefix + "_TCK", odin->triggerConfigurationKey());
  if (m_magnetSvc) ntuple->column(prefix + "_Polarity", m_magnetSvc->isDown() ? -1 : 1);

  if (m_verbose){
    ntuple->column(prefix + "_RunNumber",   odin->runNumber());
    ntuple->column(prefix + "_EventNumber", double(odin->eventNumber()));
    ntuple->column(prefix + "_BCID",        odin->bunchId());
    ntuple->column(prefix + "_BCType",      odin->bunchCrossingType());
  }
  return;
}

void DoubleJpsiSel::FillFourMuonNtuple(Tuples::Tuple        &ntuple,
				       const LHCb::Particle *fourMuon,
				       const std::string     prefix){
  ntuple->column(prefix + "_P",            fourMuon->p());
  ntuple->column(prefix + "_PT",           fourMuon->pt());
  ntuple->column(prefix + "_RAPIDITY",     fourMuon->momentum().Rapidity());
  ntuple->column(prefix + "_VCHI2",        fourMuon->endVertex()->chi2());
  ntuple->column(prefix + "_VNDOF",        fourMuon->endVertex()->nDoF());
  ntuple->column(prefix + "_VCHI2PERNDOF", fourMuon->endVertex()->chi2PerDoF());

  double minIp     = -1.0;
  double minIpChi2 = -1.0;
  ImpactParameter(fourMuon, &minIp, &minIpChi2);

  ntuple->column(prefix + "_MINIP",     minIp);
  ntuple->column(prefix + "_MINIPCHI2", minIpChi2);

  return;
}


/**************************************************************************** *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * Information of basic particles in simulation.                              *
 *                                                                            *
 **************************************************************************** */
void
DoubleJpsiSel::FillMCParticle(Tuples::Tuple          &ntuple,
                             const LHCb::MCParticle *p,
                             const std::string       prefix){
  if (!p) return;
 
  ntuple->column(prefix + "_P",               double(p->p()));
  ntuple->column(prefix + "_PT",              double(p->pt()));
  ntuple->column(prefix + "_PSEUDORAPIDITY",  double(p->pseudoRapidity()));
  ntuple->column(prefix + "_RAPIDITY",        double(p->momentum().Rapidity()));
  ntuple->column(prefix + "_RECONSTRUCTED",   int(m_reconstructed->reconstructed(p)));
  ntuple->column(prefix + "_RECONSTRUCTIBLE", int(m_reconstructible->reconstructible(p)));
 
  return;
}
/**************************************************************************** *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * Simulated particle polarization.                                           *
 *                                                                            *
 **************************************************************************** */
void
DoubleJpsiSel::FillMCPolarization(Tuples::Tuple          &ntuple,
                                 const LHCb::MCParticle *mot,
                                 const LHCb::MCParticle *dau,
                                 const std::string       prefix){
  if (!mot) return;
  if (!dau) return;
 
  Gaudi::LorentzVector mot4Mom = mot->momentum();
  Gaudi::LorentzVector dau4Mom = dau->momentum();
 
  ROOT::Math::Boost boost(mot4Mom.BoostToCM());
 
  const Gaudi::XYZVector dau4MomBoost = (boost(dau4Mom)).Vect().unit();
  const Gaudi::XYZVector Mot4MomBoost = mot4Mom.Vect().unit();
 
  double cosTheta = fabs(dau4MomBoost.Dot(Mot4MomBoost));
 
  ntuple->column(prefix + "_COSTHETA", cosTheta);
 
  return;
}

/**************************************************************************** *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * Information on the ancestors of a MC particles.                            *
 *                                                                            *
 **************************************************************************** */
void
DoubleJpsiSel::FillMCAncestor(Tuples::Tuple          &ntuple,
                             const LHCb::MCParticle *p,
                             const std::string       prefix){


  int pdgPid = -1;

  if (p){
    pdgPid = 0;
    const LHCb::MCParticle *mother = p->mother();
    while (mother != NULL){
      pdgPid = mother->particleID().pid();
      mother = mother->mother();
    }
  }

  ntuple->column(prefix + "_TOP_ANCESTOR_PID", pdgPid);

  return;
}


void DoubleJpsiSel::FillJpsiNtuple(Tuples::Tuple        &ntuple,
				   const LHCb::Particle *jpsi,
				   const std::string     prefix){
  ntuple->column(prefix + "_P",            jpsi->p());
  ntuple->column(prefix + "_PT",           jpsi->pt());
  ntuple->column(prefix + "_RAPIDITY",     jpsi->momentum().Rapidity());
  ntuple->column(prefix + "_MASS",         jpsi->measuredMass());
  ntuple->column(prefix + "_VCHI2",        jpsi->endVertex()->chi2());
  ntuple->column(prefix + "_VNDOF",        jpsi->endVertex()->nDoF());
  ntuple->column(prefix + "_VCHI2PERNDOF", jpsi->endVertex()->chi2PerDoF());
  if (m_verbose){
    ntuple->column(prefix + "_PID",        jpsi->particleID().pid());
  }
  if (1) return;
  

  // L0 Trigger TISTOS

  bool tis = false;
  bool tos = false;
  bool dec = false;

  bool l0tis = false;
  bool l0tos = false;

  m_L0TriggerTisTosTool->setOfflineInput(*jpsi);

  std::vector<std::string>::iterator s;

  for (s = m_L0TISList.begin(); s != m_L0TISList.end(); ++s){
    m_L0TriggerTisTosTool->triggerTisTos(*s, dec, tis, tos);
    l0tis = l0tis || tis;
  }

  for (s = m_L0TOSList.begin(); s != m_L0TOSList.end(); ++s){
    m_L0TriggerTisTosTool->triggerTisTos(*s, dec, tis, tos);
    l0tos = l0tos || tos;
  }

  ntuple->column(prefix + "_L0TIS", l0tis);
  ntuple->column(prefix + "_L0TOS", l0tos);

  // HLT1 Trigger TISTOS

  bool hlt1tis = false;
  bool hlt1tos = false;

  m_TriggerTisTosTool->setOfflineInput(*jpsi);

  for (s = m_HLT1TISList.begin(); s != m_HLT1TISList.end(); ++s){
    m_TriggerTisTosTool->triggerTisTos(*s, dec, tis, tos);
    hlt1tis = hlt1tis || tis;
  }

  for (s = m_HLT1TOSList.begin(); s != m_HLT1TOSList.end(); ++s){
    m_TriggerTisTosTool->triggerTisTos(*s, dec, tis, tos);
    hlt1tos = hlt1tos || tos;
  }

  ntuple->column(prefix + "_HLT1TIS", hlt1tis);
  ntuple->column(prefix + "_HLT1TOS", hlt1tos);

  // HLT2 Trigger TISTOS

  bool hlt2tis = false;
  bool hlt2tos = false;

  for (s = m_HLT2TISList.begin(); s != m_HLT2TISList.end(); ++s){
    m_TriggerTisTosTool->triggerTisTos(*s, dec, tis, tos);
    hlt2tis = hlt2tis || tis;
  }

  for (s = m_HLT2TOSList.begin(); s != m_HLT2TOSList.end(); ++s){
    m_TriggerTisTosTool->triggerTisTos(*s, dec, tis, tos);
    hlt2tos = hlt2tos || tos;
  }

  ntuple->column(prefix + "_HLT2TIS", hlt2tis);
  ntuple->column(prefix + "_HLT2TOS", hlt2tos);

  ntuple->column(prefix + "_TIS", hlt1tis && hlt2tis && l0tis);
  ntuple->column(prefix + "_TOS", hlt1tos && hlt2tos && l0tos);

  return;
}

void DoubleJpsiSel::FillMuNtuple(Tuples::Tuple        &ntuple,
				 const LHCb::Particle *muon,
				 const std::string     prefix){

  const LHCb::ProtoParticle *proto = muon->proto();
  const LHCb::Track         *track = proto->track();
  
  ntuple->column(prefix + "_P",           muon->p());
  ntuple->column(prefix + "_PT",          muon->pt());
  ntuple->column(prefix + "_RAPIDITY",    muon->momentum().Rapidity());
  ntuple->column(prefix + "_CHI2",        track->chi2());
  ntuple->column(prefix + "_NDOF",        track->nDoF());
  ntuple->column(prefix + "_CHI2PERNDOF", track->chi2PerDoF());
  ntuple->column(prefix + "_KL",          track->info(101, -1));
  ntuple->column(prefix + "_DLLMu",       proto->info(601, -1));

  if (m_verbose){
    ntuple->column(prefix + "_DLLe",        proto->info(600, -1));
    ntuple->column(prefix + "_DLLK",        proto->info(603, -1));
    ntuple->column(prefix + "_DLLp",        proto->info(604, -1));
    ntuple->column(prefix + "_DLLMu_e",     proto->info(601, -1) -
		   proto->info(600, -1));
    ntuple->column(prefix + "_DLLMu_K",     proto->info(601, -1) -
		   proto->info(603, -1));
    ntuple->column(prefix + "_DLLMu_p",     proto->info(601, -1) -
		   proto->info(604, -1));

    const LHCb::MuonPID *muPid = proto->muonPID();
    //  ntuple->column(prefix + "_ISMUON", muPid->IsMuon());
    
    // Impact parameter
    
    double minIp     = -1.0;
    double minIpChi2 = -1.0;
    ImpactParameter(muon, &minIp, &minIpChi2);
    
    ntuple->column(prefix + "_MINIP",     minIp);
    ntuple->column(prefix + "_MINIPCHI2", minIpChi2);
  }

  return;
}

void DoubleJpsiSel::FillMotherNtuple(
	   Tuples::Tuple                     &ntuple,
	   const LHCb::Particle              *mother,
	   const Gaudi::Math::ParticleParams *fitParams,
	   const std::string                   prefix){

  Gaudi::Math::ValueWithError ctau        = fitParams->ctau();
  Gaudi::Math::ValueWithError decayLength = fitParams->decayLength();

  ntuple->column(prefix + "_CHI2",           mother->endVertex()->chi2());
  ntuple->column(prefix + "_NDOF",           mother->endVertex()->nDoF()); 
  ntuple->column(prefix + "_CHI2PERNDOF",    mother->endVertex()->chi2PerDoF());

  if (m_verbose){
    ntuple->column(prefix + "_CTAU",           ctau.value());
    ntuple->column(prefix + "_CTAUERR",        ctau.error());
    ntuple->column(prefix + "_DECAYLENGTH",    decayLength.value());
    ntuple->column(prefix + "_DECAYLENGTHERR", decayLength.error());
    
    double minIp     = -1.0;
    double minIpChi2 = -1.0;
    ImpactParameter(mother, &minIp, &minIpChi2);
    
    ntuple->column(prefix + "_MINIP",     minIp);
    ntuple->column(prefix + "_MINIPCHI2", minIpChi2);
    
    const LHCb::VertexBase *pv = bestPV(mother);
    Gaudi::XYZVector        dl = pv->position() - mother->endVertex()->position();
    
    ntuple->column(prefix + "_DL", dl.R());
  }

  return;
}



const LHCb::Particle *
DoubleJpsiSel::GetParticles(const LHCb::Particle *p,
                           int charge,
                           int pid){
  LHCb::Particle::ConstVector allParticles = p->daughtersVector();
  LHCb::Particle::ConstVector chargedParticles;
  LHCb::Particle::ConstVector foundParticles;
  LHCb::Particle::ConstVector otherParticles;
 
  //
  // Select by charge
  //
  

  size_t nParticles = 0;
  if (charge > 0){
    nParticles = DaVinci::filter(allParticles,
                                 bind(&LHCb::Particle::charge,_1) > 0,
                                 chargedParticles);
  } else if (charge < 0){
    nParticles = DaVinci::filter(allParticles,
                                 bind(&LHCb::Particle::charge,_1) < 0,
                                 chargedParticles);
  } else{
    nParticles = DaVinci::filter(allParticles,
                                 bind(&LHCb::Particle::charge,_1) == 0,
                                 chargedParticles);
  }
 
  if (nParticles < 1) {
    warning() << "Number of particle with charge " << charge
              << " = " << nParticles
              << endmsg;
    return NULL;
  }
 
  //
  // Select by particle-id
  //
 
  LHCb::ParticleID pID(pid);
  nParticles = DaVinci::filter(chargedParticles,
                               bind(&LHCb::Particle::particleID,_1) == pID,
                               foundParticles, otherParticles);
 
  if (nParticles != 1) {
    warning() << "Number of particle with charge " << charge
              << " and pid " << pid
              << " = "       << foundParticles.size()
              << " of "      << chargedParticles.size()
              << endmsg;
    return NULL;
  }
 
  return *(foundParticles.begin());
}
 



void DoubleJpsiSel::FillPolarizationNtuple(Tuples::Tuple &ntuple,
					   const LHCb::Particle *jpsi,
					   const LHCb::Particle *muon,
					   const std::string     prefix){

  
  Gaudi::LorentzVector jpsi4Mom = jpsi->momentum();
  Gaudi::LorentzVector mu4Mom   = muon->momentum();

  ROOT::Math::Boost boost(jpsi4Mom.BoostToCM());

  const Gaudi::XYZVector boostedMu4Mom   = (boost(mu4Mom)).Vect().unit();
  const Gaudi::XYZVector boostedJpsi4Mom = jpsi4Mom.Vect().unit();

  double cosTheta = fabs(boostedMu4Mom.Dot(boostedJpsi4Mom));

  ntuple->column(prefix + "_COSTHETA", cosTheta);
  
  return;
}

 
const LHCb::Particle *DoubleJpsiSel::GetMuon(const LHCb::Particle *jpsi,
					     int charge){
  LHCb::Particle::ConstVector chargedMuons;
  LHCb::Particle::ConstVector allMuons = jpsi->daughtersVector();

  charge = (charge >= 0) ? 1 : -1;

  size_t nMuons = (charge > 0) ?
    DaVinci::filter(allMuons, bind(&LHCb::Particle::charge,_1) < 0,
		    chargedMuons)  :
    DaVinci::filter(allMuons, bind(&LHCb::Particle::charge,_1) > 0,
		    chargedMuons);

  if (nMuons != 1){
    warning() << "Number of mu with charge " << charge
	      << " = " << nMuons
	      << endmsg;
    return NULL;
  }

  return *(chargedMuons.begin());
}

void DoubleJpsiSel::ImpactParameter(
	   const LHCb::Particle *particle,
	   double               *minIp,
	   double               *minIpChi2){
  const LHCb::RecVertex::Range pv = m_dva->primaryVertices();

  if (!pv.empty()){
    LHCb::RecVertex::Range::const_iterator ipv;
    for (ipv = pv.begin(); ipv != pv.end(); ++ipv){
      LHCb::RecVertex newPv(**ipv);

      double ip   = 0;
      double chi2 = 0;

      LHCb::VertexBase* newPvPtr = (LHCb::VertexBase*) &newPv; 
      StatusCode scIp = m_dist->distance(particle, newPvPtr, ip, chi2);

      if (scIp){
	if ((ip   < *minIp)     || (*minIp     < 0)) *minIp     = ip;
	if ((chi2 < *minIpChi2) || (*minIpChi2 < 0)) *minIpChi2 = chi2;
      }
    }
  }

  return;
}
 

void
DoubleJpsiSel::FillAnalysisDoubleJpsiMuMu(Tuples::Tuple &ntuple,
					 const LHCb::Particle *jpsi1,
					 const LHCb::Particle *muPlus1,
					 const LHCb::Particle *muMinus1,
					 const LHCb::Particle *jpsi2,
					 const LHCb::Particle *muPlus2,
					 const LHCb::Particle *muMinus2,
					 const LHCb::Particle *mother){
  const std::string prefixJpsi1    = "Jpsi_1";
  const std::string prefixJpsi2    = "Jpsi_2";
  const std::string prefixMother   = "DoubleJpsi";
  const std::string prefixMuPlus1  = "MuPlus1";
  const std::string prefixMuPlus2  = "MuPlus2";
  const std::string prefixMuMinus1 = "MuMinus1";
  const std::string prefixMuMinus2 = "MuMinus2";

  const LHCb::Track *trackPlus1  = muPlus1->proto()->track();  
  const LHCb::Track *trackMinus1 = muMinus1->proto()->track();
  const LHCb::Track *trackPlus2  = muPlus2->proto()->track();  
  const LHCb::Track *trackMinus2 = muMinus2->proto()->track();

  FillDataEventNtuple(ntuple);
  FillDataParticle(ntuple,      jpsi1,      prefixJpsi1);
  
  FillDataParticle(ntuple,       muPlus1,    prefixMuPlus1);
  FillDataPid(ntuple,            muPlus1,    prefixMuPlus1);
  FillDataTrackGhostProb(ntuple, trackPlus1, prefixMuPlus1);
  
  //FillDataMCTruth(ntuple, jpsi1, prefixJpsi1);
  //FillDataMCTruth(ntuple, jpsi2, prefixJpsi2);
 
  FillDataParticle(ntuple,       muMinus1,    prefixMuMinus1);
  FillDataPid(ntuple,            muMinus1,    prefixMuMinus1);
  FillDataTrackGhostProb(ntuple, trackMinus1, prefixMuMinus1);
  
  FillDataParticle(ntuple,      jpsi2,       prefixJpsi2);
  
  FillDataParticle(ntuple,       muPlus2,     prefixMuPlus2);
  FillDataPid(ntuple,            muPlus2,     prefixMuPlus2);
  FillDataTrackGhostProb(ntuple, trackPlus2,  prefixMuPlus2);
  
  FillDataParticle(ntuple,       muMinus2,    prefixMuMinus2);
  FillDataPid(ntuple,            muMinus2,    prefixMuMinus2);
  FillDataTrackGhostProb(ntuple, trackMinus2, prefixMuMinus2);
  
  //  if (m_fillTrigger){
  // FillTrigger(ntuple, jpsi1, prefixJpsi1);
  // FillTrigger(ntuple, jpsi2, prefixJpsi2);
  //}
    
  FillDataPolarizationNtuple(ntuple, jpsi1, muPlus1, prefixJpsi1);
  FillDataPolarizationNtuple(ntuple, jpsi2, muPlus2, prefixJpsi2);
  FillDataPropertime(ntuple, jpsi1, prefixJpsi1);
  FillDataPropertime(ntuple, jpsi2, prefixJpsi2);
  // FillDataParticle(ntuple, mother, prefixMother);
  //FillDataTransverseAngleNtuple(ntuple, jpsi1, jpsi2, prefixMother);

  return;
}


/**************************************************************************** *
 *                                                                            *
 * FillAnalysisDoubleJpsiMuMu:                                                *
 *                                                                            *
 * Information about two J/psi in the same event and their                    *
 * mu mu decay products.                                                      *
 *                                                                            *
 **************************************************************************** */
void
DoubleJpsiSel::FillAnalysisDoubleJpsiMuMu(Tuples::Tuple &ntuple,
                                          const LHCb::Particle *jpsi1,
                                          const LHCb::Particle *muPlus1,
                                          const LHCb::Particle *muMinus1,
                                          const LHCb::Particle *jpsi2,
                                          const LHCb::Particle *muPlus2,
                                          const LHCb::Particle *muMinus2,
                                          const LHCb::Particle *fourMuon,
                                          const LHCb::Particle *mother,
                                          const LHCb::Particle *mother1,
                                          const LHCb::Particle *mother2
                                          ){
  const std::string prefixFourMuon = "FourMuons";
  const std::string prefixJpsi1    = "JPsi_1";
  const std::string prefixJpsi2    = "JPsi_2";
  const std::string prefixMother   = "DoubleJpsi";
  const std::string prefixMuPlus1  = "MuPlus1";
  const std::string prefixMuPlus2  = "MuPlus2";
  const std::string prefixMuMinus1 = "MuMinus1";
  const std::string prefixMuMinus2 = "MuMinus2";
 
  const LHCb::Track *trackPlus1  = muPlus1->proto()->track();  
  const LHCb::Track *trackMinus1 = muMinus1->proto()->track();
  const LHCb::Track *trackPlus2  = muPlus2->proto()->track();  
  const LHCb::Track *trackMinus2 = muMinus2->proto()->track();
 
  //FillDataEventNtuple(ntuple); //AC and comments out program it calls
  FillDataComposite(ntuple,      jpsi1,      prefixJpsi1);
  FillDataPropertime(ntuple, jpsi1, prefixJpsi1);
 
  FillDataParticle(ntuple,       muPlus1,    prefixMuPlus1);
  FillDataPid(ntuple,            muPlus1,    prefixMuPlus1);
  FillDataTrackGhostProb(ntuple, trackPlus1, prefixMuPlus1);

  FillDataImpactParameter(ntuple,       muPlus1,    prefixMuPlus1);
  FillDataImpactParameter(ntuple,       muMinus1,    prefixMuMinus2);

  // FillDataImpactParameter(ntuple,       jpsi1,    prefixJpsi1);
  //FillDataImpactParameter(ntuple,       jpsi2,    prefixJpsi2);
  FillDataParticle(ntuple,       muMinus1,    prefixMuMinus1);
  FillDataPid(ntuple,            muMinus1,    prefixMuMinus1);
  FillDataTrackGhostProb(ntuple, trackMinus1, prefixMuMinus1);
 
  FillDataComposite(ntuple,      jpsi2,       prefixJpsi2);

  FillDataPropertime(ntuple, jpsi2, prefixJpsi2);
 
  FillDataParticle(ntuple,       muPlus2,     prefixMuPlus2);
  FillDataPid(ntuple,            muPlus2,     prefixMuPlus2);
  FillDataTrackGhostProb(ntuple, trackPlus2,  prefixMuPlus2);
 
  FillDataParticle(ntuple,       muMinus2,    prefixMuMinus2);
  FillDataPid(ntuple,            muMinus2,    prefixMuMinus2);
  FillDataTrackGhostProb(ntuple, trackMinus2, prefixMuMinus2);
 
  FillDataImpactParameter(ntuple,       muPlus2,    prefixMuPlus2);
  FillDataImpactParameter(ntuple,       muMinus2,    prefixMuMinus2);
  
  FillDataEventNtuple(ntuple); //best track
  FillDataDTF(ntuple, mother1->endVertex(),prefixJpsi1); //jpsi 1 DTF
  FillDataDTF(ntuple, mother2->endVertex(), prefixJpsi2); //jpsi 2 DTF

  //if (m_fillTrigger){
  //  FillTrigger(ntuple, jpsi1, prefixJpsi1);
  //  FillTrigger(ntuple, jpsi2, prefixJpsi2);
  // }
 
  if (m_verbose){
    //    FillDataFourMuonNtuple(ntuple, fourMuon, prefixFourMuon);
    FillDataComposite(ntuple,    fourMuon, prefixFourMuon);
    //    FillDataMotherNtuple(ntuple, mother,   prefixMother);
  }
 
  FillDataPolarizationNtuple(ntuple, jpsi1, muPlus1, prefixJpsi1);
  FillDataPolarizationNtuple(ntuple, jpsi2, muPlus2, prefixJpsi2);
  FillDataComposite(ntuple, mother, prefixMother);
  FillDataTransverseAngleNtuple(ntuple, jpsi1, jpsi2, prefixMother);
  FillDataDTF(ntuple, mother->endVertex(), prefixMother);
  return;
}




/**************************************************************************** *
 *                                                                            *
 * SignalSelect::FillDataTransverseAngleNtuple                                *
 *                                                                            *
 * The angle in the transverse plane between two particles.                   *
 *                                                                            *
 **************************************************************************** */
void
DoubleJpsiSel::FillDataTransverseAngleNtuple(Tuples::Tuple        &ntuple,
                                            const LHCb::Particle *part1,
                                            const LHCb::Particle *part2,
                                            const std::string     prefix){
  Gaudi::LorentzVector p1 = part1->momentum();
  Gaudi::LorentzVector p2 = part2->momentum();
 
  double pp1      = sqrt(p1.X()*p1.X() + p1.Y()*p1.Y());
  double pp2      = sqrt(p2.X()*p2.X() + p2.Y()*p2.Y());
  double cosPhi   = (p1.X()*p2.X() + p1.Y()*p2.Y())/(pp1*pp2);
  double deltaPhi = acos(cosPhi);
 
  ntuple->column(prefix + "_DELTA_PHI", deltaPhi);
  return;
}


void
DoubleJpsiSel::FillDataPolarizationNtuple(Tuples::Tuple        &ntuple,
                                         const LHCb::Particle *mother,
                                         const LHCb::Particle *daughter,
                                         const std::string     prefix){
  Gaudi::LorentzVector mother4Mom   = mother->momentum();
  Gaudi::LorentzVector daughter4Mom = daughter->momentum();
 
  ROOT::Math::Boost boost(mother4Mom.BoostToCM());
 
  const Gaudi::XYZVector boostedDaughter4Mom = (boost(daughter4Mom)).Vect().unit();
  const Gaudi::XYZVector boostedMother4Mom   = mother4Mom.Vect().unit();
 
  double cosTheta = fabs(boostedDaughter4Mom.Dot(boostedMother4Mom));
 
  ntuple->column(prefix + "_COSTHETA", cosTheta);
 
  return;
}



void
DoubleJpsiSel::FillDataEventNtuple(Tuples::Tuple &ntuple){
  // TCK and beam information
  

  LHCb::ODIN *odin(NULL);
 
  if (exist<LHCb::ODIN>(LHCb::ODINLocation::Default)){
    odin = get<LHCb::ODIN>(LHCb::ODINLocation::Default);
  } else if (exist<LHCb::ODIN>(LHCb::ODINLocation::Default, false)){
    odin = get<LHCb::ODIN>(LHCb::ODINLocation::Default, false);
  } else {
    warning() << "Cannot load the ODIN data object" << endreq;
  }
 
  const std::string prefix = "Event";
 
  if (m_magnetSvc) ntuple->column(prefix + "_Polarity", m_magnetSvc->isDown() ? -1 : 1);
 
  m_tck = odin ? int(odin->triggerConfigurationKey()) : -1;
  ntuple->column(prefix + "_TCK", m_tck);
 
  ntuple->column(prefix + "_RunNumber",
                 odin ? odin->runNumber()         : -1);
  ntuple->column(prefix + "_EventNumber",
                 odin ? odin->eventNumber()       : -1);
  ntuple->column(prefix + "_BCID",
                 odin ? odin->bunchId()           : -1);
  ntuple->column(prefix + "_BCType",
                 odin ? odin->bunchCrossingType() : -1);
  //
  // Reconstruction information.
  //
 
  std::string GECLoc = "GlobalEventCounters";
 
  const Gaudi::Numbers* numbers = NULL;
  if (exist<Gaudi::Numbers>(GECLoc))
    numbers = get<Gaudi::Numbers>(GECLoc);
 
  const LHCb::RecSummary *recSummary = NULL;
  if (exist<LHCb::RecSummary>(LHCb::RecSummaryLocation::Default)){
    recSummary = get<LHCb::RecSummary>(LHCb::RecSummaryLocation::Default);
  } else if (exist<LHCb::RecSummary>(LHCb::RecSummaryLocation::Default,   false)){
    recSummary = get<LHCb::RecSummary>(LHCb::RecSummaryLocation::Default, false);
  }
 
  if (numbers){
    const Gaudi::Numbers::Map &m = numbers->numbers();
    ntuple->column(prefix + "_BestTracksMult", int(m["nLong"]));
    ntuple->column(prefix + "_nSPDHits",       int(m["nSpd"]));
   
    if (m_verbose){
      ntuple->column(prefix + "_VeloTracksMult", int(m["nVelo"]));
      ntuple->column(prefix + "_nOTClusters",    int(m["nOT"]));
      ntuple->column(prefix + "_nITClusters",    int(m["nITClusters"]));
      ntuple->column(prefix + "_nTTClusters",    int(m["nTTClusters"]));
    }
  } else if (recSummary){
    ntuple->column(prefix + "_BestTracksMult",
                   int(recSummary->info(LHCb::RecSummary::nTracks,  -1)));
    ntuple->column(prefix + "_nSPDHits",
                   int(recSummary->info(LHCb::RecSummary::nSPDhits, -1)));
   
    if (m_verbose){
      ntuple->column(prefix + "_VeloTracksMult",
                     int(recSummary->info(LHCb::RecSummary::nVeloTracks, -1)));
      ntuple->column(prefix + "_nOTClusters",
                     int(recSummary->info(LHCb::RecSummary::nOTClusters, -1)));
      ntuple->column(prefix + "_nITClusters",
                     int(recSummary->info(LHCb::RecSummary::nITClusters, -1)));
      ntuple->column(prefix + "_nTTClusters",
                     int(recSummary->info(LHCb::RecSummary::nTTClusters, -1)));
    }
  } else{
    //    warning() << "Cannot load neither the RecSummary data object "
    //        << "nor the Global Event Counters."
    //        << endreq;
  }
 
  //
  // Primary vertex multiplicity
  //
 
  //  const LHCb::RecVertex::Range pv = m_dva->primaryVertices();
  //  ntuple->column(prefix + "_nPV", int(pv.size()));
 
  if (m_inputType == "DST"){
    std::string pvLocation = LHCb::RecVertexLocation::Primary;
    const LHCb::RecVertex::Container *pv = getIfExists<LHCb::RecVertex::Container>(pvLocation);
    if (pv)
      ntuple->column(prefix + "_nPV", int(pv->size()));
  } else if (m_inputType == "MDST"){
    const LHCb::RecVertex::Range pv1 = m_dva->primaryVertices();
    ntuple->column(prefix + "_nPV", int(pv1.size()));
  }  

  return;
}


/**************************************************************************** *
 *                                                                            *
 * SignalSelect::FillDataParticle                                             *
 *                                                                            *
 * Information of basic particles in real data                                #
 *                                                                            *
 **************************************************************************** */
void
DoubleJpsiSel::FillDataParticle(Tuples::Tuple        &ntuple,
                               const LHCb::Particle *p,
                               const std::string     prefix){
 
  const LHCb::ProtoParticle *proto = p->proto();
  const LHCb::Track         *track = proto->track();
 
  ntuple->column(prefix + "_P",              double(p->p()));
  ntuple->column(prefix + "_PT",             double(p->pt()));
  ntuple->column(prefix + "_RAPIDITY",       double(p->momentum().Rapidity()));
  ntuple->column(prefix + "_PSEUDORAPIDITY", double(track->pseudoRapidity()));
  ntuple->column(prefix + "_CHI2",           double(track->chi2()));
  ntuple->column(prefix + "_NDOF",           int(track->nDoF()));
  ntuple->column(prefix + "_CHI2PERNDOF",    double(track->chi2PerDoF()));
  ntuple->column(prefix + "_KL",             double(track->info(101, -1)));
  ntuple->column(prefix + "_CHARGE",         int(p->charge()));
 

  return;
}


void
DoubleJpsiSel::FillDataImpactParameter(Tuples::Tuple        &ntuple,
                               const LHCb::Particle *p,
                               const std::string     prefix){
 
  const LHCb::ProtoParticle *proto = p->proto();
  const LHCb::Track         *track = proto->track();
 
  
 
 
    // Impact parameter
   
    double minIp     = -1.0;
    double minIpChi2 = -1.0;
    ImpactParameter(p, &minIp, &minIpChi2);
   
    ntuple->column(prefix + "_MINIP",     minIp);
    ntuple->column(prefix + "_MINIPCHI2", minIpChi2);
  
 
  return;
}


/* *************************************************************************** *
 *                                                                             *
 * FillDataComposite                                                           *
 *                                                                             *
 * Information of composite particles in real data                             *
 *                                                                             *
 ***************************************************************************** */
void
DoubleJpsiSel::FillDataComposite(Tuples::Tuple        &ntuple,
                                const LHCb::Particle *p,
                                const std::string     prefix){
  ntuple->column(prefix + "_P",            double(p->p()));
  ntuple->column(prefix + "_PT",           double(p->pt()));
  ntuple->column(prefix + "_RAPIDITY",     double(p->momentum().Rapidity()));
  ntuple->column(prefix + "_MASS",         double(p->measuredMass()));
  ntuple->column(prefix + "_VCHI2",        double(p->endVertex()->chi2()));
  ntuple->column(prefix + "_VNDOF",        double(p->endVertex()->nDoF()));
  ntuple->column(prefix + "_VCHI2PERNDOF", double(p->endVertex()->chi2PerDoF()));
  ntuple->column(prefix + "_PID",          double(p->particleID().pid()));
 
  //
  // Background category when running over simulation.
  //
 
  
  if ((true) && (!m_MCTruth)){
    IBackgroundCategory::categories cat = IBackgroundCategory::Undefined;
 
    std::vector<IBackgroundCategory*>::const_iterator it;
    for (it = m_bkgs.begin(); it != m_bkgs.end(); ++it){
      cat = (*it)->category(p);
      if (cat != IBackgroundCategory::Undefined) break;
    }
     
    ntuple->column( prefix + "_BKGCAT", int(cat));
  }
 
  const LHCb::VertexBase *pv = bestPV(p);
  Gaudi::XYZVector        dl = pv->position() - p->endVertex()->position(); 
  ntuple->column(prefix + "_DL", double(dl.R()));


  return;
}


/**************************************************************************** *
 *                                                                            *
 * SignalSelect::FillDataTrackGhostProb                                       *
 *                                                                            *
 * Information of tracks in real data                                         #
 *                                                                            *
 **************************************************************************** */
    void
      DoubleJpsiSel::FillDataTrackGhostProb(Tuples::Tuple     &ntuple,
                                           const LHCb::Track *track,
                                           const std::string  prefix){
      ntuple->column(prefix + "_TRACK_GHOSTPROB", track->ghostProbability());
      ntuple->column(prefix + "_Pseudorapidity", track->pseudoRapidity());
      return;
    }


/* *************************************************************************** *
 *                                                                             *
 * FillDataDTF                                                                 *
 *                                                                             *
 * Information of the DecayTreeFitter results.                                 *
 *                                                                             *
 ***************************************************************************** */
void
DoubleJpsiSel::FillDataDTF(Tuples::Tuple      &ntuple,
			  const LHCb::Vertex *v,
			  const std::string   prefix){
  ntuple->column(prefix + "_DTF_NDOF",        v->nDoF());
  ntuple->column(prefix + "_DTF_CHI2",        double(v->chi2()));
  ntuple->column(prefix + "_DTF_CHI2PERNDOF", double(v->chi2PerDoF()));

  return;
}



/**************************************************************************** *
 *                                                                            *
 * FillDataMCTruth                                              *
 *                                                                            *
 * Information about the MC truth in simulation                               *
 *                                                                            *
 **************************************************************************** */
void
DoubleJpsiSel::FillDataMCTruth(Tuples::Tuple        &ntuple,
                              const LHCb::Particle *p,
                              const std::string     prefix){
 

  if (!p) return;

  const LHCb::MCParticle *mcp = m_associator->relatedMCP(p);

  if (!mcp){
 debug() << "No related mcp found" << endmsg;
 //return;
  }
  
  int pid = (mcp) ? int(mcp->particleID().pid()) : -1;
  ntuple->column(prefix + "_TRUE_PID", pid);

  FillMCAncestor(ntuple, mcp, prefix);



  bool isMuMuTrue = false;

  if(mcp)
  {
   const SmartRefVector<LHCb::MCVertex> &vertices = mcp->endVertices();
   SmartRefVector<LHCb::MCVertex>::const_iterator iVert;
   for (iVert = vertices.begin(); iVert != vertices.end(); ++iVert){
     
      const LHCb::MCParticle *muPlus  = NULL;
      const LHCb::MCParticle *muMinus = NULL;
      
      const SmartRefVector<LHCb::MCParticle> &daughters = (*iVert)->products();
      SmartRefVector<LHCb::MCParticle>::const_iterator iDau;
      for (iDau = daughters.begin(); iDau != daughters.end(); ++iDau){
        if ((*iDau)->particleID().pid() ==  13) muPlus  = *iDau;
        if ((*iDau)->particleID().pid() == -13) muMinus = *iDau;
      }
      if(!muPlus || !muMinus) isMuMuTrue = false;
      else isMuMuTrue = true;
   }
  }
  
  ntuple->column("isMuMuMCTrue", isMuMuTrue);


  /*

  //---AC - check MCTruth info

  bool isJpsiTrue = false;
  bool isMuMuTrue = false;
   if(mcp->particleID().pid()!= 443){;
      isJpsiTrue = true;
    }

  const SmartRefVector<LHCb::MCVertex> &vertices = mcp->endVertices();
    SmartRefVector<LHCb::MCVertex>::const_iterator iVert;
    for (iVert = vertices.begin(); iVert != vertices.end(); ++iVert){
      
      const LHCb::MCParticle *muPlus  = NULL;
      const LHCb::MCParticle *muMinus = NULL;

      const SmartRefVector<LHCb::MCParticle> &daughters = (*iVert)->products();
      SmartRefVector<LHCb::MCParticle>::const_iterator iDau;
      for (iDau = daughters.begin(); iDau != daughters.end(); ++iDau){
        if ((*iDau)->particleID().pid() ==  13) muPlus  = *iDau;
        if ((*iDau)->particleID().pid() == -13) muMinus = *iDau;
      }
      
      //if (!muPlus)  continue;
      //if (!muMinus) continue;
      if(!muPlus && !muMinus) isMuMuTrue = false;
      else isMuMuTrue = true;
      ntuple->column("isJpsiMCTrue", isJpsiTrue);
      ntuple->column("isMuMuMCTrue", isMuMuTrue);
    }
    
  */  

  ///
       
  return;
}



/* *************************************************************************** *
 *                                                                             *
 * FillDataMCTruth                                                             *
 *                                                                             *
 * Information of the MCTruthChecks                                            *
 *                                                                             *
 ***************************************************************************** */
void
DoubleJpsiSel::FillDataMCTruth(Tuples::Tuple      &ntuple,
			  const double motherPID,
			  const std::string   prefix){
  ntuple->column(prefix + "TruthInfo",        motherPID);
  return;
}



/**************************************************************************** *
 *                                                                            *
 * SignalSelect::FillDataPid                                                  *
 *                                                                            *
 * Information on the pid of a particle.                                      #
 *                                                                            *
 **************************************************************************** */
void
DoubleJpsiSel::FillDataPid(Tuples::Tuple        &ntuple,
			  const LHCb::Particle *p,
			  const std::string     prefix){

  int dllMask = (m_verbose) ? 0x3F : m_selFillDll;
  const LHCb::ProtoParticle *proto = p->proto();

  if (dllMask & 0x1){
    ntuple->column(prefix + "_DLLMu",
		   double(proto->info(LHCb::ProtoParticle::CombDLLmu,   -1000)));
    ntuple->column(prefix + "_ProbNNMu",
		   double(proto->info(LHCb::ProtoParticle::ProbNNmu,    -1000)));
  }
  if (dllMask & 0x2){
    ntuple->column(prefix + "_DLLe",
		   double(proto->info(LHCb::ProtoParticle::CombDLLe,    -1000)));
    ntuple->column(prefix + "_ProbNNe",
		   double(proto->info(LHCb::ProtoParticle::ProbNNe,     -1000)));
  }
  if (dllMask & 0x4){
    ntuple->column(prefix + "_DLLK",
		   double(proto->info(LHCb::ProtoParticle::CombDLLk,    -1000)));
    ntuple->column(prefix + "_ProbNNK",
		   double(proto->info(LHCb::ProtoParticle::ProbNNk,     -1000)));
  }
  if (dllMask & 0x8){
    ntuple->column(prefix + "_DLLp",
		   double(proto->info(LHCb::ProtoParticle::CombDLLp,    -1000)));
    ntuple->column(prefix + "_ProbNNp",
		   double(proto->info(LHCb::ProtoParticle::ProbNNp,     -1000)));
  }
  if (dllMask & 0x10){
    ntuple->column(prefix + "_ProbNNpi",
		   double(proto->info(LHCb::ProtoParticle::ProbNNpi,    -1000)));
  }
  if (dllMask & 0x20){
    ntuple->column(prefix + "_ProbNNghost",
		   double(proto->info(LHCb::ProtoParticle::ProbNNghost, -1000)));
  }
  return;
}




//===========================================================================
//loop for 2010 MC
//===========================================================================

StatusCode DoubleJpsiSel::MCTenLoop(const LHCb::Particle::Range &muons){
  LHCb::Particle::ConstVector JpsiVector;
  counter("A. start 2010 loop")++;

  StatusCode sc = StatusCode::FAILURE;

  LHCb::Particle::ConstVector MuPlus, MuMinus;

  size_t nMuons = DaVinci::filter(muons, bind(&LHCb::Particle::charge,_1)<0, MuMinus);
  if (nMuons>1) {
    nMuons += DaVinci::filter(muons, bind(&LHCb::Particle::charge,_1)>0, MuPlus);
  }

  //
  //Loop on the muons positive and negative
  //
  for (LHCb::Particle::ConstVector::const_iterator imp = MuPlus.begin() ;
       imp != MuPlus.end() ; ++imp){
    for ( LHCb::Particle::ConstVector::const_iterator imm =  MuMinus.begin() ;
          imm != MuMinus.end() ; ++imm ){
      //  m_vertextFitter = getTool<IParticleReFitter>;

      counter("B  muon loop 2010 loop")++;
      const IParticleCombiner *combiner = particleCombiner();
      
      LHCb::Vertex MuMuVertex; 
      LHCb::Particle *Jpsi1 = new LHCb::Particle;
      Jpsi1->setParticleID( m_jPsiID );
    
      LHCb::Particle::ConstVector particlesForFit1;
      const LHCb::Particle *muPlus1  = *imp;
      const LHCb::Particle *muMinus1 = *imm;
 
      particlesForFit1.push_back(muPlus1);
      particlesForFit1.push_back(muMinus1);
   
      // StatusCode scFit1 = vertexFitter()->fit( *(*imp), *(*imm), MuMuVertex, Jpsi );
      // StatusCode scFit1 = particleCombiner()->fit( *(*imp), *(*imm), MuMuVertex, Jpsi );
      StatusCode scFit1 = combiner->combine(particlesForFit1, *Jpsi1, MuMuVertex);

      if ( !scFit1 ) {
        Warning("Fit error").ignore(); 
        continue; 
      } 
 
        const LHCb::Particle *Jpsi = this->markTree(Jpsi1);
        //const LHCb::Particle Jpsi2 = *Jpsi;
        JpsiVector.push_back(Jpsi);

        // JpsiVector.push_back(Jpsi1)
        sc = StatusCode::SUCCESS;
        
    }
  }
 
  if (!sc) return sc;
  counter("made jpsi")++;

  const LHCb::Particle::Range JpsiRange = LHCb::Particle::Range(JpsiVector.begin(),JpsiVector.end());
  //  if(m_NoAccCutMC )  sc = DaughtersInAcceptance(JpsiRange);

  sc = MakeJpsi(JpsiRange);
 
  return sc;
}
    

//=============================================================================
//Daughters in acceptance 
//============================================================================= 

StatusCode DoubleJpsiSel::DaughtersInAcceptance(const LHCb::Particle::Range &jpsi){
  StatusCode sc = StatusCode::SUCCESS;
  LHCb::Particle::ConstVector::const_iterator iJpsi;
  for (iJpsi = jpsi.begin(); iJpsi != jpsi.end(); ++iJpsi){
    const LHCb::Particle *Jpsi = (*iJpsi);
    counter("A. inDaugters In acceptance ")++;
    
    
    Tuple daughtersinacctuple = nTuple("daughersinacctuple");
   
    if ((Jpsi->particleID().pid() !=    443) &&
        (Jpsi->particleID().pid() != 100443)) continue;

    const LHCb::Particle *muPlus  = GetMuon(Jpsi, 1);
    const LHCb::Particle *muMinus = GetMuon(Jpsi, -1);
    const LHCb::ProtoParticle *muProtoPlus = muPlus->proto();
    const LHCb::Track         *muTrackPlus = muProtoPlus->track();
    const LHCb::ProtoParticle *muProtoMinus = muMinus->proto();
    const LHCb::Track         *muTrackMinus = muProtoMinus->track();

    counter("B. inDaugters In acceptance ")++;

    double motherPID = CheckMCTruth(muTrackPlus,muTrackMinus);
 
    if (motherPID == 0.0) continue;

    counter("C. inDaughters In Acceptance")++;
    
    Gaudi::LorentzVector jpsi4Mom = Jpsi->momentum();
    Gaudi::LorentzVector mu4Mom   = muPlus->momentum();
    ROOT::Math::Boost boost(jpsi4Mom.BoostToCM());
    const Gaudi::XYZVector boostedMu4Mom   = (boost(mu4Mom)).Vect().unit();
    const Gaudi::XYZVector boostedJpsi4Mom = jpsi4Mom.Vect().unit();
    double cosTheta = fabs(boostedMu4Mom.Dot(boostedJpsi4Mom));

    counter("B. before tuple filling inDaugters In acceptance ")++;

    daughtersinacctuple->fill("Jpsimass,Jpsipt,Rapidity,AbsCosTheta,truthMotherPID",
                     Jpsi->measuredMass(),Jpsi->pt(),Jpsi->momentum().Rapidity(),
                     cosTheta, motherPID);
      daughtersinacctuple->write();  
 
  }
  return sc;
}


//=============================================================
// Loop on J/psi
//=============================================================

StatusCode DoubleJpsiSel::MakeJpsi(const LHCb::Particle::Range &jpsi){
  StatusCode sc = StatusCode::SUCCESS;
  // Tuple efftuple = nTuple("efftuple");
  const IParticleCombiner *combiner = particleCombiner();
  std::vector<LHCb::Particle> jpsiVec;
  counter("A start MakeJpsi")++;

  std::vector<LHCb::Particle>          jpsivec;
  std::vector<const LHCb::VertexBase*> pvvec;
  std::vector<const LHCb::Particle*>   muonvec;

  StatEntity countpt,countpt2;

  const std::string prefixFourMuon = "FourMuons";
  const std::string prefixJpsi     = "Jpsi";
  const std::string prefixJpsi1    = "Jpsi_1";
  const std::string prefixJpsi2    = "Jpsi_2";
  const std::string prefixMother   = "DoubleJpsi";
  const std::string prefixMuPlus   = "MuPlus";
  const std::string prefixMuPlus1  = "MuPlus1";
  const std::string prefixMuPlus2  = "MuPlus2";
  const std::string prefixMuMinus  = "MuMinus";
  const std::string prefixMuMinus1 = "MuMinus1";
  const std::string prefixMuMinus2 = "MuMinus2";

  //
  // Non configurable cut values.
  //

  double maxDTFChi2Ndof     = 5.0;


  /*
   * Non configurable cut values.
   */
 
  double maxFourMuonVtxChi2 = 100;
 
  /*
   * Main loop over all J/Psi to select them.
   */
 
  const LHCb::ParticleID psi2SPID(m_pidPsi2S);
  const LHCb::ParticleID chiB2PID(m_pidChiB21P);
  

  //
  // Main loop over all J/Psi to select them.
  //
  LHCb::Particle::ConstVector::const_iterator iJpsi;
  for (iJpsi = jpsi.begin(); iJpsi != jpsi.end(); ++iJpsi){
    const LHCb::Particle *Jpsi = (*iJpsi);
    
    counter("Az. start of Jpsi loop")++;

    if ((Jpsi->particleID().pid() !=    443) &&
        (Jpsi->particleID().pid() != 100443)) continue;

    counter("B. after pid check")++;
    
    
    const LHCb::Particle *muPlus  = GetMuon(Jpsi, 1);
    const LHCb::Particle *muMinus = GetMuon(Jpsi, -1);
    
    counter("C. after muons")++;
    
    //  
    // Apply the cuts on the muons and the J/Psi
    //
    
    const LHCb::ProtoParticle *muProtoPlus = muPlus->proto();
    const LHCb::Track         *muTrackPlus = muProtoPlus->track();

    const LHCb::ProtoParticle *muProtoMinus = muMinus->proto();
    const LHCb::Track         *muTrackMinus = muProtoMinus->track();

    //==========================================
    
    //
    // Make a fake particle decay to create the vertex for the  J/Psi.
    //
  
    //
    //BELOW and cuts further below JUST FOR TEST into MCtruth Check - make sure not on normally
    //
    bool testnumbMC = true; //BELOW and cuts further below JUST FOR TEST into MCtruth Check - make sure not on normally
    if(false){

      IBackgroundCategory::categories cat = IBackgroundCategory::Undefined;
      std::vector<IBackgroundCategory*>::const_iterator it;
      for (it = m_bkgs.begin(); it != m_bkgs.end(); ++it){
        cat = (*it)->category(Jpsi);
        int check = 1;
        if (cat != IBackgroundCategory::Undefined) break;
        std::cout <<"iterator: " << check << " bkgd cat " << cat << std::endl;
        check++;
      }
      
   
    }
    

    LHCb::Vertex   motherVertex;
    LHCb::Particle mother(LHCb::ParticleID(310)); 
    
    LHCb::Particle::ConstVector particlesForFit;
    particlesForFit.push_back(muMinus);
    particlesForFit.push_back(muPlus);
    
    StatusCode scFit = combiner->combine(particlesForFit, mother, motherVertex);

    if (!scFit){
      Warning("Double J/Psi fit error").ignore(); 
      continue;
    }
      
    //
    // Run a DecayTreeFitter on the mother, first checking there is a PV
    //
    
    const LHCb::VertexBase *pv = bestPV(&mother);
    if(!pv)continue;  //put in now 
    counter("Cz. after (!pv) cut")++;
    DecayTreeFitter::Fitter dtFitter(mother, pv);

    //     dtFitter.setMassConstraint(Jpsi1->particleID());

    dtFitter.fit();

    LHCb::DecayTree  motherTree = dtFitter.getFittedTree();
    LHCb::Particle  *fitMother  = motherTree.head();

    Gaudi::XYZVector dl = pv->position() - fitMother->endVertex()->position();
   
    bool acceptDoubleJpsi = true;
    if (fitMother->endVertex()->chi2PerDoF() > maxDTFChi2Ndof) acceptDoubleJpsi = false;

    counter("D. after DTF calculation but not cut")++;

    //
    // Make sure all muons are different.
    //  
   

    bool sameMuon  = false;
    bool sameProto = false;
    bool sameTrack = false;

    
    if (muPlus->key()  == muMinus->key())  sameMuon = true;
    
    if (sameMuon)  continue;
    
    if (muProtoMinus->key()  == muProtoPlus->key())  sameProto = true;
    
    if (sameProto) continue;
    
    if (muTrackPlus->key()  == muTrackMinus->key())  sameTrack = true;
    
    if (sameTrack) continue;
    
    counter("E. after same muon ")++;

    double motherPID = CheckMCTruth(muTrackPlus,muTrackMinus);
    counter("F. after running MCTruth,but not cutting, just before tuple")++;

    Gaudi::LorentzVector jpsi4Mom = Jpsi->momentum();
    Gaudi::LorentzVector mu4Mom   = muPlus->momentum();
    ROOT::Math::Boost boost(jpsi4Mom.BoostToCM());
    const Gaudi::XYZVector boostedMu4Mom   = (boost(mu4Mom)).Vect().unit();
    const Gaudi::XYZVector boostedJpsi4Mom = jpsi4Mom.Vect().unit();
    double cosTheta = fabs(boostedMu4Mom.Dot(boostedJpsi4Mom));
    double Theta = acos(boostedMu4Mom.Dot(boostedJpsi4Mom));
    

    Gaudi::LorentzVector p1 = muPlus->momentum();
    Gaudi::LorentzVector p2 = muMinus->momentum();

    double pp1      = sqrt(p1.X()*p1.X() + p1.Y()*p1.Y());
    double pp2      = sqrt(p2.X()*p2.X() + p2.Y()*p2.Y());
    double cosPhi   = (p1.X()*p2.X() + p1.Y()*p2.Y())/(pp1*pp2);
    double deltaPhi = acos(cosPhi);

    // FillDataPid(allJpsi,muPlus,prefixMuPlus);

    double PlusProbNNMu  = muProtoPlus->info(LHCb::ProtoParticle::ProbNNmu,    -1000);
    double MinusProbNNMu = muProtoMinus->info(LHCb::ProtoParticle::ProbNNmu,    -1000);

    if(m_fillAllJpsi){
      Tuple allJpsi = nTuple("AllJpsi");
      FillDataEventNtuple(allJpsi);
      FillEventNtuple(allJpsi);
      FillJpsiNtuple(allJpsi, Jpsi,   prefixJpsi);
      FillMuNtuple(allJpsi,   muPlus,  prefixMuPlus);
      FillMuNtuple(allJpsi,   muMinus, prefixMuMinus);  
      FillPolarizationNtuple(allJpsi, Jpsi, muPlus, prefixJpsi);
      FillDataTrackGhostProb(allJpsi, muTrackPlus, prefixMuPlus);
      FillDataTrackGhostProb(allJpsi, muTrackMinus, prefixMuMinus);
      FillDataDTF(allJpsi, fitMother->endVertex(), "Jpsi");
      FillDataMCTruth(allJpsi, Jpsi, "Jpsi"); //this fills MCAncestor
      FillDataPid(allJpsi,muPlus,prefixMuPlus);
      FillDataPid(allJpsi,muMinus,prefixMuMinus);
      FillDataImpactParameter(allJpsi, muPlus, prefixMuPlus);
      FillDataImpactParameter(allJpsi, muMinus, prefixMuMinus);
      FillDataPropertime(allJpsi, Jpsi, prefixJpsi);
      //  FillMCAncestor(allJpsi, Jpsi, prefixJpsi);
      allJpsi->column("Jpsi_Theta",double(Theta));
      allJpsi->column("Jpsi_Phi",double(deltaPhi));
      allJpsi->write();
    }
    

  
    if (Jpsi->measuredMass() > m_maxJpsiMass)        continue;
    if (Jpsi->measuredMass() < m_minJpsiMass)        continue;
    counter("G.after mass cut")++;
    if (Jpsi->endVertex()->chi2() > m_maxJpsiChi2)   continue;
    counter("H.after vertex Chi2")++;
    //cut on DLL - turn OFF for acc, rec, sel efficiency
    // if (muProtoPlus->info(601, -1)  <= m_minMuDll)  continue; //need to turn off except for test
    // if (muProtoMinus->info(601, -1) <= m_minMuDll)  continue; //need to turn off except for test
    counter("I.after DLL")++; //Need to turn off except for test!
    if (muTrackPlus->pt()  < m_minMuPt)               continue;
    if (muTrackMinus->pt() < m_minMuPt)               continue;
    counter("J.after pt cut")++;
    if (muTrackPlus->info(101, -1)  != -1)            continue;
    if (muTrackMinus->info(101, -1) != -1)            continue;
    counter("K.after kullback")++;
    if (muTrackPlus->chi2PerDoF()  > m_maxMuChi2NDof) continue;
    if (muTrackMinus->chi2PerDoF() > m_maxMuChi2NDof) continue;
    counter("L.after Track Chi2")++;

   
    //  const Gaudi::Math::ParticleParams *fitParams  = dtFitter.fitParams(&mother);
    
    if (muTrackPlus->ghostProbability()>0.3) continue;
    counter("M. after MuPlus Ghost")++;
    if (muTrackMinus->ghostProbability()>0.3)continue;
    counter("N. after MuMinus Ghost")++;
    if (muTrackPlus->pseudoRapidity() <2.0) continue ;
    counter("O. after MuTrackPlus PseudoRap <2.0")++;
    if (muTrackPlus->pseudoRapidity() >4.5) continue ;
    counter("P. after MuTrackPlus PseudoRap > 4.5")++;
    if (muTrackMinus->pseudoRapidity()<2.0) continue ;
    counter("Q. after MuTrackMinus pseudoRap <2.0")++;
    if (muTrackMinus->pseudoRapidity()>4.5) continue ;
    counter("R. after MuTrackMinus pseudoRap > 4.5")++;
    if(Jpsi->pt() > 10000.0) continue;
    counter("Rx after Jpsi pt cut")++;
    if(Jpsi->momentum().Rapidity() <2.0)continue;
    if(Jpsi->momentum().Rapidity() > 4.5)continue;
    counter("Ry after rapidity cut")++;
    if(Jpsi->p()< m_minMomentum )continue;

    //
    // Apply the cut on the DecayTreeFitter Chi2/Ndof
    //

    if (!acceptDoubleJpsi) continue; 

    //
    //Apply cut on MC truth
    //
  
    if (motherPID == 0.0) continue;

    counter("S.After TruthCheck")++;
    
    counter("T.FINAL COUNT")++;
  }


  return sc;  
}



//=============================================================
// Loop on Double J/psi for Double Jpsi MC
//=============================================================

StatusCode DoubleJpsiSel::MakeDoubleJpsi(const LHCb::Particle::Range &jpsi){
  StatusCode sc = StatusCode::SUCCESS;
  Tuple efftuple = nTuple("efftuple");
  const IParticleCombiner *combiner = particleCombiner();
  std::vector<LHCb::Particle> jpsiVec;
  counter("A start MakeJpsi")++;
  
  std::vector<LHCb::Particle>          jpsivec;
  std::vector<const LHCb::VertexBase*> pvvec;
  std::vector<const LHCb::Particle*>   muonvec;

  StatEntity countpt,countpt2;

  const std::string prefixFourMuon = "FourMuons";
  const std::string prefixJpsi     = "JPsi";
  const std::string prefixJpsi1    = "JPsi_1";
  const std::string prefixJpsi2    = "JPsi_2";
  const std::string prefixMother   = "DoubleJpsi";
  const std::string prefixMuPlus   = "MuPlus";
  const std::string prefixMuPlus1  = "MuPlus1";
  const std::string prefixMuPlus2  = "MuPlus2";
  const std::string prefixMuMinus  = "MuMinus";
  const std::string prefixMuMinus1 = "MuMinus1";
  const std::string prefixMuMinus2 = "MuMinus2";

  //
  // Non configurable cut values.
  //

  double maxDTFChi2Ndof     = 5.0;
  double maxFourMuonVtxChi2 = 100;
 
  const LHCb::ParticleID psi2SPID(m_pidPsi2S);
  const LHCb::ParticleID chiB2PID(m_pidChiB21P);

  //
  // Main loop over all J/Psi to select them.
  //
  LHCb::Particle::ConstVector::const_iterator iJpsi;
  for (iJpsi = jpsi.begin(); iJpsi != jpsi.end(); ++iJpsi){
    const LHCb::Particle *Jpsi = (*iJpsi);
    
    counter("Az. start of Jpsi loop")++;

    if ((Jpsi->particleID().pid() !=    443) &&
        (Jpsi->particleID().pid() != 100443)) continue;

    counter("B. after pid check")++;
    
    //
    //Get Muons - not sure why, but GetMuons seems to assign the 'wrong way round' - so hack,below
    //is just to change sign of call
    //
    const LHCb::Particle *muPlus  = GetMuon(Jpsi, -1); 
    const LHCb::Particle *muMinus = GetMuon(Jpsi, 1);
    //std::cout << "muplus charge " << muPlus->charge() << " pid " << muPlus->particleID() << std::endl;
    // std::cout << "muminus charge " << muMinus->charge() << " pid " << muMinus->particleID() << std::endl;
    
    counter("C. after muons")++;
    
    //  
    // Apply the cuts on the muons and the J/Psi
    //
    
    const LHCb::ProtoParticle *muProtoPlus = muPlus->proto();
    const LHCb::Track         *muTrackPlus = muProtoPlus->track();

    const LHCb::ProtoParticle *muProtoMinus = muMinus->proto();
    const LHCb::Track         *muTrackMinus = muProtoMinus->track();

    //
    // Make a fake particle decay to create the vertex for the  J/Psi.
    //

    LHCb::Vertex   motherVertex;
    LHCb::Particle mother(LHCb::ParticleID(310));

    LHCb::Particle::ConstVector particlesForFit;
    particlesForFit.push_back(muMinus);
    particlesForFit.push_back(muPlus);
      
    StatusCode scFit = combiner->combine(particlesForFit, mother, motherVertex);

    if (!scFit){
      Warning("Double J/Psi fit error").ignore(); 
      continue;
    }    
    
    //
    // Run a DecayTreeFitter on the mother, first checking there is a PV
    //  
    const LHCb::VertexBase *pv = bestPV(&mother);
    if(!pv)continue;  //put in now 
    counter("Cz. after (!pv) cut")++;
    DecayTreeFitter::Fitter dtFitter(mother, pv);

    dtFitter.fit();

    LHCb::DecayTree  motherTree = dtFitter.getFittedTree();
    LHCb::Particle  *fitMother  = motherTree.head();

    Gaudi::XYZVector dl = pv->position() - fitMother->endVertex()->position();
   
    bool acceptDoubleJpsi = true;
    if (fitMother->endVertex()->chi2PerDoF() > maxDTFChi2Ndof) acceptDoubleJpsi = false;

    counter("D. after DTF calculation but not cut")++;
   
    /*
    //-----------second method
    const LHCb::Particle Jpsi2 = *Jpsi;
    const LHCb::VertexBase *pv = bestPV(&Jpsi2);
       
    DecayTreeFitter::Fitter dtFitter(Jpsi2, pv);
    dtFitter.fit();
    
    LHCb::DecayTree  motherTree = dtFitter.getFittedTree();
    LHCb::Particle  *fitMother  = motherTree.head();
    //FillDataDTF(selJpsi, fitMother->endVertex(), "JPsi");

    counter("Ab.before DTF")++;
    if (fitMother->endVertex()->chi2PerDoF() > maxDTFChi2Ndof) continue;
    counter("Ac.After DTF")++;
    //=========================================
    */

    //
    // Make sure all muons are different.
    //  
   

    bool sameMuon  = false;
    bool sameProto = false;
    bool sameTrack = false;

    
    if (muPlus->key()  == muMinus->key())  sameMuon = true;
    
    if (sameMuon)  continue;
    
    if (muProtoMinus->key()  == muProtoPlus->key())  sameProto = true;
    
    if (sameProto) continue;
    
    if (muTrackPlus->key()  == muTrackMinus->key())  sameTrack = true;
    
    if (sameTrack) continue;
    
    counter("E. after same muon ")++;

    double motherPID = CheckMCTruth(muTrackPlus,muTrackMinus);
     
    counter("F. after running MCTruth,but not cutting, just before tuple")++;

    Gaudi::LorentzVector jpsi4Mom = Jpsi->momentum();
    Gaudi::LorentzVector mu4Mom   = muPlus->momentum();
    ROOT::Math::Boost boost(jpsi4Mom.BoostToCM());
    const Gaudi::XYZVector boostedMu4Mom   = (boost(mu4Mom)).Vect().unit();
    const Gaudi::XYZVector boostedJpsi4Mom = jpsi4Mom.Vect().unit();
    double cosTheta = fabs(boostedMu4Mom.Dot(boostedJpsi4Mom));
    double Theta = acos(boostedMu4Mom.Dot(boostedJpsi4Mom));
    

    Gaudi::LorentzVector p1 = muPlus->momentum();
    Gaudi::LorentzVector p2 = muMinus->momentum();

    double pp1      = sqrt(p1.X()*p1.X() + p1.Y()*p1.Y());
    double pp2      = sqrt(p2.X()*p2.X() + p2.Y()*p2.Y());
    double cosPhi   = (p1.X()*p2.X() + p1.Y()*p2.Y())/(pp1*pp2);
    double deltaPhi = acos(cosPhi);

    // FillDataPid(allJpsi,muPlus,prefixMuPlus);

    double PlusProbNNMu  = muProtoPlus->info(LHCb::ProtoParticle::ProbNNmu,    -1000);
    double MinusProbNNMu = muProtoMinus->info(LHCb::ProtoParticle::ProbNNmu,    -1000);

    /*
     * Fill an ntuple with the data for all Jpsi and muons
     * before the cuts are applied. 
     */

    if(m_fillAllJpsi){
      Tuple allJpsi = nTuple("AllJpsi");
      FillDataEventNtuple(allJpsi);
      FillEventNtuple(allJpsi);
      FillJpsiNtuple(allJpsi, Jpsi,   prefixJpsi);
      FillMuNtuple(allJpsi,   muPlus,  prefixMuPlus);
      FillMuNtuple(allJpsi,   muMinus, prefixMuMinus);  
      FillPolarizationNtuple(allJpsi, Jpsi, muPlus, prefixJpsi);
      FillDataTrackGhostProb(allJpsi, muTrackPlus, prefixMuPlus);
      FillDataTrackGhostProb(allJpsi, muTrackMinus, prefixMuMinus);
      FillDataDTF(allJpsi, fitMother->endVertex(), "JPsi");
      //  FillDataMCTruth(allJpsi, motherPID, "JPsi");
      FillDataMCTruth(allJpsi, Jpsi, "JPsi");
      FillDataPid(allJpsi,muPlus,prefixMuPlus);
      FillDataPid(allJpsi,muMinus,prefixMuMinus);
      FillDataImpactParameter(allJpsi, muPlus, prefixMuPlus);
      FillDataImpactParameter(allJpsi, muMinus, prefixMuMinus);
      FillDataPropertime(allJpsi, Jpsi, prefixJpsi);
      allJpsi->column("Jpsi_Theta",double(Theta));
      allJpsi->column("Jpsi_Phi",double(deltaPhi));
      allJpsi->write();
    }
    
    //
    //Turn off cuts, so all data can be filled into tuple and manipulated later.
    //
    /*
    if (Jpsi->measuredMass() > m_maxJpsiMass)        continue;
    if (Jpsi->measuredMass() < m_minJpsiMass)        continue;
    counter("G.after mass cut")++;
    if (Jpsi->endVertex()->chi2() > m_maxJpsiChi2)   continue;
    counter("H.after vertex Chi2")++;
    //cut on DLL - turn OFF for acc, rec, sel efficiency
    // if (muProtoPlus->info(601, -1)  <= m_minMuDll)  continue; //need to turn off except for test!
    // if (muProtoMinus->info(601, -1) <= m_minMuDll)  continue; //need to turn off except for test!
    counter("I.after DLL")++; //Need to turn off except for test!
    if (muTrackPlus->pt()  < m_minMuPt)               continue;
    if (muTrackMinus->pt() < m_minMuPt)               continue;
    counter("J.after pt cut")++;
    if (muTrackPlus->info(101, -1)  != -1)            continue;
    if (muTrackMinus->info(101, -1) != -1)            continue;
    counter("K.after kullback")++;
    if (muTrackPlus->chi2PerDoF()  > m_maxMuChi2NDof) continue;
    if (muTrackMinus->chi2PerDoF() > m_maxMuChi2NDof) continue;
    counter("L.after Track Chi2")++;
    if (muTrackPlus->ghostProbability()>0.3) continue;
    counter("M. after MuPlus Ghost")++;
    if (muTrackMinus->ghostProbability()>0.3)continue;
    counter("N. after MuMinus Ghost")++;
    if (muTrackPlus->pseudoRapidity() <2.0) continue ;
    counter("O. after MuTrackPlus PseudoRap <2.0")++;
    if (muTrackPlus->pseudoRapidity() >4.5) continue ;
    counter("P. after MuTrackPlus PseudoRap > 4.5")++;
    if (muTrackMinus->pseudoRapidity()<2.0) continue ;
    counter("Q. after MuTrackMinus pseudoRap <2.0")++;
    if (muTrackMinus->pseudoRapidity()>4.5) continue ;
    counter("R. after MuTrackMinus pseudoRap > 4.5")++;
    if(Jpsi->pt() > 10000.0) continue;
    counter("Rx after Jpsi pt cut")++;
    if(Jpsi->momentum().Rapidity() <2.0)continue;
    if(Jpsi->momentum().Rapidity() > 4.5)continue;
    counter("Ry after rapidity cut")++;
    if(Jpsi->p()< m_minMomentum )continue;
   
    //
    // Apply the cut on the DecayTreeFitter Chi2/Ndof
    //

    if (!acceptDoubleJpsi) continue;

    //
    //Apply cut on MC truth
    //  
    if (motherPID == 0.0) continue;

    */

    counter("T.FINAL single COUNT")++;
    
    const LHCb::Particle jpsi = *Jpsi;
    jpsiVec.push_back(jpsi);
  }


  /*
   * Bail out if less than 2 J/psi.
   */
 
  int size = jpsiVec.size();
  //  std::cout << "number of Jpsi"  << jpsiVec.size()   << std::endl;
  
  if (size < 2) return sc;
  counter("MORE THAN ONE Jpsi!")++;

  for (int j = 1; j <= (size-1) ; j=(j+1)){
    for (int k = (j+1); k <= size; k=(k+1) ){        
      bool acceptDoubleJpsi = true;
 
      const LHCb::Particle *jpsi1 = NULL;
      const LHCb::Particle *jpsi2 = NULL;
      counter("Z 1")++;

      if ((jpsiVec[j-1]).pt() > (jpsiVec[k-1]).pt()){
        jpsi1 = (&jpsiVec[j-1]);
        jpsi2 = (&jpsiVec[k-1]);
      } else{
        jpsi1 = (&jpsiVec[k-1]);
        jpsi2 = (&jpsiVec[j-1]);
      }
      
      /*
       * Make sure all muons are different.
       */
      counter("Z 2")++;
      const LHCb::Particle *muPlus1  = GetParticles(jpsi1,  1,  m_pidMuon);
      if (muPlus1  == NULL) continue;
      counter("Z 2a")++;
      const LHCb::Particle *muMinus1 = GetParticles(jpsi1, -1, -m_pidMuon);
      if (muMinus1 == NULL) continue;
      counter("Z 2b")++;
      const LHCb::Particle *muPlus2  = GetParticles(jpsi2,  1,  m_pidMuon);
      if (muPlus2  == NULL) continue;
      counter("Z 2c")++;
      const LHCb::Particle *muMinus2 = GetParticles(jpsi2, -1, -m_pidMuon);
      if (muMinus2 == NULL) continue;
      counter("Z 3")++;
      StatusCode sc1 = DoubleJpsiCheckMuons(muPlus1, muMinus1, muPlus2, muMinus2);
      if (sc1.isFailure()){
        counter("Shared muons")++;
        continue;
      } 
      counter("Z 4")++;
      /*
       * Make a vertex with the 4 muons.
       */
 
      LHCb::Particle::ConstVector fourMuonsForFit;
      fourMuonsForFit.push_back(muPlus1);
      fourMuonsForFit.push_back(muMinus1);
      fourMuonsForFit.push_back(muPlus2);
      fourMuonsForFit.push_back(muMinus2);
 
      LHCb::Vertex   fourMuonVertex;
      LHCb::Particle fourMuon(chiB2PID);
 
      const IParticleCombiner *combiner = particleCombiner();
      StatusCode scFit = combiner->combine(fourMuonsForFit, fourMuon,
                                           fourMuonVertex);
      if (!scFit){
        warning() << "Four muons fit errror. jpsi1 id = "
                  << jpsi1->particleID().pid()
                  << " jpsi2 id = "
                  << jpsi2->particleID().pid()
                  << endmsg;
        counter("Four muons fit error")++;
        continue;
      }
 
      if (fourMuonVertex.chi2() > maxFourMuonVtxChi2) acceptDoubleJpsi = false;
     
      /*
       * Make a fake particle decay to create the vertex for the 2 J/Psi.
        */
      counter("Z 5")++;
      LHCb::Vertex   motherVertex;
      LHCb::Particle mother(chiB2PID);
 
      LHCb::Particle::ConstVector particlesForFit;
 
      particlesForFit.push_back(jpsi1);
      particlesForFit.push_back(jpsi2);
 
      scFit = combiner->combine(particlesForFit, mother, motherVertex);
 
      if (!scFit){
        Warning("Double J/psi fit error").ignore();
        counter("Double J/psi fit error")++;
        continue;
      }
     
      /*
       * Run a DecayTreeFitter on the mother
       */
     
      const LHCb::VertexBase *pv = bestPV(&mother);
      if (!pv) continue;
      counter("Z 6")++;
       const ITrackStateProvider *stateProvider = m_stateProvider.empty() ?
        NULL : &(*m_stateProvider);
 
      DecayTreeFitter::Fitter dtFitter(mother, pv, stateProvider);
      //      DecayTreeFitter::Fitter dtFitter(mother, pv);
      dtFitter.fit();
     
      LHCb::DecayTree  motherTree = dtFitter.getFittedTree();
      LHCb::Particle  *fitMother  = motherTree.head();
     
      //
      //Do DTF info for individual Jpsi vertexes
      //
      
      //
      // Make a fake particle decay to create the vertex for the  J/Psi.
      //
      
      //
      //Jpsi1
      //
      
      LHCb::Vertex   motherVertex1;
      LHCb::Particle mother1(LHCb::ParticleID(310));
      
      LHCb::Particle::ConstVector particlesForFit1;
      particlesForFit1.push_back(muMinus1);
      particlesForFit1.push_back(muPlus1);
      
      StatusCode scFit1 = combiner->combine(particlesForFit1, mother1, motherVertex1);

      if (!scFit1){
        Warning("Double J/Psi fit error").ignore(); 
        continue;
      }    
      
      //
      // Run a DecayTreeFitter on the mother, first checking there is a PV
      //  
      const LHCb::VertexBase *pv1 = bestPV(&mother1);
      if(!pv1)continue;  //put in now 
      
      DecayTreeFitter::Fitter dtFitter1(mother1, pv1);
      
      dtFitter1.fit();

      LHCb::DecayTree  motherTree1 = dtFitter1.getFittedTree();
      const LHCb::Particle  *fitMother1  = motherTree1.head();
      
      Gaudi::XYZVector dl1 = pv1->position() - fitMother1->endVertex()->position();
      
      
      //
      //Jpsi2
      //
      LHCb::Vertex   motherVertex2;
      LHCb::Particle mother2(LHCb::ParticleID(310));
      
      LHCb::Particle::ConstVector particlesForFit2;
      particlesForFit2.push_back(muMinus2);
      particlesForFit2.push_back(muPlus2);
      
      StatusCode scFit2 = combiner->combine(particlesForFit2, mother2, motherVertex2);

      if (!scFit2){
        Warning("Double J/Psi fit error").ignore(); 
        continue;
      }    
      
      //
      // Run a DecayTreeFitter on the mother, first checking there is a PV
      //  
      const LHCb::VertexBase *pv2 = bestPV(&mother2);
      if(!pv2)continue;  //put in now 
      
      DecayTreeFitter::Fitter dtFitter2(mother2, pv2);
      
      dtFitter2.fit();

      LHCb::DecayTree  motherTree2 = dtFitter2.getFittedTree();
      const LHCb::Particle  *fitMother2  = motherTree2.head();
      
      Gaudi::XYZVector dl2 = pv2->position() - fitMother2->endVertex()->position();
      //
      //Theta and phi info            
      //
 
      //
      //Jpsi1
      //

      Gaudi::LorentzVector jpsi4Mom1 = jpsi1->momentum();
      Gaudi::LorentzVector mu4Mom1   = muPlus1->momentum();
      ROOT::Math::Boost boost1(jpsi4Mom1.BoostToCM());
      const Gaudi::XYZVector boostedMu4Mom1   = (boost1(mu4Mom1)).Vect().unit();
      const Gaudi::XYZVector boostedJpsi4Mom1 = jpsi4Mom1.Vect().unit();
   
      double cosTheta1 = fabs(boostedMu4Mom1.Dot(boostedJpsi4Mom1));
      double Theta1 = acos(boostedMu4Mom1.Dot(boostedJpsi4Mom1));

      Gaudi::LorentzVector p11 = muPlus1->momentum();
      Gaudi::LorentzVector p21 = muMinus1->momentum();

      double pp11      = sqrt(p11.X()*p11.X() + p11.Y()*p11.Y());
      double pp21      = sqrt(p21.X()*p21.X() + p21.Y()*p21.Y());
      double cosPhi1   = (p11.X()*p21.X() + p11.Y()*p21.Y())/(pp11*pp21);
      double deltaPhi1 = acos(cosPhi1);

      //
      //Jpsi2
      //

      Gaudi::LorentzVector jpsi4Mom2 = jpsi2->momentum();
      Gaudi::LorentzVector mu4Mom2   = muPlus2->momentum();
      ROOT::Math::Boost boost2(jpsi4Mom2.BoostToCM());
      const Gaudi::XYZVector boostedMu4Mom2   = (boost2(mu4Mom2)).Vect().unit();
      const Gaudi::XYZVector boostedJpsi4Mom2 = jpsi4Mom2.Vect().unit();
     
      double cosTheta2 = fabs(boostedMu4Mom2.Dot(boostedJpsi4Mom2));
      double Theta2 = acos(boostedMu4Mom2.Dot(boostedJpsi4Mom2));

      Gaudi::LorentzVector p12 = muPlus2->momentum();
      Gaudi::LorentzVector p22 = muMinus2->momentum();

      double pp12      = sqrt(p12.X()*p12.X() + p12.Y()*p12.Y());
      double pp22      = sqrt(p22.X()*p22.X() + p22.Y()*p22.Y());
      double cosPhi2   = (p12.X()*p22.X() + p12.Y()*p22.Y())/(pp12*pp22);
      double deltaPhi2 = acos(cosPhi2);

      //
      //truth info 
      //
      const LHCb::ProtoParticle *muProtoPlus1 = muPlus1->proto();
      const LHCb::Track         *muTrackPlus1 = muProtoPlus1->track();
      const LHCb::ProtoParticle *muProtoMinus1 = muMinus1->proto();
      const LHCb::Track         *muTrackMinus1 = muProtoMinus1->track();
      
      const LHCb::ProtoParticle *muProtoPlus2 = muPlus2->proto();
      const LHCb::Track         *muTrackPlus2 = muProtoPlus2->track();
      const LHCb::ProtoParticle *muProtoMinus2 = muMinus2->proto();
      const LHCb::Track         *muTrackMinus2 = muProtoMinus2->track();
      

      double motherPID1 = CheckMCTruth(muTrackPlus1,muTrackMinus1);
      double motherPID2 = CheckMCTruth(muTrackPlus2,muTrackMinus2);
      
      /*
       * Fill a ntuple with the data of all possible double J/Psi.
       */
      counter("All double Jpsi")++;
      //if (m_fillAllDoubleJpsi){ // AC
      if (true){
        Tuple allDoubleJpsi = nTuple("AllDoubleJpsi");
 
        FillAnalysisDoubleJpsiMuMu(allDoubleJpsi, jpsi1, muPlus1, muMinus1,
                                   jpsi2, muPlus2, muMinus2, &fourMuon, fitMother,fitMother1,fitMother2);
        allDoubleJpsi->column("Jpsi_1_Theta",double(Theta1));
        allDoubleJpsi->column("Jpsi_1_Phi",double(deltaPhi1));
        allDoubleJpsi->column("Jpsi_2_Theta",double(Theta2));
        allDoubleJpsi->column("Jpsi_2_Phi",double(deltaPhi2));
        FillDataMCTruth(allDoubleJpsi, jpsi1, "JPsi_1");
        FillDataMCTruth(allDoubleJpsi, jpsi2, "JPsi_2");
        //  FillDataMCTruth(allDoubleJpsi, motherPID1, "JPsi_1");
        //FillDataMCTruth(allDoubleJpsi, motherPID2, "JPsi_2");
         allDoubleJpsi->write();
      }

    }
    
  }

  return sc;  
}




StatusCode DoubleJpsiSel::FromBDTFTest(const LHCb::Particle::Range &jpsi){
  StatusCode sc = StatusCode::SUCCESS;
  Tuple efftuple = nTuple("efftuple");
  const IParticleCombiner *combiner = particleCombiner();
  std::vector<LHCb::Particle> jpsiVec;
  counter("A start MakeJpsi")++;
  
  std::vector<LHCb::Particle>          jpsivec;
  std::vector<const LHCb::VertexBase*> pvvec;
  std::vector<const LHCb::Particle*>   muonvec;

  StatEntity countpt,countpt2;

  const std::string prefixFourMuon = "FourMuons";
  const std::string prefixJpsi     = "JPsi";
  const std::string prefixJpsi1    = "JPsi_1";
  const std::string prefixJpsi2    = "JPsi_2";
  const std::string prefixMother   = "DoubleJpsi";
  const std::string prefixMuPlus   = "MuPlus";
  const std::string prefixMuPlus1  = "MuPlus1";
  const std::string prefixMuPlus2  = "MuPlus2";
  const std::string prefixMuMinus  = "MuMinus";
  const std::string prefixMuMinus1 = "MuMinus1";
  const std::string prefixMuMinus2 = "MuMinus2";

  //
  // Non configurable cut values.
  //

  double maxDTFChi2Ndof     = 5.0;


  /*
   * Non configurable cut values.
   */
 
  double maxFourMuonVtxChi2 = 100;
 
  /*
   * Main loop over all J/Psi to select them.
   */
 
  const LHCb::ParticleID psi2SPID(m_pidPsi2S);
  const LHCb::ParticleID chiB2PID(m_pidChiB21P);
  

  //
  // Main loop over all J/Psi to select them.
  //
  LHCb::Particle::ConstVector::const_iterator iJpsi;
  for (iJpsi = jpsi.begin(); iJpsi != jpsi.end(); ++iJpsi){
    const LHCb::Particle *Jpsi = (*iJpsi);
    
    counter("Az. start of Jpsi loop")++;

    if ((Jpsi->particleID().pid() !=    443) &&
        (Jpsi->particleID().pid() != 100443)) continue;

    counter("B. after pid check")++;
    
    
    const LHCb::Particle *muPlus  = GetMuon(Jpsi, 1);
    const LHCb::Particle *muMinus = GetMuon(Jpsi, -1);
    
    counter("C. after muons")++;
    
    //  
    // Apply the cuts on the muons and the J/Psi
    //
    
    const LHCb::ProtoParticle *muProtoPlus = muPlus->proto();
    const LHCb::Track         *muTrackPlus = muProtoPlus->track();

    const LHCb::ProtoParticle *muProtoMinus = muMinus->proto();
    const LHCb::Track         *muTrackMinus = muProtoMinus->track();

    
    //
    // Make a fake particle decay to create the vertex for the  J/Psi.
    //  
    
    bool testnumbMC = true; //BELOW and cuts further below JUST FOR TEST into MCtruth Check - make sure not on normally
    if(false){

      IBackgroundCategory::categories cat = IBackgroundCategory::Undefined;
      std::vector<IBackgroundCategory*>::const_iterator it;
      for (it = m_bkgs.begin(); it != m_bkgs.end(); ++it){
        cat = (*it)->category(Jpsi);
        int check = 1;
        if (cat != IBackgroundCategory::Undefined) break;
        std::cout <<"iterator: " << check << " bkgd cat " << cat << std::endl;
        check++;
      }
      
    }
    

    LHCb::Vertex   motherVertex;
    LHCb::Particle mother(LHCb::ParticleID(310)); 
    

    LHCb::Particle::ConstVector particlesForFit;
    particlesForFit.push_back(muMinus);
    particlesForFit.push_back(muPlus);

      
    StatusCode scFit = combiner->combine(particlesForFit, mother, motherVertex);

    if (!scFit){
      Warning("Double J/Psi fit error").ignore(); 
      continue;
    }
      
    
    //
    // Run a DecayTreeFitter on the mother
    //
    
    const LHCb::VertexBase *pv = bestPV(&mother);
    if(!pv)continue;  //put in now 
    counter("Cz. after (!pv) cut")++;
    DecayTreeFitter::Fitter dtFitter(mother, pv);

    //     dtFitter.setMassConstraint(Jpsi1->particleID());

    dtFitter.fit();

    LHCb::DecayTree  motherTree = dtFitter.getFittedTree();
    LHCb::Particle  *fitMother  = motherTree.head();

    Gaudi::XYZVector dl = pv->position() - fitMother->endVertex()->position();
   
    bool acceptDoubleJpsi = true;
    if (fitMother->endVertex()->chi2PerDoF() > maxDTFChi2Ndof) acceptDoubleJpsi = false;

    counter("D. after DTF calculation but not cut")++;
    
  
    //
    // Make sure all muons are different.
    //  
   

    bool sameMuon  = false;
    bool sameProto = false;
    bool sameTrack = false;

    
    if (muPlus->key()  == muMinus->key())  sameMuon = true;
    
    if (sameMuon)  continue;
    
    if (muProtoMinus->key()  == muProtoPlus->key())  sameProto = true;
    
    if (sameProto) continue;
    
    if (muTrackPlus->key()  == muTrackMinus->key())  sameTrack = true;
    
    if (sameTrack) continue;
    
    counter("E. after same muon ")++;

    double motherPID = CheckMCTruth(muTrackPlus,muTrackMinus);
    counter("F. after running MCTruth,but not cutting, just before tuple")++;


    //
    // Apply the cut on the DecayTreeFitter Chi2/Ndof
    //

    //////counter("D.before DTF")++;
    //if (!acceptDoubleJpsi) continue; TURNED OFF FOR TEST

    // counter("Rz.After DTF")++;

    //
    //Apply cut on MC truth
    //


    //NOTE!!!!BELOW TAKEN OUT FOR B TESTS - need to be put back in
    // if (motherPID == 0.0) continue;

    counter("S.After TruthCheck")++;

    Gaudi::LorentzVector jpsi4Mom = Jpsi->momentum();
    Gaudi::LorentzVector mu4Mom   = muPlus->momentum();
    ROOT::Math::Boost boost(jpsi4Mom.BoostToCM());
    const Gaudi::XYZVector boostedMu4Mom   = (boost(mu4Mom)).Vect().unit();
    const Gaudi::XYZVector boostedJpsi4Mom = jpsi4Mom.Vect().unit();
    double cosTheta = fabs(boostedMu4Mom.Dot(boostedJpsi4Mom));
    double Theta = acos(boostedMu4Mom.Dot(boostedJpsi4Mom));
    

    Gaudi::LorentzVector p1 = muPlus->momentum();
    Gaudi::LorentzVector p2 = muMinus->momentum();

    double pp1      = sqrt(p1.X()*p1.X() + p1.Y()*p1.Y());
    double pp2      = sqrt(p2.X()*p2.X() + p2.Y()*p2.Y());
    double cosPhi   = (p1.X()*p2.X() + p1.Y()*p2.Y())/(pp1*pp2);
    double deltaPhi = acos(cosPhi);

    // FillDataPid(allJpsi,muPlus,prefixMuPlus);

    double PlusProbNNMu  = muProtoPlus->info(LHCb::ProtoParticle::ProbNNmu,    -1000);
    double MinusProbNNMu = muProtoMinus->info(LHCb::ProtoParticle::ProbNNmu,    -1000);



    efftuple->fill(
"Jpsimass,Jpsipt,Rap,AbsCosTheta,Theta,deltaPhi,ghostMuPlus,ghostMuMinus,pseudorapMuPlus,pseudorapMuMinus,NNMuPlus,NNMuMinus",
                   Jpsi->measuredMass(),Jpsi->pt(),Jpsi->momentum().Rapidity(),
                   cosTheta,Theta,deltaPhi,muTrackPlus->ghostProbability(),muTrackMinus->ghostProbability(),
muTrackPlus->pseudoRapidity(),muTrackMinus->pseudoRapidity(),PlusProbNNMu,MinusProbNNMu);

    efftuple->write();
    
    counter("T.FINAL COUNT")++;

    const LHCb::Particle jpsi = *Jpsi;
    jpsiVec.push_back(jpsi);
    
   
  }

  /*
   * Bail out if less than 2 J/psi.
   */
 
  int size = jpsiVec.size();
  std::cout << "number of Jpsi"  << jpsiVec.size()   << std::endl;
  
  if (size < 2) return sc;
  counter("MORE THAN ONE Jpsi!")++;

  for (int j = 1; j <= (size-1) ; j=(j+1)){
    for (int k = (j+1); k <= size; k=(k+1) ){        
      bool acceptDoubleJpsi = true;
 
      const LHCb::Particle *jpsi1 = NULL;
      const LHCb::Particle *jpsi2 = NULL;
      counter("Z 1")++;

      if ((jpsiVec[j-1]).pt() > (jpsiVec[k-1]).pt()){
        jpsi1 = (&jpsiVec[j-1]);
        jpsi2 = (&jpsiVec[k-1]);
      } else{
        jpsi1 = (&jpsiVec[k-1]);
        jpsi2 = (&jpsiVec[j-1]);
      }
      
      /*
       * Make sure all muons are different.
       */
      counter("Z 2")++;
      const LHCb::Particle *muPlus1  = GetParticles(jpsi1,  1,  m_pidMuon);
      if (muPlus1  == NULL) continue;
      counter("Z 2a")++;
      const LHCb::Particle *muMinus1 = GetParticles(jpsi1, -1, -m_pidMuon);
      if (muMinus1 == NULL) continue;
      counter("Z 2b")++;
      const LHCb::Particle *muPlus2  = GetParticles(jpsi2,  1,  m_pidMuon);
      if (muPlus2  == NULL) continue;
      counter("Z 2c")++;
      const LHCb::Particle *muMinus2 = GetParticles(jpsi2, -1, -m_pidMuon);
      if (muMinus2 == NULL) continue;
      counter("Z 3")++;
      StatusCode sc1 = DoubleJpsiCheckMuons(muPlus1, muMinus1, muPlus2, muMinus2);
      if (sc1.isFailure()){
        counter("Shared muons")++;
        continue;
      } 
      counter("Z 4")++;
      /*
       * Make a vertex with the 4 muons.
       */
 
      LHCb::Particle::ConstVector fourMuonsForFit;
      fourMuonsForFit.push_back(muPlus1);
      fourMuonsForFit.push_back(muMinus1);
      fourMuonsForFit.push_back(muPlus2);
      fourMuonsForFit.push_back(muMinus2);
 
      LHCb::Vertex   fourMuonVertex;
      LHCb::Particle fourMuon(chiB2PID);
 
      const IParticleCombiner *combiner = particleCombiner();
      StatusCode scFit = combiner->combine(fourMuonsForFit, fourMuon,
                                           fourMuonVertex);
      if (!scFit){
        warning() << "Four muons fit errror. jpsi1 id = "
                  << jpsi1->particleID().pid()
                  << " jpsi2 id = "
                  << jpsi2->particleID().pid()
                  << endmsg;
        counter("Four muons fit error")++;
        continue;
      }
 
      if (fourMuonVertex.chi2() > maxFourMuonVtxChi2) acceptDoubleJpsi = false;
     
      /*
       * Make a fake particle decay to create the vertex for the 2 J/Psi.
       */
      counter("Z 5")++;
      LHCb::Vertex   motherVertex;
      LHCb::Particle mother(chiB2PID);
 
      LHCb::Particle::ConstVector particlesForFit;
 
      particlesForFit.push_back(jpsi1);
      particlesForFit.push_back(jpsi2);
 
      scFit = combiner->combine(particlesForFit, mother, motherVertex);
 
      if (!scFit){
        Warning("Double J/psi fit error").ignore();
        counter("Double J/psi fit error")++;
        continue;
      }
     
      /*
       * Run a DecayTreeFitter on the mother
       */
     
      const LHCb::VertexBase *pv = bestPV(&mother);
      if (!pv) continue;
      counter("Z 6")++;
       const ITrackStateProvider *stateProvider = m_stateProvider.empty() ?
        NULL : &(*m_stateProvider);
 
      DecayTreeFitter::Fitter dtFitter(mother, pv, stateProvider);
      //      DecayTreeFitter::Fitter dtFitter(mother, pv);
      dtFitter.fit();
     
      LHCb::DecayTree  motherTree = dtFitter.getFittedTree();
      LHCb::Particle  *fitMother  = motherTree.head();
     
      /*
       * Fill a ntuple with the data of all possible double J/Psi.
       */
      counter("All double Jpsi")++;
      //if (m_fillAllDoubleJpsi){ // AC
      if (true){
        Tuple allDoubleJpsi = nTuple("AllDoubleJpsi");
 
      }
 
      /*
       * Apply the cut on the DecayTreeFitter Chi2/Ndof
       */
     
      if (fitMother->endVertex()->chi2PerDoF() > m_maxDTFChi2Ndof)
        acceptDoubleJpsi = false;
     
      /*
       * Skip non selected double J/Psi.
       */
 
      if (!acceptDoubleJpsi) continue;
 
      /*
       * Fill a ntuple with the data of all selected double J/Psi.
       */
     
      if (m_fillSelectedDoubleJpsi){
        Tuple selDoubleJpsi = nTuple("SelectedDoubleJpsi");
 
        // FillAnalysisDoubleJpsiMuMu(selDoubleJpsi, jpsi1, muPlus1, muMinus1,
        //                           jpsi2, muPlus2, muMinus2, &fourMuon, fitMother);
       
        //selDoubleJpsi->write();
      }
    }
    
  }
  

  return sc;  
}





//=============================================================================
//  Finalize
//=============================================================================
StatusCode DoubleJpsiSel::finalize() {

  if (msgLevel(MSG::DEBUG)) debug() << "==> Finalize" << endmsg;
  return DaVinciTupleAlgorithm::finalize();

}


//=====================================================================================================
// MCTruth efficiency of single J/psi
//======================================================================================================


StatusCode DoubleJpsiSel::MCTruth(){
  Tuple efftuple = nTuple("effitupletruth");
  Tuples::Tuple allMCJpsi = nTuple("AllMCTrueJpsi");
  StatusCode sc = StatusCode::FAILURE;
  //Tuple allMCJpsi = nTuple("AllMCTrueJpsi");
  const std::string prefixJpsi    = "Jpsi";
  const std::string prefixMuPlus  = "MuPlus";
  const std::string prefixMuMinus = "MuMinus";
  // Gaudi::LorentzVector p2;


  LHCb::MCParticles* mcParts;
  if (exist<LHCb::MCParticles>(LHCb::MCParticleLocation::Default)){
    mcParts = get<LHCb::MCParticles>(LHCb::MCParticleLocation::Default);
  } else{
    counter("No default location")++;
    return sc;
  }



  LHCb::MCParticles::const_iterator iPart;
  for (iPart = mcParts->begin(); iPart != mcParts->end(); ++iPart){
    counter("TruthA. 1st MC loop")++;

    const LHCb::MCParticle *jpsi = *iPart;

    //===============
    const LHCb::MCVertex * jpsiOrgVertex = jpsi->originVertex();

    bool fromPV = false;
    bool oneFeed = false;
    bool twoFeed = false;
    bool nonPrompt = false;
    bool isSignal = false;
    bool isJpsiAndMu = false;
    
      
    //
    //get rid of non-prompt (i.e. 1st or second feed-down)
    //
    if (jpsiOrgVertex->isPrimary()){
      fromPV = true;
      // counter(" zzz. TruthPV  ")++;
    }
    else if (jpsi->mother()->originVertex()->isPrimary()){
      oneFeed = true;
      // counter(" zzz. TruthoneFeed  ")++;
    }
    else if (jpsi->mother()->mother()->originVertex()->isPrimary()){
      twoFeed = true;
      //counter(" zzz. Truth twoFeed  ")++;
    }
    else{
      //counter("TruthB. Truth 1st reject non-prompt  ")++;
      //continue;
      nonPrompt = true;
    }
    counter("TruthB. After 1st reject non-prompt  ")++;
    //
    //feeddown from charmonium family and cc-bar fragments
    //
    if (oneFeed){
      if (jpsi->mother()->particleID().pid()==(9900441))oneFeed=true; //cc-bar[1S08]
      else if (jpsi->mother()->particleID().pid()==(9900443)) oneFeed=true; ////cc-bar[3S08]
      else if(jpsi->mother()->particleID().pid()==(9910441)) oneFeed = true;//cc-bar[3P08]
      //else if(Jpsi1->particleID().pid()==(91)) fromjpsi = true //
      else if ((jpsi)->mother()->particleID().pid()==(100443)) oneFeed=true; //psi(2s)
      else if ((jpsi)->mother()->particleID().pid()==(445)) oneFeed=true; //Chic2(1P)
      else if ((jpsi)->mother()->particleID().pid()==(20443)) oneFeed=true; //Chic1(1P)
      else if ((jpsi)->mother()->particleID().pid()==(10441)) oneFeed=true; //Chic0(1P)
      else if ((jpsi)->mother()->particleID().pid()==(100445)) oneFeed=true; //chic2(2P)
      else if ((jpsi)->mother()->particleID().pid()==(10443)) oneFeed=true; //hc(1P)
      else if ((jpsi)->mother()->particleID().pid()==(30443)) oneFeed=true; //psi(3770) ???
      else {
        //counter ("TruthC. Truth rejected onefeed")++;
        //continue;
        oneFeed = false;
      }
    }
    
    counter ("TruthC. after rejected onefeed")++;

    //
    //get two feed from psi(2s) and potentially cc(bar) fragments
    //
    if (twoFeed){
      if ((jpsi)->mother()->mother()->particleID().pid()==(100443))twoFeed = true; //psi(2s)
      else  if (jpsi->mother()->mother()->particleID().pid()==(30443))twoFeed=true; //cc-bar[1S08]
      else  if (jpsi->mother()->mother()->particleID().pid()==(9900441))twoFeed=true; //cc-bar[1S08]
      else if (jpsi->mother()->mother()->particleID().pid()==(9900443)) twoFeed=true; ////cc-bar[3S08]
      else if(jpsi->mother()->mother()->particleID().pid()==(9910441)) twoFeed = true;//cc-bar[3P08]
      else {
        //counter ("TruthD. Truth rejected twofeed")++;
        //  continue;
        twoFeed = false;
      }
    }
    counter ("TruthD. after rejected twofeed")++;
    //if (fromPV == false && oneFeed ==false && twoFeed == false){
      //counter ("zzz. Truth 2nd rejected non-prompt")++;
    //  continue; 
    //  }
    counter ("TruthE. after rejected non-prompt")++;
    //===============
    
    if(jpsi->particleID().pid()!= 443){
      // std::cout << "pid of Jpsi not a Jpsi: "  << jpsi->particleID().pid()  << std::endl;
      //counter("Jpsi not a Jpsi")++;
      continue;
      //isJpsiAndMu = true;
    }
    else isJpsiAndMu = true;
    
    counter("TruthF. Jpsi not a Jpsi")++;

    
    const SmartRefVector<LHCb::MCVertex> &vertices = jpsi->endVertices();
    SmartRefVector<LHCb::MCVertex>::const_iterator iVert;
    for (iVert = vertices.begin(); iVert != vertices.end(); ++iVert){
      
      const LHCb::MCParticle *muPlus  = NULL;
      const LHCb::MCParticle *muMinus = NULL;

      const SmartRefVector<LHCb::MCParticle> &daughters = (*iVert)->products();
      SmartRefVector<LHCb::MCParticle>::const_iterator iDau;
      for (iDau = daughters.begin(); iDau != daughters.end(); ++iDau){
        if ((*iDau)->particleID().pid() ==  13) muPlus  = *iDau;
        if ((*iDau)->particleID().pid() == -13) muMinus = *iDau;
      }
      
      if (!muPlus)  continue;
      if (!muMinus) continue;
      if(isJpsiAndMu)if(!muPlus && !muMinus) isJpsiAndMu = false;
         
      counter("TruthG. after MC true muons")++;
      
      // if(isJpsiAndMu && !nonPrompt && (fromPV || oneFeed || twoFeed )) isSignal = true;
      //removed above as just contunue out of loop if not JpsiAndMu
      if(!nonPrompt && (fromPV || oneFeed || twoFeed )) isSignal = true;

      counter("TruthGv.")++;
      //
      //get variables to read into tuple
      //
      double rap =  jpsi->momentum().Rapidity();
      Gaudi::LorentzVector JpsiTrue4Mom  =  jpsi->momentum();
      Gaudi::LorentzVector muTrue4Mom  = muPlus->momentum();
      counter("TruthGw.")++;
      ROOT::Math::Boost   boost ( JpsiTrue4Mom.BoostToCM() ) ;
      counter("TruthGx.")++;
      const Gaudi::XYZVector boostedMu = (boost( muTrue4Mom )).Vect().unit();
      const Gaudi::XYZVector boostedJpsi = JpsiTrue4Mom.Vect().unit();
      counter("TruthGy.")++;
      double cosTTrue = boostedMu.Dot(boostedJpsi) ;
      double abscosT = fabs(cosTTrue);
      //double motherPID = (*mother)->pdg_id();
      double ptrans = JpsiTrue4Mom.Pt();
       counter("TruthGz.")++;
      double TruthTheta = acos(boostedMu.Dot(boostedJpsi));
      
      Gaudi::LorentzVector p1 = muPlus->momentum();
      Gaudi::LorentzVector p2 = muMinus->momentum();
      double pp1      = sqrt(p1.X()*p1.X() + p1.Y()*p1.Y());
      double pp2      = sqrt(p2.X()*p2.X() + p2.Y()*p2.Y());
      double cosPhi   = (p1.X()*p2.X() + p1.Y()*p2.Y())/(pp1*pp2);
      double TruthPhi = acos(cosPhi);
      
      counter("TruthH")++;
      
      plot3D(ptrans,rap,abscosT,"Truth","Truth",
             0.0,10000.0,2.0,4.5,0.0,1.0,500,0.5,0.1);
      efftuple->fill(
                     "TruthJpsimass,TruthJpsipt,TruthRapidity,TruhtAbsCosTheta,TruthTheta,TruthPhi",
                     jpsi->virtualMass(),JpsiTrue4Mom.Pt(),
                     jpsi->momentum().Rapidity(),abscosT,TruthTheta,TruthPhi);
      
      efftuple->write();
      sc = StatusCode::SUCCESS;


      // Tuple allMCJpsi = nTuple("AllMCTrueJpsi");
 
      // FillDataEventNtuple(allMCJpsi);
      FillDataEventNtuple(allMCJpsi);
      FillMCComposite(allMCJpsi, jpsi,    prefixJpsi);
      /*
      int pdgPid = 0;
      const LHCb::MCParticle *motherTop = jpsi->mother();
      while (motherTop != NULL){
        pdgPid = motherTop->particleID().pid();
        motherTop = motherTop->mother();
      }
      allMCJpsi->column(prefixJpsi + "_TOP_ANCESTOR_PID", pdgPid);
      */
      counter("TruthI")++;
      
      FillMCAncestor(allMCJpsi,  jpsi,    prefixJpsi);
      
      FillMCParticle(allMCJpsi,  muPlus,  prefixMuPlus);
      FillMCParticle(allMCJpsi,  muMinus, prefixMuMinus);
 
      FillMCPolarization(allMCJpsi, jpsi, muPlus, prefixJpsi);
      allMCJpsi->column("isSignal", isSignal);
      allMCJpsi->column("Jpsi_Theta",double(TruthTheta));
      allMCJpsi->column("Jpsi_Phi",double(TruthPhi));
      allMCJpsi->write();
      counter("TruthJ")++;
      
    }
  }
  return sc;
}



StatusCode DoubleJpsiSel::MCTruthTest(){
  Tuple efftuple = nTuple("effitupletruth");
  StatusCode sc = StatusCode::FAILURE;
  //Tuple allMCJpsi = nTuple("AllMCTrueJpsi");
  const std::string prefixJpsi    = "Jpsi";
  const std::string prefixMuPlus  = "MuPlus";
  const std::string prefixMuMinus = "MuMinus";
  // Gaudi::LorentzVector p2;


  LHCb::MCParticles* mcParts;
  if (exist<LHCb::MCParticles>(LHCb::MCParticleLocation::Default)){
    mcParts = get<LHCb::MCParticles>(LHCb::MCParticleLocation::Default);
  } else{
    counter("No default location")++;
    return sc;
  }



  LHCb::MCParticles::const_iterator iPart;
  for (iPart = mcParts->begin(); iPart != mcParts->end(); ++iPart){
    counter("TruthA. 1st MC loop")++;

    const LHCb::MCParticle *jpsi = *iPart;

    //===============
    const LHCb::MCVertex * jpsiOrgVertex = jpsi->originVertex();

    bool fromPV = false;
    bool oneFeed = false;
    bool twoFeed = false;
        
    if(jpsi->particleID().pid()!= 443){
      // std::cout << "pid of Jpsi not a Jpsi: "  << jpsi->particleID().pid()  << std::endl;
      //counter("Jpsi not a Jpsi")++;
      continue;
    }
    counter("TruthA. Jpsi not a Jpsi")++;

      
    const SmartRefVector<LHCb::MCVertex> &vertices = jpsi->endVertices();
    SmartRefVector<LHCb::MCVertex>::const_iterator iVert;
    for (iVert = vertices.begin(); iVert != vertices.end(); ++iVert){
      
      const LHCb::MCParticle *muPlus  = NULL;
      const LHCb::MCParticle *muMinus = NULL;

      const SmartRefVector<LHCb::MCParticle> &daughters = (*iVert)->products();
      SmartRefVector<LHCb::MCParticle>::const_iterator iDau;
      for (iDau = daughters.begin(); iDau != daughters.end(); ++iDau){
        if ((*iDau)->particleID().pid() ==  13) muPlus  = *iDau;
        if ((*iDau)->particleID().pid() == -13) muMinus = *iDau;
      }
      
      if (!muPlus)  continue;
      if (!muMinus) continue;


      counter("TruthB. after MC true muons")++;


    //
    //get rid of non-prompt (i.e. 1st or second feed-down)
    //
    if (jpsiOrgVertex->isPrimary()){
      fromPV = true;
      // counter(" zzz. TruthPV  ")++;
    }
    else if (jpsi->mother()->originVertex()->isPrimary()){
      oneFeed = true;
      // counter(" zzz. TruthoneFeed  ")++;
    }
    else if (jpsi->mother()->mother()->originVertex()->isPrimary()){
      twoFeed = true;
      //counter(" zzz. Truth twoFeed  ")++;
    }
    else{
      //counter("TruthB. Truth 1st reject non-prompt  ")++;
      continue;
    }
    counter("TruthBz. After 1st reject non-prompt  ")++;
    //
    //feeddown from charmonium family and cc-bar fragments
    //
    if (oneFeed){
      if (jpsi->mother()->particleID().pid()==(9900441))oneFeed=true; //cc-bar[1S08]
      else if (jpsi->mother()->particleID().pid()==(9900443)) oneFeed=true; ////cc-bar[3S08]
      else if(jpsi->mother()->particleID().pid()==(9910441)) oneFeed = true;//cc-bar[3P08]
      //else if(Jpsi1->particleID().pid()==(91)) fromjpsi = true //
      else if ((jpsi)->mother()->particleID().pid()==(100443)) oneFeed=true; //psi(2s)
      else if ((jpsi)->mother()->particleID().pid()==(445)) oneFeed=true; //Chic2(1P)
      else if ((jpsi)->mother()->particleID().pid()==(20443)) oneFeed=true; //Chic1(1P)
      else if ((jpsi)->mother()->particleID().pid()==(10441)) oneFeed=true; //Chic0(1P)
      else if ((jpsi)->mother()->particleID().pid()==(100445)) oneFeed=true; //chic2(2P)
      else if ((jpsi)->mother()->particleID().pid()==(10443)) oneFeed=true; //hc(1P)
      else if ((jpsi)->mother()->particleID().pid()==(30443)) oneFeed=true; //psi(3770) ???
      else {
        //counter ("TruthC. Truth rejected onefeed")++;
        continue;
      }
    }
    
    counter ("TruthC. after rejected onefeed")++;

    //
    //get two feed from psi(2s) and potentially cc(bar) fragments
    //
    //const LHCb::MCParticle * twoFeedOrigin = (jpsi)->mother()->mother();
    //   int twoFeedOrigPid = (jpsi)->mother()->mother()->particleID().pid();
    if (twoFeed){
      if ((jpsi)->mother()->mother()->particleID().pid()==(100443))twoFeed = true; //psi(2s)
      else  if (jpsi->mother()->mother()->particleID().pid()==(9900441))twoFeed=true; //cc-bar[1S08]
      else if (jpsi->mother()->mother()->particleID().pid()==(9900443)) twoFeed=true; ////cc-bar[3S08]
      else if(jpsi->mother()->mother()->particleID().pid()==(9910441)) twoFeed = true;//cc-bar[3P08]
      //  else if(twoFeedOrigPid == 445 || twoFeedOrigPid == 20443 ||twoFeedOrigPid ==10441 
      //|| twoFeedOrigPid ==100445 || twoFeedOrigPid ==10443 ||twoFeedOrigPid == 30443 ) {
      // else{ 
      else if(jpsi->mother()->mother()->particleID().pid() == 445 || jpsi->mother()->mother()->particleID().pid() == 20443 || 
              jpsi->mother()->mother()->particleID().pid() == 10441 || jpsi->mother()->mother()->particleID().pid() == 100445 ||
              jpsi->mother()->mother()->particleID().pid() == 10443 || jpsi->mother()->mother()->particleID().pid() == 30443)
      {
        counter ("TruthD. Truth *rejected* twofeed charm at top")++;
        std::cout << " origin " << (jpsi)->mother()->mother()->particleID().pid()<< " -> " << 
          (jpsi)->mother()->particleID().pid() <<  " -> " << (jpsi)->particleID().pid()  <<std::endl;
        continue;
      }
    }
    counter ("TruthD. after rejected twofeed")++;
    if (fromPV == false && oneFeed ==false && twoFeed == false){
      //counter ("zzz. Truth 2nd rejected non-prompt")++;
      continue; 
      }
    counter ("TruthE. after rejected non-prompt")++;
    
      //
      //get variables to read into tuple
      //
      double rap =  jpsi->momentum().Rapidity();
      Gaudi::LorentzVector JpsiTrue4Mom  =  jpsi->momentum();
      Gaudi::LorentzVector muTrue4Mom  = muPlus->momentum(); 
      ROOT::Math::Boost   boost ( JpsiTrue4Mom.BoostToCM() ) ;
      const Gaudi::XYZVector boostedMu = (boost( muTrue4Mom )).Vect().unit();
      const Gaudi::XYZVector boostedJpsi = JpsiTrue4Mom.Vect().unit();
      double cosTTrue = boostedMu.Dot(boostedJpsi) ;
      double abscosT = fabs(cosTTrue);
      //double motherPID = (*mother)->pdg_id();
      double ptrans = JpsiTrue4Mom.Pt();
      
      double TruthTheta = acos(boostedMu.Dot(boostedJpsi));
      
      Gaudi::LorentzVector p1 = muPlus->momentum();
      Gaudi::LorentzVector p2 = muMinus->momentum();
      double pp1      = sqrt(p1.X()*p1.X() + p1.Y()*p1.Y());
      double pp2      = sqrt(p2.X()*p2.X() + p2.Y()*p2.Y());
      double cosPhi   = (p1.X()*p2.X() + p1.Y()*p2.Y())/(pp1*pp2);
      double TruthPhi = acos(cosPhi);
      
      
      plot3D(ptrans,rap,abscosT,"Truth","Truth",
             0.0,10000.0,2.0,4.5,0.0,1.0,500,0.5,0.1);
      efftuple->fill(
                     "TruthJpsimass,TruthJpsipt,TruthRapidity,TruhtAbsCosTheta,TruthTheta,TruthPhi",
                     jpsi->virtualMass(),JpsiTrue4Mom.Pt(),
                     jpsi->momentum().Rapidity(),abscosT,TruthTheta,TruthPhi);
      
      efftuple->write();
      sc = StatusCode::SUCCESS;


      Tuple allMCJpsi = nTuple("AllMCTrueJpsi");
 
      // FillDataEventNtuple(allMCJpsi);
      FillDataEventNtuple(allMCJpsi);
      FillMCComposite(allMCJpsi, jpsi,    prefixJpsi);
      FillMCAncestor(allMCJpsi,  jpsi,    prefixJpsi);
      
      FillMCParticle(allMCJpsi,  muPlus,  prefixMuPlus);
      FillMCParticle(allMCJpsi,  muMinus, prefixMuMinus);
 
      FillMCPolarization(allMCJpsi, jpsi, muPlus, prefixJpsi);
      allMCJpsi->column("Jpsi_Theta",double(TruthTheta));
      allMCJpsi->column("Jpsi_Phi",double(TruthPhi));
      allMCJpsi->write();
      
      
    }
  }
  return sc;
}




//=====================================================================================================
// MCTruth efficiency of 2x J/psi
//======================================================================================================


StatusCode DoubleJpsiSel::MCTruthDoubleJpsi(){

  const LHCb::ParticleID psi2SPID(m_pidPsi2S);
  const LHCb::ParticleID chiB2PID(m_pidChiB21P);

  // std::vector<const LHCb::MCParticle> jpsiVecMC;
  std::vector<const LHCb::MCParticle *> jpsiVec;
  //// std::vector<const LHCb::VertexBase*> pvvec;
  //std::vector<const LHCb::MCParticle*>   muonvec;

  Tuple efftuple = nTuple("effitupletruth");
  StatusCode sc = StatusCode::FAILURE;
  

  const std::string prefixFourMuon = "FourMuons";
  const std::string prefixJpsi1    = "JPsi_1";
  const std::string prefixJpsi2    = "JPsi_2";
  const std::string prefixMother   = "DoubleJpsi";
  const std::string prefixMuPlus1  = "MuPlus1";
  const std::string prefixMuPlus2  = "MuPlus2";
  const std::string prefixMuMinus1 = "MuMinus1";
  const std::string prefixMuMinus2 = "MuMinus2";

  //Tuple allMCJpsi = nTuple("AllMCTrueJpsi");
  const std::string prefixJpsi    = "Jpsi";
  const std::string prefixMuPlus  = "MuPlus";
  const std::string prefixMuMinus = "MuMinus";
  // Gaudi::LorentzVector p2;


  LHCb::MCParticles* mcParts;
  if (exist<LHCb::MCParticles>(LHCb::MCParticleLocation::Default)){
    mcParts = get<LHCb::MCParticles>(LHCb::MCParticleLocation::Default);
  } else{
    counter("No default location")++;
    return sc;
  }



  LHCb::MCParticles::const_iterator iPart;
  for (iPart = mcParts->begin(); iPart != mcParts->end(); ++iPart){
    counter("TruthA. 1st MC loop")++;

    const LHCb::MCParticle *jpsi = *iPart;


    //LHCb::MCParticle *jpsi = *iPart;
    //===============
    const LHCb::MCVertex * jpsiOrgVertex = jpsi->originVertex();
    // LHCb::MCVertex * jpsiOrgVertex = jpsi->originVertex();

    bool fromPV = false;
    bool oneFeed = false;
    bool twoFeed = false;
    
    
    //
    //get rid of non-prompt (i.e. 1st or second feed-down)
    //
    if (jpsiOrgVertex->isPrimary()){
      fromPV = true;
      // counter(" zzz. TruthPV  ")++;
    }
    else if (jpsi->mother()->originVertex()->isPrimary()){
      oneFeed = true;
      // counter(" zzz. TruthoneFeed  ")++;
    }
    else if (jpsi->mother()->mother()->originVertex()->isPrimary()){
      twoFeed = true;
      //counter(" zzz. Truth twoFeed  ")++;
    }
    else{
      //counter("TruthB. Truth 1st reject non-prompt  ")++;
      continue;
    }
    counter("TruthB. After 1st reject non-prompt  ")++;
    //
    //feeddown from charmonium family and cc-bar fragments
    //
    if (oneFeed){
      if (jpsi->mother()->particleID().pid()==(9900441))oneFeed=true; //cc-bar[1S08]
      else if (jpsi->mother()->particleID().pid()==(9900443)) oneFeed=true; ////cc-bar[3S08]
      else if(jpsi->mother()->particleID().pid()==(9910441)) oneFeed = true;//cc-bar[3P08]
      //else if(Jpsi1->particleID().pid()==(91)) fromjpsi = true //
      else if ((jpsi)->mother()->particleID().pid()==(100443)) oneFeed=true; //psi(2s)
      else if ((jpsi)->mother()->particleID().pid()==(445)) oneFeed=true; //Chic2(1P)
      else if ((jpsi)->mother()->particleID().pid()==(20443)) oneFeed=true; //Chic1(1P)
      else if ((jpsi)->mother()->particleID().pid()==(10441)) oneFeed=true; //Chic0(1P)
      else if ((jpsi)->mother()->particleID().pid()==(100445)) oneFeed=true; //chic2(2P)
      else if ((jpsi)->mother()->particleID().pid()==(100443)) oneFeed=true; //hc(1P)
      else if ((jpsi)->mother()->particleID().pid()==(30443)) oneFeed=true; //psi(3770) ???
      else {
        //counter ("TruthC. Truth rejected onefeed")++;
        continue;
      }
    }
    
    counter ("TruthC. after rejected onefeed")++;

    //
    //get two feed from psi(2s) and potentially cc(bar) fragments
    //
    if (twoFeed){
      if ((jpsi)->mother()->mother()->particleID().pid()==(100443))twoFeed = true; //psi(2s)
      else  if (jpsi->mother()->mother()->particleID().pid()==(9900441))twoFeed=true; //cc-bar[1S08]
      else if (jpsi->mother()->mother()->particleID().pid()==(9900443)) twoFeed=true; ////cc-bar[3S08]
      else if(jpsi->mother()->mother()->particleID().pid()==(9910441)) twoFeed = true;//cc-bar[3P08]
      else {
        //counter ("TruthD. Truth rejected twofeed")++;
        continue;
      }
    }
    counter ("TruthD. after rejected twofeed")++;
    if (fromPV == false && oneFeed ==false && twoFeed == false){
      //counter ("zzz. Truth 2nd rejected non-prompt")++;
      continue; 
      }
    counter ("TruthE. after rejected non-prompt")++;
    //===============
    
    if(jpsi->particleID().pid()!= 443){
      // std::cout << "pid of Jpsi not a Jpsi: "  << jpsi->particleID().pid()  << std::endl;
      //counter("Jpsi not a Jpsi")++;
      continue;
    }
    counter("TruthF. Jpsi not a Jpsi")++;
    
    
    const SmartRefVector<LHCb::MCVertex> &vertices = jpsi->endVertices();
    //SmartRefVector<LHCb::MCVertex> &vertices = jpsi->endVertices();
    SmartRefVector<LHCb::MCVertex>::const_iterator iVert;
    for (iVert = vertices.begin(); iVert != vertices.end(); ++iVert){
      
      const LHCb::MCParticle *muPlus  = NULL;
      const LHCb::MCParticle *muMinus = NULL;
      //LHCb::MCParticle *muPlus  = NULL;
      //LHCb::MCParticle *muMinus = NULL;

      const SmartRefVector<LHCb::MCParticle> &daughters = (*iVert)->products();
      // SmartRefVector<LHCb::MCParticle> &daughters = (*iVert)->products();
      SmartRefVector<LHCb::MCParticle>::const_iterator iDau;
      for (iDau = daughters.begin(); iDau != daughters.end(); ++iDau){
        if ((*iDau)->particleID().pid() ==  13) muPlus  = *iDau;
        if ((*iDau)->particleID().pid() == -13) muMinus = *iDau;
      }
      
      if (!muPlus)  continue;
      if (!muMinus) continue;


      counter("TruthG. after MC true muons")++;
              
      
      //
      //get variables to read into tuple
      //
      double rap =  jpsi->momentum().Rapidity();
      Gaudi::LorentzVector JpsiTrue4Mom  =  jpsi->momentum();
      Gaudi::LorentzVector muTrue4Mom  = muPlus->momentum(); 
      ROOT::Math::Boost   boost ( JpsiTrue4Mom.BoostToCM() ) ;
      const Gaudi::XYZVector boostedMu = (boost( muTrue4Mom )).Vect().unit();
      const Gaudi::XYZVector boostedJpsi = JpsiTrue4Mom.Vect().unit();
      double cosTTrue = boostedMu.Dot(boostedJpsi) ;
      double abscosT = fabs(cosTTrue);
      //double motherPID = (*mother)->pdg_id();
      double ptrans = JpsiTrue4Mom.Pt();
      
      double TruthTheta = acos(boostedMu.Dot(boostedJpsi));
      
      Gaudi::LorentzVector p1 = muPlus->momentum();
      Gaudi::LorentzVector p2 = muMinus->momentum();
      double pp1      = sqrt(p1.X()*p1.X() + p1.Y()*p1.Y());
      double pp2      = sqrt(p2.X()*p2.X() + p2.Y()*p2.Y());
      double cosPhi   = (p1.X()*p2.X() + p1.Y()*p2.Y())/(pp1*pp2);
      double TruthPhi = acos(cosPhi);
      
      
      plot3D(ptrans,rap,abscosT,"Truth","Truth",
             0.0,10000.0,2.0,4.5,0.0,1.0,500,0.5,0.1);
      efftuple->fill(
                     "TruthJpsimass,TruthJpsipt,TruthRapidity,TruhtAbsCosTheta,TruthTheta,TruthPhi",
                     jpsi->virtualMass(),JpsiTrue4Mom.Pt(),
                     jpsi->momentum().Rapidity(),abscosT,TruthTheta,TruthPhi);
      
      efftuple->write();
      sc = StatusCode::SUCCESS;


      Tuple allMCJpsi = nTuple("AllMCTrueJpsi");
 
      // FillDataEventNtuple(allMCJpsi);
      FillDataEventNtuple(allMCJpsi);
      FillMCComposite(allMCJpsi, jpsi,    prefixJpsi);
      FillMCAncestor(allMCJpsi,  jpsi,    prefixJpsi);
      
      FillMCParticle(allMCJpsi,  muPlus,  prefixMuPlus);
      FillMCParticle(allMCJpsi,  muMinus, prefixMuMinus);
 
      FillMCPolarization(allMCJpsi, jpsi, muPlus, prefixJpsi);
      allMCJpsi->column("Jpsi_Theta",double(TruthTheta));
      allMCJpsi->column("Jpsi_Phi",double(TruthPhi));
      allMCJpsi->write();

    
      jpsiVec.push_back(jpsi);
      
    }
  }

  
  //
  // Bail out if less than 2 J/psi.
  //
 
  int size = jpsiVec.size();
  //  std::cout << "number of Jpsi"  << jpsiVec.size()   << std::endl;
  
  if (size < 2) return sc;
  counter("TruthH. MORE THAN ONE Jpsi!")++;

  for (int j = 1; j <= (size-1) ; j=(j+1)){
    for (int k = (j+1); k <= size; k=(k+1) ){        
      bool acceptDoubleJpsi = true;
 
      const LHCb::MCParticle *jpsi1 = NULL;
      const LHCb::MCParticle *jpsi2 = NULL;
      counter("Z 1")++;

      if ((jpsiVec[j-1])->pt() > (jpsiVec[k-1])->pt()){
        jpsi1 = (jpsiVec[j-1]);
        jpsi2 = (jpsiVec[k-1]);
        // jpsi1 = (&jpsiVec[j-1]);
        //jpsi2 = (&jpsiVec[k-1]);
      } else{
        jpsi1 = (jpsiVec[k-1]);
        jpsi2 = (jpsiVec[j-1]);
        //jpsi1 = (&jpsiVec[k-1]);
        //jpsi2 = (&jpsiVec[j-1]);
      }
      ///$$$
      //
      //  Make sure all muons are different.
      //


      const LHCb::MCParticle *muPlus1  = NULL;
      const LHCb::MCParticle *muMinus1 = NULL;
      const LHCb::MCParticle *muPlus2  = NULL;
      const LHCb::MCParticle *muMinus2 = NULL;

      const SmartRefVector<LHCb::MCVertex> &vertices1 = jpsi1->endVertices();
      SmartRefVector<LHCb::MCVertex>::const_iterator iVert1;
      for (iVert1 = vertices1.begin(); iVert1 != vertices1.end(); ++iVert1){
      

        const SmartRefVector<LHCb::MCParticle> &daughtersJpsi1 = (*iVert1)->products();
  
        SmartRefVector<LHCb::MCParticle>::const_iterator iDauJpsi1;
        for (iDauJpsi1 = daughtersJpsi1.begin(); iDauJpsi1 != daughtersJpsi1.end(); ++iDauJpsi1){
          if ((*iDauJpsi1)->particleID().pid() ==  13)   muPlus1  = *iDauJpsi1;
          if ((*iDauJpsi1)->particleID().pid() == -13)  muMinus1 = *iDauJpsi1;
        } 
      }
      if (muPlus1  == NULL) continue;
      if (muMinus1 == NULL) continue;

      const SmartRefVector<LHCb::MCVertex> &vertices2 = jpsi2->endVertices();
      SmartRefVector<LHCb::MCVertex>::const_iterator iVert2;
      for (iVert2 = vertices2.begin(); iVert2 != vertices2.end(); ++iVert2){
      

        const SmartRefVector<LHCb::MCParticle> &daughtersJpsi2 = (*iVert2)->products();
  
        SmartRefVector<LHCb::MCParticle>::const_iterator iDauJpsi2;
        for (iDauJpsi2 = daughtersJpsi2.begin(); iDauJpsi2 != daughtersJpsi2.end(); ++iDauJpsi2){
          if ((*iDauJpsi2)->particleID().pid() ==  13)   muPlus2  = *iDauJpsi2;
          if ((*iDauJpsi2)->particleID().pid() == -13)  muMinus2 = *iDauJpsi2;
        } 
      }
      if (muPlus2  == NULL) continue;
      if (muMinus2 == NULL) continue;

      //
      //check muons
      //

      /*
      if (muPlus1->key()  == muPlus2->key()) {
        std::cout << "Cut first cut"  << std::endl;
        counter("Cut 1")++;
        return StatusCode::FAILURE;
      }
      
      if (muMinus1->key() == muMinus2->key()) {
        std::cout << "Cut second cut"  << std::endl;
        counter("Cut 2")++;
        return StatusCode::FAILURE;
      }
      
      */

      //
      //phi and theta variables
      //

      //
      //Jpsi1
      //
      double rap1 =  jpsi1->momentum().Rapidity();
      Gaudi::LorentzVector JpsiTrue4Mom1  =  jpsi1->momentum();
      Gaudi::LorentzVector muTrue4Mom1  = muPlus1->momentum(); 
      ROOT::Math::Boost   boost1 ( JpsiTrue4Mom1.BoostToCM() ) ;
      const Gaudi::XYZVector boostedMu1 = (boost1( muTrue4Mom1 )).Vect().unit();
      const Gaudi::XYZVector boostedJpsi1 = JpsiTrue4Mom1.Vect().unit();
      double cosTTrue1 = boostedMu1.Dot(boostedJpsi1) ;
      double abscosT1 = fabs(cosTTrue1);
      //double motherPID = (*mother)->pdg_id();
      double ptrans1 = JpsiTrue4Mom1.Pt();
      
      double TruthTheta1 = acos(boostedMu1.Dot(boostedJpsi1));
      
      Gaudi::LorentzVector p11 = muPlus1->momentum();
      Gaudi::LorentzVector p21 = muMinus1->momentum();
      double pp11      = sqrt(p11.X()*p11.X() + p11.Y()*p11.Y());
      double pp21      = sqrt(p21.X()*p21.X() + p21.Y()*p21.Y());
      double cosPhi1   = (p11.X()*p21.X() + p11.Y()*p21.Y())/(pp11*pp21);
      double TruthPhi1 = acos(cosPhi1);

      //
      //Jpsi2
      //
      double rap2 =  jpsi2->momentum().Rapidity();
      Gaudi::LorentzVector JpsiTrue4Mom2  =  jpsi2->momentum();
      Gaudi::LorentzVector muTrue4Mom2  = muPlus2->momentum(); 
      ROOT::Math::Boost   boost2 ( JpsiTrue4Mom2.BoostToCM() ) ;
      const Gaudi::XYZVector boostedMu2 = (boost2( muTrue4Mom2 )).Vect().unit();
      const Gaudi::XYZVector boostedJpsi2 = JpsiTrue4Mom2.Vect().unit();
      double cosTTrue2 = boostedMu2.Dot(boostedJpsi2) ;
      double abscosT2 = fabs(cosTTrue2);
      //double motherPID = (*mother)->pdg_id();
      double ptrans2 = JpsiTrue4Mom2.Pt();
      
      double TruthTheta2 = acos(boostedMu2.Dot(boostedJpsi2));
      
      Gaudi::LorentzVector p12 = muPlus2->momentum();
      Gaudi::LorentzVector p22 = muMinus2->momentum();
      double pp12      = sqrt(p12.X()*p12.X() + p12.Y()*p12.Y());
      double pp22      = sqrt(p22.X()*p22.X() + p22.Y()*p22.Y());
      double cosPhi2   = (p12.X()*p22.X() + p12.Y()*p22.Y())/(pp12*pp22);
      double TruthPhi2 = acos(cosPhi2);


      Tuple allDoubleMCJpsi = nTuple("AllDoubleMCTrueJpsi");
 
      FillDataEventNtuple(allDoubleMCJpsi);
      FillDataEventNtuple(allDoubleMCJpsi);
      FillMCComposite(allDoubleMCJpsi, jpsi1,    prefixJpsi1);
      FillMCAncestor(allDoubleMCJpsi,  jpsi1,    prefixJpsi1);
      FillMCParticle(allDoubleMCJpsi,  muPlus1,  prefixMuPlus1);
      FillMCParticle(allDoubleMCJpsi,  muMinus1, prefixMuMinus1);
      FillMCPolarization(allDoubleMCJpsi, jpsi1, muPlus1, prefixJpsi1);
      FillMCComposite(allDoubleMCJpsi, jpsi2,    prefixJpsi2);
      FillMCAncestor(allDoubleMCJpsi,  jpsi2,    prefixJpsi2);
      FillMCParticle(allDoubleMCJpsi,  muPlus2,  prefixMuPlus2);
      FillMCParticle(allDoubleMCJpsi,  muMinus2, prefixMuMinus2);
      FillMCPolarization(allDoubleMCJpsi, jpsi2, muPlus2, prefixJpsi2);

      allDoubleMCJpsi->column("Jpsi_1_Theta",double(TruthTheta1));
      allDoubleMCJpsi->column("Jpsi_1_Phi",double(TruthPhi1));
      allDoubleMCJpsi->column("Jpsi_2_Theta",double(TruthTheta2));
      allDoubleMCJpsi->column("Jpsi_2_Phi",double(TruthPhi2));

      allDoubleMCJpsi->write();

  


    }
    
  }

  return sc;
}





/*

StatusCode DoubleJpsiSel::MCTruth(){
  Tuple efftuple = nTuple("effitupletruth");
  StatusCode sc = StatusCode::FAILURE;
  //Tuple allMCJpsi = nTuple("AllMCTrueJpsi");
  const std::string prefixJpsi    = "Jpsi";
  const std::string prefixMuPlus  = "MuPlus";
  const std::string prefixMuMinus = "MuMinus";
  // Gaudi::LorentzVector p2;


  LHCb::MCParticles* mcParts;
  if (exist<LHCb::MCParticles>(LHCb::MCParticleLocation::Default)){
    mcParts = get<LHCb::MCParticles>(LHCb::MCParticleLocation::Default);
  } else{
    counter("No default location")++;
    return sc;
  }



  LHCb::MCParticles::const_iterator iPart;
  for (iPart = mcParts->begin(); iPart != mcParts->end(); ++iPart){
    counter("A. 1st MC loop")++;

    const LHCb::MCParticle *jpsi = *iPart;

    // Select only J/psi

    if (jpsi->particleID().pid() != 443) continue;
    counter("Found MC true J/psi")++;
    
    const SmartRefVector<LHCb::MCVertex> &vertices = jpsi->endVertices();
    SmartRefVector<LHCb::MCVertex>::const_iterator iVert;
    for (iVert = vertices.begin(); iVert != vertices.end(); ++iVert){
      
      const LHCb::MCParticle *muPlus  = NULL;
      const LHCb::MCParticle *muMinus = NULL;

      const SmartRefVector<LHCb::MCParticle> &daughters = (*iVert)->products();
      SmartRefVector<LHCb::MCParticle>::const_iterator iDau;
      for (iDau = daughters.begin(); iDau != daughters.end(); ++iDau){
        if ((*iDau)->particleID().pid() ==  13) muPlus  = *iDau;
        if ((*iDau)->particleID().pid() == -13) muMinus = *iDau;
      }
      
      if (!muPlus)  continue;
      if (!muMinus) continue;


      counter("Found MC true muons")++;
              
      
      //
      //get variables to read into tuple
      //
      double rap =  jpsi->momentum().Rapidity();
      Gaudi::LorentzVector JpsiTrue4Mom  =  jpsi->momentum();
      Gaudi::LorentzVector muTrue4Mom  = muPlus->momentum(); 
      ROOT::Math::Boost   boost ( JpsiTrue4Mom.BoostToCM() ) ;
      const Gaudi::XYZVector boostedMu = (boost( muTrue4Mom )).Vect().unit();
      const Gaudi::XYZVector boostedJpsi = JpsiTrue4Mom.Vect().unit();
      double cosTTrue = boostedMu.Dot(boostedJpsi) ;
      double abscosT = fabs(cosTTrue);
      //double motherPID = (*mother)->pdg_id();
      double ptrans = JpsiTrue4Mom.Pt();
      
      double TruthTheta = acos(boostedMu.Dot(boostedJpsi));
      
      Gaudi::LorentzVector p1 = muPlus->momentum();
      Gaudi::LorentzVector p2 = muMinus->momentum();
      double pp1      = sqrt(p1.X()*p1.X() + p1.Y()*p1.Y());
      double pp2      = sqrt(p2.X()*p2.X() + p2.Y()*p2.Y());
      double cosPhi   = (p1.X()*p2.X() + p1.Y()*p2.Y())/(pp1*pp2);
      double TruthPhi = acos(cosPhi);
      
      
      plot3D(ptrans,rap,abscosT,"Truth","Truth",
             0.0,10000.0,2.0,4.5,0.0,1.0,500,0.5,0.1);
      efftuple->fill(
                     "TruthJpsimass,TruthJpsipt,TruthRapidity,TruhtAbsCosTheta,TruthTheta,TruthPhi",
                     jpsi->virtualMass(),JpsiTrue4Mom.Pt(),
                     jpsi->momentum().Rapidity(),abscosT,TruthTheta,TruthPhi);
      
      efftuple->write();
      sc = StatusCode::SUCCESS;
      //
      //read into a tuple
      //
      //below commented out as not needed and would have to convert to new reco
      // if(m_fillMCTruthTuple){
      // FillMCTrueMuNtuple(allMCJpsi, mcpartp, prefixMuPlus);
      //FillMCTrueJpsiNtuple(allMCJpsi, mcpartjpsi, prefixJpsi);
      //FillMCTruePolarizationNtuple(allMCJpsi, mcpartjpsi, mcpartp, prefixJpsi);
      //allMCJpsi->write();
      //}
      
      
    }
  }
  return sc;
}


*/

//================================

// the below no longer works on the latest sim08 versions of MC - converted to the above
StatusCode DoubleJpsiSel::MCTruthOld(){ 
  
  
  Tuple efftuple = nTuple("effitupletruth");
  StatusCode sc = StatusCode::FAILURE;
  const LHCb::HepMCEvents* hepmc;
  Tuple allMCJpsi = nTuple("AllMCTrueJpsi");
  const std::string prefixJpsi    = "Jpsi";
  const std::string prefixMuPlus  = "MuPlus";
  const std::string prefixMuMinus = "MuMinus";
  // Gaudi::LorentzVector p2;
  if (exist<LHCb::HepMCEvents>(LHCb::HepMCEventLocation::Default)){
    hepmc = get<LHCb::HepMCEvents>(LHCb::HepMCEventLocation::Default);
  } 
  else {
    std::cout << "::MCTruth - no HepMCEvents" << std::endl;
    return sc;
  }
  
  for (LHCb::HepMCEvents::const_iterator ie = hepmc->begin(); hepmc->end() !=ie;++ie) {
    counter("A. 1st MC loop")++;
    const LHCb::HepMCEvent* event = *ie;
    const HepMC::GenEvent* gevt = event->pGenEvt();
    if (!gevt){
      warning() << "::MCTruth - GenEvent* is not defined" << endmsg;
      continue;
    }
    //
    //loop over particles
    //
    for(HepMC::GenEvent::particle_const_iterator ip=gevt->particles_begin();gevt->particles_end()!=ip;++ip){
      
      counter("B. 2nd  MC loop")++;
      const HepMC::GenParticle* mcpart = *ip;
      int pdg_id = mcpart->pdg_id();
      const LHCb::MCParticle* mcpartp = NULL;//mc true mu positive - defined here for use in cos theta loop etc
      const LHCb::MCParticle* mcpartn = NULL;//mc true mu negative - defined here for use in cos theta loop etc
      
      //
      //filter out J/psi
      //
      
      if (pdg_id == 443){
        HepMC::GenVertex * ev = mcpart->end_vertex();
        if (!ev){
          warning() << "::MCTruth - GenVertex is null" << endmsg;
          continue;
        }
        //  std::cout << "**********************" << std::endl;
        
        bool isrecp = false; 
        bool isrecn = false;
        counter("C. ")++;
        //
        //loop over decay vertex
        //
        // Gaudi::LorentzVector p2;
        for(HepMC::GenVertex::particles_out_const_iterator 
              iv=ev->particles_out_const_begin();ev->particles_out_const_end()!=iv;++iv){          
          int pdg_idmu = (*iv)->pdg_id();
          
          const IHepMC2MC::HepMC2MC* mctable;
          const IHepMC2MC* m_hepmc2mc;
          
          counter("D. ")++;

          if (pdg_idmu==13){ // muon
            // const IHepMC2MC::HepMC2MC* mctable;//====Variables to turn HepMC to (LHCb::MC) to use isReconstructible
            //const IHepMC2MC* m_hepmc2mc;
            m_hepmc2mc = tool<IHepMC2MC >("LoKi::HepMC2MC");
            if (m_hepmc2mc->hepMC2MC() != NULL) mctable = m_hepmc2mc->hepMC2MC();
            else return sc;
            
            HepMC::GenParticle* hepmcpart;
            hepmcpart = *iv;
            const IHepMC2MC::HepMC2MC::Range mclinksp = mctable->relations(hepmcpart);
            
            mcpartp = mclinksp.back().to();
         
          }
          counter("E. ")++;
            
          //======
          if (pdg_idmu==-13){ // anti-muon
            // const IHepMC2MC::HepMC2MC* mctablen;//====Variables to turn HepMC to (LHCb::MC) to use isReconstructible
            //const IHepMC2MC* m_hepmc2mc;
            m_hepmc2mc = tool<IHepMC2MC >("LoKi::HepMC2MC");
            if (m_hepmc2mc->hepMC2MC() != NULL) mctable = m_hepmc2mc->hepMC2MC();
            else return sc;
            
            HepMC::GenParticle* hepmcpartn;
            hepmcpartn = *iv;
            const IHepMC2MC::HepMC2MC::Range mclinksn = mctable->relations(hepmcpartn);
            
            mcpartn = mclinksn.back().to();
            //p2 = mcpartn->momentum();
            //std::cout << " mcpart*N* just made, pt is: " << p2  << std::endl;
          }
          
             
          counter("F. ")++;    
              //========
               

          if (mcpartp != NULL){ //====check if muon is reconstructible
            //Note: the isreconstructable is not actually needed now, so is commented out.  
            if (pdg_idmu == 13){
              isrecp = true;
              //  isrecp = m_recible->isReconstructibleAs((m_recible->reconstructible(mcpartp)),mcpartp);
              if (isrecp) counter("Y. MC Truth number of mu+ reconstructible")++; 
            } 
          }
          
          if (mcpartn != NULL){
            if (pdg_idmu == -13){
              isrecn = true;
           
              // isrecn = m_recible->isReconstructibleAs((m_recible->reconstructible(mcpartp)),mcpartp);
              if (isrecn) counter("Y. MC Truth number of mu- reconstructible")++; 
            }
          } 

          //  std::cout << "isrecp is: " << isrecp << " isrecn: " << isrecn  << "pid is: " << pdg_idmu << std::endl;
          
          counter("G. ")++;
          if(isrecp && isrecn){ 
            
            counter("Z. MC Truth J/psi from reconstructable muons")++;
            //
            //Prompt MC truth cut for DOUBLE Jpsi
            //
            for ( HepMC::GenVertex::particle_iterator mother = mcpart->production_vertex()->
                    particles_begin(HepMC::parents);
                  mother != mcpart->production_vertex()->particles_end(HepMC::parents); ++mother ){
              bool prompt = false; 
              const IHepMC2MC::HepMC2MC::Range mclinksjpsi = mctable->relations(mcpart);
              const LHCb::MCParticle* mcpartjpsi = mclinksjpsi.back().to();
              
              if(m_NoAccCutMC){
                bool inAcceptance = AccCheckMCTruth(mcpart);
                if(!inAcceptance) continue;
              }
               counter("H. ")++;
              // const LHCb::MCVertex * primaryVertex = mcpartjpsi->primaryVertex();
              //const LHCb::MCVertex * jpsiOrgVertex = mcpartjpsi->originVertex();
              //bool fromPV = false;
              //bool fromFeed = false;
              //if (jpsiOrgVertex->isPrimary()) fromPV = true;
              //else if((jpsiOrgVertex->mother()->originVertex())==0){
              //if ((jpsiOrgVertex->mother()->originVertex())->isPrimary()) fromFeed = true;
              //}
              //if (fromPV) prompt = true;
              //else if ((*mother)->pdg_id()==(100443) && fromFeed) prompt=true; //psi(2s)
              //else if ((*mother)->pdg_id()==(445) && fromFeed) prompt=true; //Chic2(1P)
              // else if ((*mother)->pdg_id()==(100445) && fromFeed) prompt=true; //chic2(2P)
              //else prompt = false;
              /////  if ((*mother)->pdg_id()==(9900441))prompt=true;
              ///// else if ((*mother)->pdg_id()==(9900443)) prompt=true;
              ///// else if((*mother)->pdg_id()==(9910441)) prompt = true;
              ///// else if((*mother)->pdg_id()==(91)) prompt = true;
              //  if (prompt) {
              
              counter("ZZ. FINAL MC Truth")++;
              
              
              //
              //get variables to read into tuple
              //
              double rap =  mcpartjpsi->momentum().Rapidity();
              Gaudi::LorentzVector JpsiTrue4Mom  =  mcpartjpsi->momentum();
              Gaudi::LorentzVector muTrue4Mom  = mcpartp->momentum(); 
              ROOT::Math::Boost   boost ( JpsiTrue4Mom.BoostToCM() ) ;
              const Gaudi::XYZVector boostedMu = (boost( muTrue4Mom )).Vect().unit();
              const Gaudi::XYZVector boostedJpsi = JpsiTrue4Mom.Vect().unit();
              double cosTTrue = boostedMu.Dot(boostedJpsi) ;
              double abscosT = fabs(cosTTrue);
              double motherPID = (*mother)->pdg_id();
              double ptrans = JpsiTrue4Mom.Pt();
              
              double TruthTheta = acos(boostedMu.Dot(boostedJpsi));
              
              Gaudi::LorentzVector p1 = mcpartp->momentum();
              Gaudi::LorentzVector p2 = mcpartn->momentum();
              double pp1      = sqrt(p1.X()*p1.X() + p1.Y()*p1.Y());
              double pp2      = sqrt(p2.X()*p2.X() + p2.Y()*p2.Y());
              double cosPhi   = (p1.X()*p2.X() + p1.Y()*p2.Y())/(pp1*pp2);
              double TruthPhi = acos(cosPhi);
              

              plot3D(ptrans,rap,abscosT,"Truth","Truth",
                     0.0,10000.0,2.0,4.5,0.0,1.0,500,0.5,0.1);
                efftuple->fill(
"TruthJpsimass,TruthJpsipt,TruthRapidity,TruhtAbsCosTheta,TruthTheta,TruthPhi, truthMotherPID",
                                 mcpartjpsi->virtualMass(),JpsiTrue4Mom.Pt(),
mcpartjpsi->momentum().Rapidity(),abscosT,TruthTheta,TruthPhi, motherPID);

                efftuple->write();
                  sc = StatusCode::SUCCESS;
                  //
                  //read into a tuple
                  //
                  if(m_fillMCTruthTuple){
                    FillMCTrueMuNtuple(allMCJpsi, mcpartp, prefixMuPlus);
                    FillMCTrueJpsiNtuple(allMCJpsi, mcpartjpsi, prefixJpsi);
                    FillMCTruePolarizationNtuple(allMCJpsi, mcpartjpsi, mcpartp, prefixJpsi);
                    allMCJpsi->write();
                  }
                  isrecp = false; 
                  isrecn = false;
                  // }
            }            
          }
        }          
      }
    }
  }
  //}
return sc;
}




//=================================================================================================
//
//For both muons simultaniously checks sign, that both muons aren't same and come from same Jpsi.
//==================================================================================================



double DoubleJpsiSel::CheckMCTruth(const LHCb::Track* part_trkP,const LHCb::Track* part_trkN)const{
  LHCb::Track2MC2D::Range mclinksN,mclinksP;//changed to const
  double mothercode = 0.0;
  const LHCb::Track2MC2D* table2d;
  counter("start of DoubleTrackToHepMC")++;
  if (exist<LHCb::Track2MC2D>(LHCb::Track2MCLocation::Default)){
    table2d = get<LHCb::Track2MC2D>(LHCb::Track2MCLocation::Default);}
  else{
    std::cout << " no HepMCEvents" << std::endl;
    return 0.0;    
  }
  
  mclinksP = table2d->relations(part_trkP);
  mclinksN = table2d->relations(part_trkN);
  
  for( LHCb::Track2MC2D::Range::iterator mclinkP=mclinksP.begin();mclinksP.end()!=mclinkP;++mclinkP){
    for( LHCb::Track2MC2D::Range::iterator mclinkN=mclinksN.begin();mclinksN.end()!=mclinkN;++mclinkN){
      const LHCb::MCParticle* mcpartP = mclinkP->to();
      const LHCb::MCParticle* mcpartN = mclinkN->to();     
      
      counter("TruthA. start of MCtruthcheck loop on muons")++;

      double mupluspid = mcpartP->particleID().pid();
      double muminuspid = mcpartN->particleID().pid();
      // std::cout << "muplus " << mupluspid <<" muminus " << muminuspid << std::endl;
      counter("TruthB. after check pid of muons")++;
      if(mcpartP->key() == mcpartN->key())return 0.0; ///check to see if same muon
       
      counter("TruthC. after same muon check ")++;
      //NOTE - I HAVE SWAPPED the below round - not sure, but it seems to have always been wrong
      //I guess there must be an error elsewhere - need to track down if can!!!

      //if (mupluspid!=13)return 0.0;//=====check to see if real muon of correct sign==========
      //if (muminuspid!=-13)return 0.0;
      
      if (mupluspid!=-13)return 0.0;//=====check to see if real muon of correct sign==========
      if (muminuspid!=13)return 0.0;


      counter("TruthD. second pid check??")++;

      const LHCb::MCParticle* Jpsi1 =  mcpartP->mother();
      const LHCb::MCParticle* Jpsi2 =  mcpartN->mother();

      const LHCb::MCVertex * jpsiOrgVertex = Jpsi1->originVertex();

      bool fromPV = false;
      bool oneFeed = false;
      bool twoFeed = false;
      
      
      //
      //get rid of non-prompt (i.e. 1st or second feed-down)
      //
      if (jpsiOrgVertex->isPrimary()){
        fromPV = true;
        //    counter(" bbb. number from PV  ")++;
      }
      else if (Jpsi1->mother()->originVertex()->isPrimary()){
        oneFeed = true;
        // counter(" bbbb. number of oneFeed  ")++;
      }
      else if (Jpsi1->mother()->mother()->originVertex()->isPrimary()){
        twoFeed = true;
        // counter(" bbbc. number twoFeed  ")++;
      }
      else{
        // counter(" zzzd. 1st reject non-prompt  ")++;
        return 0.0;
      }
      counter("TruthE. After 1st reject non-prompt  ")++;
      //
      //feeddown from charmonium family and cc-bar fragments
      //
      if (oneFeed){
        if (Jpsi1->mother()->particleID().pid()==(9900441))oneFeed=true; //cc-bar[1S08]
        else if (Jpsi1->mother()->particleID().pid()==(9900443)) oneFeed=true; ////cc-bar[3S08]
        else if(Jpsi1->mother()->particleID().pid()==(9910441)) oneFeed = true;//cc-bar[3P08]
        //else if(Jpsi1->particleID().pid()==(91)) fromjpsi = true //
        else if ((Jpsi1)->mother()->particleID().pid()==(100443)) oneFeed=true; //psi(2s)
        else if ((Jpsi1)->mother()->particleID().pid()==(445)) oneFeed=true; //Chic2(1P)
        else if ((Jpsi1)->mother()->particleID().pid()==(20443)) oneFeed=true; //Chic1(1P)
        else if ((Jpsi1)->mother()->particleID().pid()==(10441)) oneFeed=true; //Chic0(1P)
        else if ((Jpsi1)->mother()->particleID().pid()==(100445)) oneFeed=true; //chic2(2P)
        else if ((Jpsi1)->mother()->particleID().pid()==(10443)) oneFeed=true; //hc(1P)
        else if ((Jpsi1)->mother()->particleID().pid()==(30443)) oneFeed=true; //psi(3770) ???
        else {
          // counter ("zzze. rejected onefeed")++;
          return 0.0;
        }
      }
      counter ("TruthF. after rejected onefeed")++;
      //
      //get two feed from psi(2s) and potentially cc(bar) fragments
      //
      if (twoFeed){
        if ((Jpsi1)->mother()->mother()->particleID().pid()==(100443))twoFeed = true; //psi(2s)
        else  if (Jpsi1->mother()->mother()->particleID().pid()==(30443))twoFeed=true; //psi(3770)
        else  if (Jpsi1->mother()->mother()->particleID().pid()==(9900441))twoFeed=true; //cc-bar[1S08]
        else if (Jpsi1->mother()->mother()->particleID().pid()==(9900443)) twoFeed=true; ////cc-bar[3S08]
        else if(Jpsi1->mother()->mother()->particleID().pid()==(9910441)) twoFeed = true;//cc-bar[3P08]
        else {
          //counter ("zzzzf. rejected twofeed")++;
          return 0.0;
        }
      }
      counter ("TruthH. after rejected twofeed")++;
      if (fromPV == false && oneFeed ==false && twoFeed == false){
        //counter ("zzzg. 2nd rejected non-prompt")++;
        return 0.0; 
      }

        counter ("TruthI. after rejected non-prompt")++;
      if(Jpsi1->particleID().pid()!= 443){
        // std::cout << "pid of Jpsi not a Jpsi: "  << jpsi->particleID().pid()  << std::endl;
        //counter("zzzh. Jpsi not a Jpsi")++;
        return 0.0;
      }

   

      double Jpsi1pid = Jpsi1->particleID().pid();
      double Jpsi2pid = Jpsi2->particleID().pid();
     
      counter("TruthJ. Jpsi not a Jpsi")++;
      
      if(Jpsi1->key() != Jpsi2->key() ) return 0.0;
      
      counter("TruthK. last")++;
      
      mothercode = Jpsi1pid;     
    }
  }
  
  return mothercode; 
}




double DoubleJpsiSel::DoubleTrackToHEPMCOld(const LHCb::Track* part_trkP,const LHCb::Track* part_trkN)const{
  bool mclookup = true;
  bool fromjpsi = false;
  const IHepMC2MC* m_mc2hepmc;
  LHCb::Track2MC2D::Range mclinksN,mclinksP;//changed to const
  const IHepMC2MC::MC2HepMC* mctable;
  double mothercode = 0.0;
  const LHCb::Track2MC2D* table2d;
  counter("start of DoubleTrackToHepMC")++;
  if (exist<LHCb::Track2MC2D>(LHCb::Track2MCLocation::Default)){
    table2d = get<LHCb::Track2MC2D>(LHCb::Track2MCLocation::Default);}
  else{
    std::cout << "::DoubleTrackToHEPMC - no HepMCEvents" << std::endl;
    return 0.0;    
  }



  //  table2d = get<LHCb::Track2MC2D>(LHCb::Track2MCLocation::Default);
  m_mc2hepmc = tool<IHepMC2MC >("LoKi::HepMC2MC");
  if( m_mc2hepmc->mc2HepMC()!= NULL)mctable = m_mc2hepmc->mc2HepMC();
  mclinksP = table2d->relations(part_trkP);
  mclinksN = table2d->relations(part_trkN);
  //  std::cout << "start of HEP LOOP" <<std::endl;

  if (0 == mctable){ //check to see if there is data=============
    mclookup = false;
    std::cout << "::DoubleTrackToHEPMC -  MC2HepMC problematic " << std::endl;
  }
  for( LHCb::Track2MC2D::Range::iterator mclinkP=mclinksP.begin();mclinksP.end()!=mclinkP;++mclinkP){
    for( LHCb::Track2MC2D::Range::iterator mclinkN=mclinksN.begin();mclinksN.end()!=mclinkN;++mclinkN){
      const LHCb::MCParticle* mcpartP = mclinkP->to();
      const LHCb::MCParticle* mcpartN = mclinkN->to();
      if(mcpartP->key() == mcpartN->key())return 0.0; ///check to see if same muon
      if (mclookup){
        const IHepMC2MC::MC2HepMC::Range mclinksP = mctable->relations(mcpartP);
        const IHepMC2MC::MC2HepMC::Range mclinksN = mctable->relations(mcpartN);
        if((!mclinksN.empty())&&(!mclinksP.empty())){
          //generating particle and production vertex - =================
          const HepMC::GenParticle* irP = mclinksP.back().to();
          const HepMC::GenParticle* irN = mclinksN.back().to();
          HepMC::GenVertex* pvP = irP->production_vertex();
          HepMC::GenVertex* pvN = irN->production_vertex();
          //std::cout << "before check on pdg of muon "<<std::endl;
          // std::cout << "PID of muon in DoubleTrack"<<irN->pdg_id() <<std::endl;
          if ((irN->pdg_id())!=-13)return 0.0;//=====check to see if real muon of correct sign==========
          if ((irP->pdg_id())!=13)return 0.0;
          if ((pvN == NULL)&&(pvP == NULL))return 0.0;
          for(HepMC::GenVertex::particles_in_const_iterator ivP=pvP->particles_in_const_begin(); \
              pvP->particles_in_const_end()!=ivP;++ivP){
            for(HepMC::GenVertex::particles_in_const_iterator ivN=pvN->particles_in_const_begin(); \
                pvN->particles_in_const_end()!=ivN;++ivN){
                    
              const HepMC::GenParticle* mc2partP = *ivP; 
              const HepMC::GenParticle* mc2partN = *ivN;
              if(mc2partP->barcode() != mc2partN->barcode())return 0.0; //check to see both from same particle
              int pdg_id2 = mc2partP->pdg_id();
              pdg_id2 = abs(pdg_id2);
              //std::cout <<"before 443 cut" <<std::endl;
              //   std::cout << "PID of Jpsi in DoubleTrack"<< pdg_id2  <<std::endl;
              if(pdg_id2 != 443)return 0.0; //=====check to see if pv particle (of muon) is j/psi========
              counter("I. after J/psi truth ")++;

              //=================
              //  if(m_NoAccCutMC){
              // bool inAcceptance = AccCheckMCTruth(mc2partP);
              //  if(!inAcceptance) continue;
              // }
              

              //=================

              //==================check to see if comes from c-c(bar)=====
              if ( mc2partP->production_vertex() ) {
                for ( HepMC::GenVertex::particle_iterator mother = mc2partP->production_vertex()->
                        particles_begin(HepMC::parents);
                      mother != mc2partP->production_vertex()->particles_end(HepMC::parents); ++mother )
                {
   
                  // const IHepMC2MC::HepMC2MC::Range mclinksjpsi = mctable->relations(mc2partP);
                  //const LHCb::MCParticle* mcpartjpsi = mclinksjpsi.back().to();

                  

                  // const LHCb::MCVertex * jpsiOrgVertex = mcpartjpsi->originVertex();
                  //bool fromPV = false;
                  //bool fromFeed = false;
                  //if (jpsiOrgVertex->isPrimary()) fromPV = true;
                  //else if((jpsiOrgVertex->mother()->originVertex())==0){
                  //  if ((jpsiOrgVertex->mother()->originVertex())->isPrimary()) fromFeed = true;
                  // }
                  //if (fromPV) prompt = true;
                  //else if ((*mother)->pdg_id()==(100443) && fromFeed) prompt=true; //psi(2s)
                  //else if ((*mother)->pdg_id()==(445) && fromFeed) prompt=true; //Chic2(1P)
                  //else if ((*mother)->pdg_id()==(100445) && fromFeed) prompt=true; //chic2(2P)
                  //else prompt = false;


                  ///// if ((*mother)->pdg_id()==(9900441))fromjpsi=true;
                  /////else if ((*mother)->pdg_id()==(9900443)) fromjpsi=true;
                  /////else if((*mother)->pdg_id()==(9910441)) fromjpsi = true;
                  /////else if((*mother)->pdg_id()==(91)) fromjpsi = true;
                  /////else fromjpsi = false;

                  // double rap =  mcpartjpsi->momentum().Rapidity();
                  //if (!((rap >= 2.0) && (rap <= 4.5))) continue;


                  fromjpsi = true;
                  mothercode =(*mother)->pdg_id();
                  if (fromjpsi) {
                    //counter("G. after J/psi TRUTH prompt")++;
                    Tuple tuple = nTuple("HEPMCevttuple");
                    tuple->column("HEPMCGenmass",mc2partP->generated_mass());  
                    tuple->column("HEPMCmass",mc2partP->momentum().m());   
                    tuple->write();
                  }
                } 
              }
            }
          }      
        }
      }
    }
  }
  return mothercode; 
}



//===========================================================================
//check on MCTruth to see if in acceptance (to account for generator level MC)
//===========================================================================

bool DoubleJpsiSel::AccCheckMCTruth(const HepMC::GenParticle* theSignal)const
{
  counter("start of AccCheckMCTruth")++;
HepMC::GenVertex * EV = theSignal -> end_vertex() ;
  if ( 0 == EV ) return true ;

  typedef std::vector< HepMC::GenParticle * > Particles ;
  Particles stables ;
  HepMC::GenVertex::particle_iterator iter ;

  for ( iter = EV -> particles_begin( HepMC::descendants ) ; 
        iter != EV -> particles_end( HepMC::descendants ) ; ++iter ) {
    if ( 0 == (*iter) -> end_vertex() ) stables.push_back( *iter ) ;
  }  

  if ( stables.empty() )
    Exception( "Signal has no stable daughters !" ) ;

  double angle( 0. ) ;
  double firstpz = stables.front() -> momentum().pz() ;

  debug() << "New event" << endmsg ;

  for ( Particles::const_iterator it = stables.begin() ; it != stables.end() ;
        ++it ) {

    debug() << "Check particle " << (*it) -> pdg_id() << " with angle " 
            << (*it) -> momentum().theta() / Gaudi::Units::mrad 
            << " mrad." << endmsg ;
   
    // Remove neutrinos
    if ( ( 12 == abs( (*it) -> pdg_id() ) ) || 
         ( 14 == abs( (*it) -> pdg_id() ) ) || 
         ( 16 == abs( (*it) -> pdg_id() ) ) ) continue ;
 
    // Don't use daughters of Lambda and KS:
    HepMC::GenParticle * theParent ;
    theParent = 
      *( (*it) -> production_vertex() -> particles_in_const_begin() ) ;
    if ( 3122 == abs( theParent -> pdg_id() ) ) continue ;
    if ( 310 == theParent -> pdg_id() ) continue ;

    // Consider only gammas from pi0 and eta
    if ( 22 == (*it) -> pdg_id() ) {
      if ( ( 111 != theParent -> pdg_id() ) &&
           ( 221 != theParent -> pdg_id() ) ) continue ;
    }

    // All particles in same direction
    if ( 0 > ( firstpz * ( (*it) -> momentum().pz() ) ) ) return false ;
      
    angle = (*it) -> momentum().theta() ;
     counter("point 2 of AccCheckMCTruth")++;
    LHCb::ParticleID pid( (*it) -> pdg_id() ) ;
    if ( 0 == pid.threeCharge() ) {
      if ( fabs( sin( angle ) ) > fabs( sin( m_neutralThetaMax ) ) ) 
        return false ;
      if ( fabs( sin( angle ) ) < fabs( sin( m_neutralThetaMin ) ) ) 
        return false ;
    } else {
      if ( fabs( sin( angle ) ) > fabs( sin( m_chargedThetaMax ) ) ) 
        return false ;
      if ( fabs( sin( angle ) ) < fabs( sin( m_chargedThetaMin ) ) ) 
        return false ;
    }
  }

  counter("point 3 of AccCheckMCTruth")++;
  debug() << "Event passed !" << endmsg ;
  
  return true ;
}


