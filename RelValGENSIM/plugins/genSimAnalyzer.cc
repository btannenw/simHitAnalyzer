// -*- C++ -*-
//
// Package:    L1Trigger/genSimAnalyzer
// Class:      genSimAnalyzer
//
/**\class genSimAnalyzer genSimAnalyzer.cc L1Trigger/genSimAnalyzer/plugins/genSimAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Benjamin Tannenwald
//         Created:  Tue, 20 Nov 2018 16:40:40 GMT
//
//


// system include files
#include <memory>
#include <utility>

// Math Include
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// Includes for SimHits
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

// Includes for SimTracks
#include "SimDataFormats/Track/interface/SimTrack.h"

// includes for global positioning
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h" 
#include "DataFormats/DetId/interface/DetId.h" 
#include "FWCore/Framework/interface/ESWatcher.h"


//#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h" 
//#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h" 

// Includes for output file
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class genSimAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit genSimAnalyzer(const edm::ParameterSet&);
      ~genSimAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT< std::vector<PSimHit> > fastTimeHitsBarrelToken_; // BBT 11-20-18
      edm::EDGetTokenT< std::vector<PSimHit> > fastTimeHitsEndcapToken_; // BBT 11-20-18
      edm::EDGetTokenT< std::vector<reco::GenParticle> > genParticlesToken_; // BBT 11-20-18
      edm::EDGetTokenT< std::vector<SimTrack> > simTracksToken_; // BBT 11-20-18
      edm::EDGetTokenT< float > genParticlesTimeToken_;
      edm::InputTag fastTimeBarrelHitsTag; // BBT 11-20-18
      edm::InputTag fastTimeEndcapHitsTag; // BBT 11-20-18
      edm::InputTag genParticlesTag; // BBT 11-20-18
      edm::InputTag genParticlesTimeTag; // BBT 11-20-18
      edm::InputTag simTracksTag; // BBT 11-20-18
      int pdgIdTest;
  /*edm::InputTag primarySimVertexPosTag; // BBT 10-26-18
      edm::InputTag primarySimVertexTimeTag; // BBT 10-26-18
      edm::EDGetTokenT< ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> > primarySimVertexPosToken_;
      edm::EDGetTokenT< float > primarySimVertexTimeToken_;
  */
  // --------------- output tree and members ---------------
      TTree* simHitTree;
  
      std::vector<double> v_barrel_E_loss; 
      std::vector<double> v_barrel_moduleType; 
      std::vector<double> v_barrel_particleType; 
      std::vector<double> v_barrel_DetUnitID; 
      std::vector<double> v_barrel_trackID; 
      std::vector<double> v_barrel_absP; 
      std::vector<double> v_barrel_localX; 
      std::vector<double> v_barrel_localY; 
      std::vector<double> v_barrel_localZ; 

      std::vector<double> v_endcap_E_loss; 
      std::vector<double> v_endcap_particleType; 
      std::vector<double> v_endcap_DetUnitID; 
      std::vector<double> v_endcap_trackID; 
      std::vector<double> v_endcap_localX; 
      std::vector<double> v_endcap_localY; 
      std::vector<double> v_endcap_localZ; 
      std::vector<double> v_endcap_globalX; 
      std::vector<double> v_endcap_globalY; 
      std::vector<double> v_endcap_globalZ; 
      std::vector<double> v_endcap_globalPt; 
      std::vector<double> v_endcap_globalEta; 
      std::vector<double> v_endcap_globalPhi; 

      std::vector<double> v_endcap_absP; 
      std::vector<double> v_endcap_tof;
      std::vector<double> v_endcap_pt;  
      std::vector<double> v_endcap_localX_entry; 
      std::vector<double> v_endcap_localY_entry; 
      std::vector<double> v_endcap_localZ_entry; 
      std::vector<double> v_endcap_localX_exit; 
      std::vector<double> v_endcap_localY_exit; 
      std::vector<double> v_endcap_localZ_exit; 

      std::vector<double> v_simTrack_X; 
      std::vector<double> v_simTrack_Y; 
      std::vector<double> v_simTrack_Z; 
      std::vector<double> v_simTrack_pt; 
      std::vector<double> v_simTrack_eta; 
      std::vector<double> v_simTrack_phi; 
      std::vector<double> v_simTrack_trackID; 
      std::vector<double> v_simTrack_pdgID; 

  /*
      // Vertex 10-26-18 BBT
      double d_vertexX; 
      double d_vertexY; 
      double d_vertexZ; 
      double d_vertexTime; 

      // SimHit 10-26-18 BBT
      std::vector<double> v_simHitX; 
      std::vector<double> v_simHitY; 
      std::vector<double> v_simHitZ; 
      std::vector<double> v_simHitTime; 
  */

      // --------------- histograms ---------------
      TH1F* track_pt;
      TH1F* barrel_E_loss;
      TH1F* barrel_localX;
      TH1F* barrel_localY;
      TH1F* barrel_localZ;
      TH1F* barrel_localX_modType1;
      TH1F* barrel_localX_modType2;
      TH1F* barrel_localX_modType3;
      TH1F* barrel_localY_modType1;
      TH1F* barrel_localY_modType2;
      TH1F* barrel_localY_modType3;
      TH1F* barrel_localZ_modType1;
      TH1F* barrel_localZ_modType2;
      TH1F* barrel_localZ_modType3;
      TH1F* endcap_E_loss;
      TH1F* endcap_localX;
      TH1F* endcap_localY;
      TH1F* endcap_localZ;


     edm::ESWatcher<MTDDigiGeometryRecord> geomwatcher_;
     edm::ESWatcher<MTDGeometryRecord> geomwatcher2_;
     const MTDGeometry* geom_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
genSimAnalyzer::genSimAnalyzer(const edm::ParameterSet& cfg)//:
  //fastTimeHitsBarrelToken_( consumes< std::vector<PSimHit> >(cfg.getParameter("fastTimeBarrelHits"))),
  //fastTimeHitsEndcapToken_( consumes< std::vector<PSimHit> >(cfg.getParameter("fastTimeEndcapHits"))),
  //genParticlesToken_( consumes< std::vector< reco::GenParticle> >(cfg.getParameter("genParticles_xyz"))),
  //genParticlesTimeToken_( consumes< float >(cfg.getParameter("genParticles_t")))
  //)
{
  // *** 0. Tags
  fastTimeBarrelHitsTag = cfg.getParameter<edm::InputTag>("fastTimeBarrelHits"); // BBT 11-20-18
  fastTimeEndcapHitsTag = cfg.getParameter<edm::InputTag>("fastTimeEndcapHits");; // BBT 11-20-18
  simTracksTag = cfg.getParameter<edm::InputTag>("simTracks");; // BBT 11-20-18
  genParticlesTag = cfg.getParameter<edm::InputTag>("genParticles"); // BBT 11-20-18
  genParticlesTimeTag = cfg.getParameter<edm::InputTag>("genParticles_t"); // BBT 11-20-18
  pdgIdTest    = cfg.getParameter<int>("pdgIdTest"); // BBT 11-20-18
  
  fastTimeHitsBarrelToken_ = consumes< std::vector<PSimHit> >(fastTimeBarrelHitsTag);
  fastTimeHitsEndcapToken_ = consumes< std::vector<PSimHit> >(fastTimeEndcapHitsTag);
  simTracksToken_ = consumes< std::vector<SimTrack> >(simTracksTag);
  //genParticlesTimeToken_   = consumes< float >(genParticlesTimeTag);
  //genParticlesToken_       = consumes< std::vector<reco::GenParticle> >(genParticlesTag);

   // output tree and members
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  simHitTree = fs->make<TTree>("simHitTree", "GEN-SIM Sim Hit Information");
  simHitTree->Branch("barrel_E_loss",    &v_barrel_E_loss);
  simHitTree->Branch("barrel_moduleType",    &v_barrel_moduleType);
  simHitTree->Branch("barrel_DetUnitID", &v_barrel_DetUnitID);
  simHitTree->Branch("barrel_trackID", &v_barrel_trackID);
  simHitTree->Branch("barrel_absP", &v_barrel_absP);
  simHitTree->Branch("barrel_pdgID", &v_barrel_particleType);
  simHitTree->Branch("barrel_localX",         &v_barrel_localX);
  simHitTree->Branch("barrel_localY",         &v_barrel_localY);
  simHitTree->Branch("barrel_localZ",         &v_barrel_localZ);

  simHitTree->Branch("endcap_E_loss",    &v_endcap_E_loss);
  simHitTree->Branch("endcap_DetUnitID", &v_endcap_DetUnitID);
  simHitTree->Branch("endcap_trackID", &v_endcap_trackID);
  simHitTree->Branch("endcap_pdgID", &v_endcap_particleType);
  simHitTree->Branch("endcap_localX",         &v_endcap_localX);
  simHitTree->Branch("endcap_localY",         &v_endcap_localY);
  simHitTree->Branch("endcap_localZ",         &v_endcap_localZ);
  simHitTree->Branch("endcap_globalX",         &v_endcap_globalX);
  simHitTree->Branch("endcap_globalY",         &v_endcap_globalY);
  simHitTree->Branch("endcap_globalZ",         &v_endcap_globalZ);
  simHitTree->Branch("endcap_globalPt",         &v_endcap_globalPt);
  simHitTree->Branch("endcap_globalEta",         &v_endcap_globalEta);
  simHitTree->Branch("endcap_globalPhi",         &v_endcap_globalPhi);

  simHitTree->Branch("endcap_absP",    &v_endcap_absP);
  simHitTree->Branch("endcap_tof",    &v_endcap_tof);
  simHitTree->Branch("endcap_pt",    &v_endcap_pt);
  simHitTree->Branch("endcap_localX_entry",         &v_endcap_localX_entry);
  simHitTree->Branch("endcap_localY_entry",         &v_endcap_localY_entry);
  simHitTree->Branch("endcap_localZ_entry",         &v_endcap_localZ_entry);
  simHitTree->Branch("endcap_localX_exit",         &v_endcap_localX_exit);
  simHitTree->Branch("endcap_localY_exit",         &v_endcap_localY_exit);
  simHitTree->Branch("endcap_localZ_exit",         &v_endcap_localZ_exit);

  simHitTree->Branch("simTrack_trackID",         &v_simTrack_trackID);
  simHitTree->Branch("simTrack_pdgID",         &v_simTrack_pdgID);
  simHitTree->Branch("simTrack_X",         &v_simTrack_X);
  simHitTree->Branch("simTrack_Y",         &v_simTrack_Y);
  simHitTree->Branch("simTrack_Z",         &v_simTrack_Z);
  simHitTree->Branch("simTrack_pt",         &v_simTrack_pt);
  simHitTree->Branch("simTrack_eta",         &v_simTrack_eta);
  simHitTree->Branch("simTrack_phi",         &v_simTrack_phi);

  /*simHitTree->Branch("track_eta",    &v_trackEta);
  simHitTree->Branch("track_time",    &v_trackTime);
  simHitTree->Branch("track_time_smeared",    &v_trackTime_smeared);
  simHitTree->Branch("track_time_PVcorrected",    &v_trackTime_PVcorrected);

  simHitTree->Branch("vertex_X",    &d_vertexX);
  simHitTree->Branch("vertex_Y",    &d_vertexY);
  simHitTree->Branch("vertex_Z",    &d_vertexZ);
  simHitTree->Branch("vertex_time", &d_vertexTime);
  */
  // histograms
  track_pt             = fs->make<TH1F>( "track_pt"  , "Sim Track p_{t}", 42,  0., 10.5 );
  barrel_E_loss        = fs->make<TH1F>( "barrel_E_loss"  , "E_loss (Barrel)", 48,  0., 12. );
  barrel_localX        = fs->make<TH1F>( "barrel_localX"  , "Module X-Coordinate (Barrel)", 30,  -1.5, 1.5 );
  barrel_localY        = fs->make<TH1F>( "barrel_localY"  , "Module Y-Coordinate (Barrel)", 30,  -30, 30 );
  barrel_localZ        = fs->make<TH1F>( "barrel_localZ"  , "Module Z-Coordinate (Barrel)", 30,  -2, 2 );
  endcap_E_loss        = fs->make<TH1F>( "endcap_E_loss"  , "E_loss (Endcap)", 48,  0., 0.4 );
  endcap_localX        = fs->make<TH1F>( "endcap_localX"  , "Module X-Coordinate (Endcap)", 30,  -25, 25 );
  endcap_localY        = fs->make<TH1F>( "endcap_localY"  , "Module Y-Coordinate (Endcap)", 30,  -50, 50 );
  endcap_localZ        = fs->make<TH1F>( "endcap_localZ"  , "Module Z-Coordinate (Endcap)", 30,  -0.15, 0.15 );

  barrel_localX_modType1 = fs->make<TH1F>( "barrel_localX_modType1"  , "Module X-Coordinate (Barrel)", 30,  -1.5, 1.5 );
  barrel_localY_modType1 = fs->make<TH1F>( "barrel_localY_modType1"  , "Module Y-Coordinate (Barrel)", 30,  -30, 30 );
  barrel_localZ_modType1 = fs->make<TH1F>( "barrel_localZ_modType1"  , "Module Z-Coordinate (Barrel)", 30,  -2, 2 );

  barrel_localX_modType2 = fs->make<TH1F>( "barrel_localX_modType2"  , "Module X-Coordinate (Barrel)", 30,  -1.5, 1.5 );
  barrel_localY_modType2 = fs->make<TH1F>( "barrel_localY_modType2"  , "Module Y-Coordinate (Barrel)", 30,  -30, 30 );
  barrel_localZ_modType2 = fs->make<TH1F>( "barrel_localZ_modType2"  , "Module Z-Coordinate (Barrel)", 30,  -2, 2 );

  barrel_localX_modType3 = fs->make<TH1F>( "barrel_localX_modType3"  , "Module X-Coordinate (Barrel)", 30,  -1.5, 1.5 );
  barrel_localY_modType3 = fs->make<TH1F>( "barrel_localY_modType3"  , "Module Y-Coordinate (Barrel)", 30,  -30, 30 );
  barrel_localZ_modType3 = fs->make<TH1F>( "barrel_localZ_modType3"  , "Module Z-Coordinate (Barrel)", 30,  -2, 2 );

}
		  

genSimAnalyzer::~genSimAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
genSimAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //std::cout << " !!! Running genSimAnalyzer::analyze() !!!" << std::endl;

   using namespace edm;
   
   // Rel-Val testing for GEN-SIM
   // *** 1. Gen particles
   //edm::Handle< float > genParticlesTimeHandle;
   //iEvent.getByToken(genParticlesTimeToken_, genParticlesTimeHandle);
   //const auto& genParticlesTime = *genParticlesTimeHandle;
   //edm::Handle< std::vector< reco::GenParticle> > genParticlesHandle;
   //iEvent.getByToken(genParticlesToken_, genParticlesHandle);
   //const auto& genParticles = *genParticlesHandle;

   // *** 2. MTD Sim Hits
   edm::Handle< std::vector< PSimHit > > fastTimeHitsBarrelHandle;
   iEvent.getByToken(fastTimeHitsBarrelToken_, fastTimeHitsBarrelHandle);
   const auto& fastTimeHitsBarrel = *fastTimeHitsBarrelHandle;
   edm::Handle< std::vector< PSimHit > > fastTimeHitsEndcapHandle;
   iEvent.getByToken(fastTimeHitsEndcapToken_, fastTimeHitsEndcapHandle);
   const auto& fastTimeHitsEndcap = *fastTimeHitsEndcapHandle;
   /*
   // access the MTD
   std::cout << "1" << std::endl;
   edm::ESHandle<MTDGeometry> theMTDGeometry;
   std::cout << "1" << std::endl;
   iSetup.get<MTDDigiGeometryRecord>().get(theMTDGeometry);
   std::cout << "1" << std::endl;
   //iSetup.get<MTDGeometryRecord>().get(theMTDGeometry);
   std::cout << "1" << std::endl;
   const MTDGeometry& theMTD(*theMTDGeometry);
   std::cout << "1" << std::endl;
   */
   /*
   edm::ESHandle<MTDGeometry> geom;
   std::cout << "1" << std::endl;
   if( geomwatcher_.check(iSetup) || geom_ == nullptr ) {
     iSetup.get<MTDDigiGeometryRecord>().get(geom);
     std::cout << "1" << std::endl;
     geom_ = geom.product();
     std::cout << "1" << std::endl;
     const MTDGeometry& theMTD(*geom);
   }
   */
   // *** 3. Sim Tracks
   edm::Handle< std::vector< SimTrack > > simTracksHandle;
   iEvent.getByToken(simTracksToken_, simTracksHandle);
   const auto& simTracks = *simTracksHandle;

   /*
   // Timing
   edm::Handle<edm::ValueMap<float> > timingValues;
   edm::Handle<edm::ValueMap<float> > timingValues_smeared;
   iEvent.getByToken(timingValuesToken_,timingValues);
   iEvent.getByToken(timingValuesSmearedToken_,timingValues_smeared);

   // pair stuff
   std::vector< std::pair < TTTrack< Ref_Phase2TrackerDigi_ >, double > > l1Tracks_pair;

   // Get vertex info, BBT 10-26-18
   edm::Handle< ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> > primarySimVertexPos;
   iEvent.getByToken(primarySimVertexPosToken_, primarySimVertexPos);
   edm::Handle< float > primarySimVertexTime;
   iEvent.getByToken(primarySimVertexTimeToken_, primarySimVertexTime);
   */
   // -----------------------------------------------------------------------------------------------

   // variables for output
   v_barrel_E_loss.clear();
   v_barrel_moduleType.clear();
   v_barrel_DetUnitID.clear(); 
   v_barrel_trackID.clear(); 
   v_barrel_absP.clear(); 
   v_barrel_particleType.clear(); 
   v_barrel_localX.clear(); 
   v_barrel_localY.clear(); 
   v_barrel_localZ.clear(); 
   
   v_endcap_E_loss.clear(); 
   v_endcap_DetUnitID.clear(); 
   v_endcap_trackID.clear(); 
   v_endcap_particleType.clear(); 
   v_endcap_localX.clear(); 
   v_endcap_localY.clear(); 
   v_endcap_localZ.clear(); 
   v_endcap_globalX.clear(); 
   v_endcap_globalY.clear(); 
   v_endcap_globalZ.clear(); 
   v_endcap_globalPt.clear(); 
   v_endcap_globalEta.clear(); 
   v_endcap_globalPhi.clear(); 

   v_endcap_absP.clear(); 
   v_endcap_pt.clear(); 
   v_endcap_tof.clear(); 
   v_endcap_localX_entry.clear(); 
   v_endcap_localY_entry.clear(); 
   v_endcap_localZ_entry.clear(); 
   v_endcap_localX_exit.clear(); 
   v_endcap_localY_exit.clear(); 
   v_endcap_localZ_exit.clear(); 

   v_simTrack_X.clear(); 
   v_simTrack_Y.clear(); 
   v_simTrack_Z.clear(); 
   v_simTrack_pt.clear(); 
   v_simTrack_eta.clear(); 
   v_simTrack_phi.clear(); 
   v_simTrack_trackID.clear(); 
   v_simTrack_pdgID.clear(); 
   

   /*
   //l1Tracks_time.clear();
   v_trackPt.clear();
   v_trackEta.clear();
   v_trackTime.clear();
   v_trackTime_smeared.clear();
   v_trackTime_PVcorrected.clear();
   d_vertexTime = -999;
   d_vertexX = -999;
   d_vertexY = -999;
   d_vertexZ = -999;
   */
/*
   // VERTEX stuff, BBT 10-26-18
   d_vertexTime = (*primarySimVertexTime) ;
   d_vertexX    = (*primarySimVertexPos).X() ;
   d_vertexY    = (*primarySimVertexPos).Y() ;
   d_vertexZ    = (*primarySimVertexPos).Z() ;
   */

   // ----------------------------------------------------------------------------------------------
   // loop over simTracks
   // ----------------------------------------------------------------------------------------------
    for( const auto& simTrack : simTracks ) {  
     //std::cout << "NEW TRACK" << std::endl;
      //if ( abs(simTrack.type()) == pdgIdTest) { // only proceed if particle has the pdgID of the particle we're testing
      if ( simTrack.trackId() == 1 || simTrack.trackId() == 2 ) { // only proceed if particle has the pdgID of the particle we're testing
       v_simTrack_X.push_back( simTrack.momentum().x() ); 
       v_simTrack_Y.push_back( simTrack.momentum().y() ); 
       v_simTrack_Z.push_back( simTrack.momentum().z() ); 
       v_simTrack_pt.push_back( simTrack.momentum().pt() ); 
       v_simTrack_eta.push_back( simTrack.momentum().eta() ); 
       v_simTrack_phi.push_back( simTrack.momentum().phi() ); 
       v_simTrack_trackID.push_back( simTrack.trackId() ); 
       v_simTrack_pdgID.push_back( simTrack.type() ); 
       track_pt->Fill( simTrack.momentum().pt() );
     }
   }

   // ----------------------------------------------------------------------------------------------
   // loop over barrel timing hits
   // ----------------------------------------------------------------------------------------------
   for( const auto& barrelHit : fastTimeHitsBarrel ) {  
     //std::cout << "NEW TRACK" << std::endl;
     //if ( abs(barrelHit.particleType()) == pdgIdTest) { // only proceed if particle has the pdgID of the particle we're testing
     if ( barrelHit.trackId() == 1 || barrelHit.trackId() == 2 ) { // only proceed if particle has the pdgID of the particle we're testing
       v_barrel_E_loss.push_back( 1e3*barrelHit.energyLoss() ); // GeV -> MeV
       v_barrel_DetUnitID.push_back( barrelHit.detUnitId() );
       v_barrel_particleType.push_back( barrelHit.particleType() );
       v_barrel_localX.push_back( barrelHit.localPosition().x() );
       v_barrel_localY.push_back( barrelHit.localPosition().y() );
       v_barrel_localZ.push_back( barrelHit.localPosition().z() );
       v_barrel_trackID.push_back( barrelHit.trackId() );
       v_barrel_absP.push_back( barrelHit.pabs() );

       BTLDetId hitId( barrelHit.detUnitId() );
       //DetId geoId = BTLDetId(hitId.mtdSide(),hitId.mtdRR(),hitId.module()+18*(hitId.modType()-1),0,1);   
       v_barrel_moduleType.push_back( hitId.modType() );

       barrel_E_loss->Fill(1e3* barrelHit.energyLoss() );
       barrel_localX->Fill( barrelHit.localPosition().x() );
       barrel_localY->Fill( barrelHit.localPosition().y() );
       barrel_localZ->Fill( barrelHit.localPosition().z() );
       if (hitId.modType() == 1) {
	 barrel_localX_modType1->Fill( barrelHit.localPosition().x() );
	 barrel_localY_modType1->Fill( barrelHit.localPosition().y() );
	 barrel_localZ_modType1->Fill( barrelHit.localPosition().z() );
       }
       if (hitId.modType() == 2) {
	 barrel_localX_modType2->Fill( barrelHit.localPosition().x() );
	 barrel_localY_modType2->Fill( barrelHit.localPosition().y() );
	 barrel_localZ_modType2->Fill( barrelHit.localPosition().z() );
       }
       if (hitId.modType() == 3) {
	 barrel_localX_modType3->Fill( barrelHit.localPosition().x() );
	 barrel_localY_modType3->Fill( barrelHit.localPosition().y() );
	 barrel_localZ_modType3->Fill( barrelHit.localPosition().z() );
       }

     }
   }


   // ----------------------------------------------------------------------------------------------
   // loop over endcap timing hits
   // ----------------------------------------------------------------------------------------------
   for( const auto& endcapHit : fastTimeHitsEndcap ) {  
     //std::cout << "NEW TRACK" << std::endl;
     if ( endcapHit.trackId() == 1 || endcapHit.trackId() == 2 ) { // only proceed if particle has the pdgID of the particle we're testing
       //if ( abs(endcapHit.particleType()) == pdgIdTest) { // only proceed if particle has the pdgID of the particle we're testing
       v_endcap_E_loss.push_back( 1e3*endcapHit.energyLoss() ); // GeV -> MeV
       v_endcap_DetUnitID.push_back( endcapHit.detUnitId() );
       v_endcap_particleType.push_back( endcapHit.particleType() );
       v_endcap_localX.push_back( endcapHit.localPosition().x() );
       v_endcap_localY.push_back( endcapHit.localPosition().y() );
       v_endcap_localZ.push_back( endcapHit.localPosition().z() );

       v_endcap_localX_entry.push_back( endcapHit.entryPoint().x() );
       v_endcap_localY_entry.push_back( endcapHit.entryPoint().y() );
       v_endcap_localZ_entry.push_back( endcapHit.entryPoint().z() );
       v_endcap_localX_exit.push_back( endcapHit.exitPoint().x() );
       v_endcap_localY_exit.push_back( endcapHit.exitPoint().y() );
       v_endcap_localZ_exit.push_back( endcapHit.exitPoint().z() );
       v_endcap_tof.push_back( endcapHit.timeOfFlight() );
       v_endcap_absP.push_back( endcapHit.pabs() );
       v_endcap_pt.push_back( endcapHit.momentumAtEntry().perp() );


       /*
       std::cout << "2" << std::endl;
       ETLDetId hitId( endcapHit.detUnitId() );
       std::cout << "2" << std::endl;
       DetId geoId = ETLDetId(hitId.mtdSide(),hitId.mtdRR(),hitId.module(),0);
       //ROOT::Math::XYZVector endcapGlobalPos(geoId.position().x(),geoId.position().y(),geoId.position().z());
       std::cout << "2" << std::endl;
       const GeomDet *theDet = theMTD.idToDet(geoId);
       std::cout << "2" << std::endl;
       */
       // create a DetId from the detUnitId
       //DetId theDetUnitId( endcapHit.detUnitId());
       // get the DetUnit via the DetUnitId and cast it to a StripGeomDetUnit
       //const GeomDet *theDet = theMTD.idToDet(theDetUnitId);

       v_endcap_globalX.push_back( endcapHit.momentumAtEntry().x() );
       v_endcap_globalY.push_back( endcapHit.momentumAtEntry().y() );
       //v_endcap_globalY.push_back( theDet->surface().toGlobal(endcapHit.localPosition()).x() );
       //v_endcap_globalY.push_back( geoId.surface().toGlobal(endcapHit.localPosition()).x() );
       v_endcap_globalZ.push_back( endcapHit.momentumAtEntry().z() );
       v_endcap_globalPt.push_back( endcapHit.momentumAtEntry().perp() );
       v_endcap_globalEta.push_back( endcapHit.momentumAtEntry().eta() );
       v_endcap_globalPhi.push_back( endcapHit.momentumAtEntry().phi() );

       v_endcap_trackID.push_back( endcapHit.trackId() );

       endcap_E_loss->Fill( 1e3*endcapHit.energyLoss() );
       endcap_localX->Fill( endcapHit.localPosition().x() ); // phi
       endcap_localY->Fill( endcapHit.localPosition().y() ); // radial
       endcap_localZ->Fill( endcapHit.localPosition().z() ); // z --> should be thin

     }
   }
   
   // ***. Fill tree
   simHitTree->Fill();


   
}


// ------------ method called once each job just before starting event loop  ------------
void
genSimAnalyzer::beginJob()
{
  std::cout << "STARTING JOB!!!" << std::endl;

  //m_allstub_x = new std::vector<float>;
  //m_allstub_y = new std::vector<float>;
  //m_allstub_z = new std::vector<float>;

}

// ------------ method called once each job just after ending the event loop  ------------
void
genSimAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
genSimAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(genSimAnalyzer);
