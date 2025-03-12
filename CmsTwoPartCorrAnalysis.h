// system include files
#include <memory>

// CMSSW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "Analyzers/CmsTwoPartCorrAnalysis/data/EFF/trackingEfficiency2018PbPb.h"
// user include files
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TVector3.h"
#include <string>
#include "Analyzers/CmsTwoPartCorrAnalysis/interface/DiHadronCorrelationEvt.h"
//
// class declaration
//

class CmsTwoPartCorrAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit CmsTwoPartCorrAnalysis(const edm::ParameterSet&);
      ~CmsTwoPartCorrAnalysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      void   LoopVertices(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void   LoopTracksMC(const edm::Event& iEvent, const edm::EventSetup& iSetup,
			  bool istrg, int evtclass, bool isMC );
      void   LoopTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup,
                        bool istrg, int evtclass);
      double GetEffWeight(double eta, double pt, int evtclass);
      double getEff(const TLorentzVector a, int evtclass);
      void   AssignpTbins(double pt,  double eta,
                          double phi, int charge,
                          double eff1, double eff2, bool istrg, bool isGen);
            


      double deltaeta(const TVector3 a, const TVector3 b);
      double deltaphi(const TVector3 a, const TVector3 b);


      double GetDeltaEta(double eta_trg, double eta_ass);
      double GetDeltaPhi(double phi_trg, double phi_ass);
      double GetDPhiStar(double phi1, double pt1, int charge1, double phi2, double pt2, int charge2, double radius, double bSign);

    //   double Getqinv(double qinv , double Px_trg, double Px_ass,  double Py_trg, double Py_ass, double Pz_trg, double Pz_ass , double Energy, double P);

      void FillHistsSignal( int ievt, bool isGen );
      void FillHistsBackground(int ievt_trg, int jevt_ass, bool isGen );

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


      TVector3 trigger;
      TVector3 associated;

      // ----------member data ---------------------------
      // ## tracks ##
      // used to select what tracks to read from configuration file
      edm::EDGetTokenT<reco::TrackCollection> trackTags;
      edm::EDGetTokenT<reco::GenParticleCollection> ftrackMC;

      // ## vertex ##
      // used to select what vertex to read from configuration file
      edm::EDGetTokenT<reco::VertexCollection> vtxTags;

      // ## centrality ##
      // used to select what centrality collection to read from configuration file
      edm::EDGetTokenT<reco::Centrality> centralityTags;
      // used to access centrality bins
      edm::EDGetTokenT<int> centralityBinTags;
      
      //edm::Handle<std::vector<float>> mvaoutput_;
      //edm::InputTag mvaSrcLabel_;
      edm::EDGetTokenT<std::vector<float>> mvaSrc;
    
      
      // ## event classifier ##
      int evtclassifier;//Switch to centrality dependent correction instead of mult
      int centmin;
      int centmax;
      int noffmin;
      int noffmax;
      bool fIsMC; //for MC run---
      bool IsEffCorrection;
      //# Switch for QA histogram-------
      bool isQAhisto;
      // ## Eff/Fake coorection
      bool cweight; //Apply (True) or not (False) the eff/fake correction
      edm::InputTag fname; //Name of the file that contains the eff corrections                                                                                                                            
      edm::InputTag fhpt;
      edm::InputTag fmb; 
      edm::InputTag fpix; 
      edm::InputTag fplus; 
      edm::InputTag fminus; 
      std::vector<int> effCorrBinMin;
      std::vector<int> effCorrBinMax;
      TFile* feff; //ROOT File that contains histograms used for corrections
      std::vector<TH2D*> heff; //Histograms used for corrections (eta, pT)
      
      // ## Vertex variables
      double zminVtx; //min value for Z cut on vtx position
      double zmaxVtx; //max value for Z cut on vtx position
      double rhomaxVtx; //max value for XY cut on vtx position
      unsigned int nTrkAssoToVtx; //number of track requiered to be associated to a vertex
      bool selectVtxByMult;       //False: sel best vtx by sum pT^2 True: sel best vtx with highest multiplicity
      int nVtx;
      double xBestVtx; //Best vertex X position
      double yBestVtx; //Best vertex Y position
      double zBestVtx; //Best vertex Z position
      double rhoBestVtx; //Best vertex XY position
      double xBestVtxError; //Best vertex X error
      double yBestVtxError; //Best vertex Y error
      double zBestVtxError; //Best vertex Z error
      
      double pTmin_trg; //min pt of the trigger tracks
      double pTmax_trg; //max pt of the trigger tracks
      double pTmin_ass; //min pt of the associated tracks
      double pTmax_ass; //max pt of the associated tracks
      double etamin_trg; //min eta of the trigger tracks
      double etamax_trg; //max eta of the trigger tracks
      double etamin_ass; //min eta of the associated tracks
      double etamax_ass; //max eta of the associated tracks
      unsigned int bkgFactor;
      int nEtaBins;
      int nPhiBins;

      int nTrkTot_trgGen;
      int nTrkTot_assGen;
      int nTrkTot_trg;
      int nTrkTot_ass;
      double nTrkTot_trgCorrtd;
      double nTrkTot_assCorrtd;
      
      // ## Dihardon corr events ##
      DiHadronCorrelationEvt* evt;
      std::vector< DiHadronCorrelationEvt > evtVec;
      
      //## use correction file
      TrkEff2018PbPb* trkEFF1;
      TrkEff2018PbPb* trkEFF2;
      
      // ## histograms ##
      //~~~> Vertex
      TH1F* hZBestVtx;

      //For MC Gen---
      TH1I* hMultTrgGen;
      TH1I* hMultAssoGen;
     
      //~~~> Global
      TH1I* hCent;
      TH1I* hMultTrg;
      TH1F* hMultTrgCorrtd;
      TH1I* hMultAsso;
      TH1F* hMultAssoCorrtd;
     
      TH1F* hPhidistGen;
      TH1D* hEtadistGen;
      TH1D* hPTdistGen;
      TH1D* hPTdistPlusGen;
      TH1D* hPTdistMinusGen;
      TH1D* hetadistPlusGen;
      TH1D* hetadistMinusGen;
      
      TH1F* hPhidist;
      TH1D* hEtadist;
      TH1D* hPTdist;
      TH1D* hPTdistPlus;
      TH1D* hPTdistMinus;
      TH1D* hetadistPlus;
      TH1D* hetadistMinus;

      TH1F* hchi2;
      TH1F* hnhits;
      TH1F* hpterr;
      TH1F* hdzvtx;

      //~~~> Trigger tracks RAW
      TH1F* hEtaTrk_trg;
      TH1F* hPtTrk_trg;
      TH1F* hPhiTrk_trg;
      
      //~~~> Trigger tracks CORRected
      TH1F* hEtaTrk_trgCorrtd;
      TH1F* hPtTrk_trgCorrtd;
      TH1F* hPhiTrk_trgCorrtd;
      
      //~~~> Associated tracks RAW
      TH1F* hEtaTrk_ass;
      TH1F* hPtTrk_ass;
      TH1F* hPhiTrk_ass;
      //~~~> Associated tracks CORRected
      TH1F* hEtaTrk_assCorrtd;
      TH1F* hPtTrk_assCorrtd;
      TH1F* hPhiTrk_assCorrtd;
    
      //~~~> Signal-----
      //For MC Gen---
      TH2D* hSignalPPGen;
      TH2D* hSignalMMGen;
      TH2D* hSignalPMGen;
      TH2D* hSignalMPGen;
      
      TH2D* hSignalPP;
      TH2D* hSignalMM;
      TH2D* hSignalPM;
      TH2D* hSignalMP;
      
      //Background-------------
      //For MC Gen---
      TH2D* hBackgroundPPGen;
      TH2D* hBackgroundMMGen;
      TH2D* hBackgroundPMGen;
      TH2D* hBackgroundMPGen;
      
      TH2D* hBackgroundPP;
      TH2D* hBackgroundMM;
      TH2D* hBackgroundPM;
      TH2D* hBackgroundMP;
      
      TH1D* hDeltaZvtx;
      

};

//
