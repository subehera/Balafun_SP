//For TwoParticle correlation analysis--
//Clean version---
//By Subash Behera---

// system include files

// CMSSW include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

// user include files
#include "Analyzers/CmsTwoPartCorrAnalysis/interface/CmsTwoPartCorrAnalysis.h"

//
// constructors and destructor
//
CmsTwoPartCorrAnalysis::CmsTwoPartCorrAnalysis(const edm::ParameterSet& iConfig) :

  trackTags (consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  ftrackMC ( consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("trackMC"))  ),
  vtxTags (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),
  centralityTags (consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"))),
  centralityBinTags (consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinSrc"))),
  mvaSrc (consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("mvaSrc"))),
  evtclassifier (iConfig.getUntrackedParameter<int>("evtclassifier")),
  centmin (iConfig.getUntrackedParameter<int>("centmin")),
  centmax (iConfig.getUntrackedParameter<int>("centmax")),
  noffmin(iConfig.getUntrackedParameter<int>("noffmin")),
  noffmax(iConfig.getUntrackedParameter<int>("noffmax")),
  fIsMC (iConfig.getUntrackedParameter<bool>("IsMC") ),
  IsEffCorrection (iConfig.getUntrackedParameter<bool>("IsEffCorrection") ),
  isQAhisto(iConfig.getUntrackedParameter<bool>("isQAhisto")),
  cweight (iConfig.getUntrackedParameter<bool>("cweight")),
  fname (iConfig.getUntrackedParameter<edm::InputTag>("fname")),
  fhpt (iConfig.getUntrackedParameter<edm::InputTag>("fhpt")),
  fmb (iConfig.getUntrackedParameter<edm::InputTag>("fmb")),
  fpix (iConfig.getUntrackedParameter<edm::InputTag>("fpix")),
  fplus (iConfig.getUntrackedParameter<edm::InputTag>("fplus")),
  fminus (iConfig.getUntrackedParameter<edm::InputTag>("fminus")),
  effCorrBinMin(iConfig.getUntrackedParameter< std::vector< int > >("effCorrBinMin")),
  effCorrBinMax(iConfig.getUntrackedParameter< std::vector< int > >("effCorrBinMax")),
  zminVtx (iConfig.getUntrackedParameter<double>("zminVtx")),
  zmaxVtx (iConfig.getUntrackedParameter<double>("zmaxVtx")),
  rhomaxVtx (iConfig.getUntrackedParameter<double>("rhomaxVtx")),
  nTrkAssoToVtx (iConfig.getUntrackedParameter<unsigned int>("nTrkAssoToVtx")),
  selectVtxByMult (iConfig.getUntrackedParameter<bool>("selectVtxByMult")),
  pTmin_trg (iConfig.getUntrackedParameter<double>("pTminTrk_trg")),
  pTmax_trg (iConfig.getUntrackedParameter<double>("pTmaxTrk_trg")),
  pTmin_ass (iConfig.getUntrackedParameter<double>("pTminTrk_ass")),
  pTmax_ass (iConfig.getUntrackedParameter<double>("pTmaxTrk_ass")),
  etamin_trg (iConfig.getUntrackedParameter<double>("etaminTrk_trg")),
  etamax_trg (iConfig.getUntrackedParameter<double>("etamaxTrk_trg")),
  etamin_ass (iConfig.getUntrackedParameter<double>("etaminTrk_ass")),
  etamax_ass (iConfig.getUntrackedParameter<double>("etamaxTrk_ass")),
  bkgFactor(iConfig.getUntrackedParameter<unsigned int>("bkgFactor")),
  nEtaBins(iConfig.getUntrackedParameter<int>("nEtaBins")),
  nPhiBins(iConfig.getUntrackedParameter<int>("nPhiBins"))
{
  //file acc & eff
  TString filename(fname.label().c_str());
  feff = 0x0;
  
  if(cweight && !filename.IsNull())
    {
      edm::FileInPath fip(Form("Analyzers/CmsTwoPartCorrAnalysis/data/EFF/gentrk/%s",filename.Data()));
      feff = new TFile(fip.fullPath().c_str(),"READ");
      heff.resize(feff->GetNkeys());
      for(unsigned int ik = 0; ik < heff.size(); ++ik)
	{
          heff[ik] = (TH2D*) feff->Get(feff->GetListOfKeys()->At(ik)->GetName());
	}
    }
  else
    {
      cweight = false;
      edm::LogWarning ("Cannot open file") <<"Invalid efficiency correction file!";
    }
  
  evt = new DiHadronCorrelationEvt(2, 2);
  //----------------------for gentrk------
  TString f_pt(fhpt.label().c_str());
  edm::FileInPath f1(Form("Analyzers/CmsTwoPartCorrAnalysis/data/EFF/effHpT/%s",f_pt.Data()));
  
  TString f_mb(fmb.label().c_str());
  edm::FileInPath f2(Form("Analyzers/CmsTwoPartCorrAnalysis/data/EFF/effMB/%s",f_mb.Data()));
  
  TString f_pix(fpix.label().c_str());
  edm::FileInPath f3(Form("Analyzers/CmsTwoPartCorrAnalysis/data/EFF/effPix/%s",f_pix.Data()));
  
  TString f_plus(fplus.label().c_str());
  edm::FileInPath f4(Form("Analyzers/CmsTwoPartCorrAnalysis/data/EFF/plus/%s",f_plus.Data()));
  
  TString f_minus(fminus.label().c_str());
  edm::FileInPath f5(Form("Analyzers/CmsTwoPartCorrAnalysis/data/EFF/minus/%s",f_minus.Data()));

  //by Subas-------
  trkEFF1 = new TrkEff2018PbPb(  "generalMB+",  false , f1.fullPath(), f2.fullPath(), f4.fullPath(), f5.fullPath(), f3.fullPath());
  trkEFF2 = new TrkEff2018PbPb(  "generalMB-",  false , f1.fullPath(), f2.fullPath(), f4.fullPath(), f5.fullPath(), f3.fullPath());
  //trkEff = new TrkEff2018PbPb(  "general",  false , f1.fullPath(), f2.fullPath(), f4.fullPath(), f5.fullPath(), f3.fullPath());
  
  
  //--------------------------------------
  TH1::SetDefaultSumw2();
  
  // Now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  // Histograms
  TFileDirectory fVtxHist  = fs->mkdir("Vertex");
  hZBestVtx   = fVtxHist.make<TH1F>("hZvtx", "", 600, -30., 30.);
  
  TFileDirectory fGlobalHist  = fs->mkdir("Global");
  
  const Int_t nptBin = 15;
  Double_t ptbining[nptBin+1] = {0.5, 0.6,  0.7, 0.8, 0.9,1.0, 1.1, 1.2, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.0};

  int nMultbin = 10000;
  int multbinLow = -0.5;
  int multbinHigh = 9999.5;
  
  hCent               = fGlobalHist.make<TH1I>("hCent", "",  200, -0.5, 199.5);
  
  if(fIsMC){
    hMultTrgGen            = fGlobalHist.make<TH1I>("hMult_trgGen",
						    Form("%1.1f < p_{T} < %1.1f GeV/c", pTmin_trg, pTmax_trg),
						    nMultbin, multbinLow, multbinHigh);
    hMultAssoGen           = fGlobalHist.make<TH1I>("hMultAssoGen",
						    Form("%1.1f < p_{T} < %1.1f GeV/c", pTmin_trg, pTmax_trg),
						    nMultbin, multbinLow, multbinHigh);
  }
  
  hMultTrg            = fGlobalHist.make<TH1I>("hMult_trg",
					     Form("%1.1f < p_{T} < %1.1f GeV/c", pTmin_trg, pTmax_trg),
					       nMultbin, multbinLow, multbinHigh);
  hMultTrgCorrtd      = fGlobalHist.make<TH1F>("hMultTrgCorrtd",
					       Form("%1.1f < p_{T} < %1.1f GeV/c", pTmin_trg, pTmax_trg),
					       nMultbin, multbinLow, multbinHigh);
  hMultAsso           = fGlobalHist.make<TH1I>("hMultAsso",
					     Form("%1.1f < p_{T} < %1.1f GeV/c", pTmin_trg, pTmax_trg),
					     nMultbin, multbinLow, multbinHigh);
  hMultAssoCorrtd     = fGlobalHist.make<TH1F>("hMultAssoCorrtd",
					       Form("%1.1f < p_{T} < %1.1f GeV/c", pTmin_trg, pTmax_trg),
					       nMultbin, multbinLow, multbinHigh);
  if(isQAhisto){
    
    if(fIsMC){     
      hPhidistGen          = fGlobalHist.make<TH1F>("hPhidistGen", "", nPhiBins, -TMath::Pi(),  TMath::Pi() );
      hEtadistGen          = fGlobalHist.make<TH1D>("hEtadistGen", "#eta all", nEtaBins, etamin_trg, etamax_trg );
      hPTdistGen           = fGlobalHist.make<TH1D>("hPTdistGen", "p_{T} all", nptBin, ptbining);
      hPTdistPlusGen       = fGlobalHist.make<TH1D>("hPTdistPlusGen", "p_{T}, h^{+}", nptBin, ptbining);
      hPTdistMinusGen      = fGlobalHist.make<TH1D>("hPTdistMinusGen", "p_{T}, h^{-}", nptBin, ptbining);
      hetadistPlusGen      = fGlobalHist.make<TH1D>("hetadistPlusGen", "#eta h^{+}", nEtaBins, etamin_trg, etamax_trg);
      hetadistMinusGen     = fGlobalHist.make<TH1D>("hetadistMinusGen", "#eta h^{-}", nEtaBins, etamin_trg, etamax_trg);
    }
    
    hPhidist          = fGlobalHist.make<TH1F>("hPhidist", "", nPhiBins, -TMath::Pi(),  TMath::Pi() );
    hEtadist          = fGlobalHist.make<TH1D>("hEtadist", "#eta all", nEtaBins, etamin_trg, etamax_trg );
    hPTdist           = fGlobalHist.make<TH1D>("hPTdist", "p_{T} all", nptBin, ptbining);
    hPTdistPlus       = fGlobalHist.make<TH1D>("hPTdistPlus", "p_{T}, h^{+}", nptBin, ptbining);
    hPTdistMinus      = fGlobalHist.make<TH1D>("hPTdistMinus", "p_{T}, h^{-}", nptBin, ptbining);
    hetadistPlus      = fGlobalHist.make<TH1D>("hetadistPlus", "#eta h^{+}", nEtaBins, etamin_trg, etamax_trg);
    hetadistMinus     = fGlobalHist.make<TH1D>("hetadistMinus", "#eta h^{-}", nEtaBins, etamin_trg, etamax_trg);
    
    //Ttrigger track raw histo--------------------------------------
    TFileDirectory fTrkTrgHist  = fs->mkdir("TrgTracksRaw");
    
    int nbins_trg = (pTmax_trg - pTmin_trg )/0.01;
    
    hEtaTrk_trg = fTrkTrgHist.make<TH1F>("hEtaTrk_trg", Form("%1.1f < p_{T} < %1.1f GeV/c", pTmin_trg, pTmax_trg),
					 300, -3., 3.);
    hPtTrk_trg  = fTrkTrgHist.make<TH1F>("hPtTrk_trg", Form("%1.1f<p_{T}<%1.1f GeV/c", pTmin_trg, pTmax_trg ),
					 nbins_trg,  pTmin_trg, pTmax_trg );
    hPhiTrk_trg  = fTrkTrgHist.make<TH1F>("hPhiTrk_trg",
					  Form("%1.1f<p_{T}<%1.1f GeV/c", pTmin_trg, pTmax_trg ),
					  640, -3.2, 3.2 );
    
    //Trigger track corrected histos-------------------------------
    TFileDirectory fTrkCorrTrgHist  = fs->mkdir("TrgTracksCorr");
    hEtaTrk_trgCorrtd = fTrkCorrTrgHist.make<TH1F>("hEtaTrk_trgCorrtd",
						   Form("%1.1f<p_{T}<%1.1f GeV/c", pTmin_trg, pTmax_trg ),
						   300, -3., 3.);
    hPtTrk_trgCorrtd  = fTrkCorrTrgHist.make<TH1F>("hPtTrk_trgCorrtd",
						   Form("%1.1f<p_{T}<%1.1f GeV/c", pTmin_trg, pTmax_trg ),
						   nbins_trg, pTmin_trg, pTmax_trg );
    hPhiTrk_trgCorrtd = fTrkCorrTrgHist.make<TH1F>("hPhiTrk_trgCorrtd",
						   Form("%1.1f<p_{T}<%1.1f GeV/c", pTmin_trg, pTmax_trg ),
						   640, -3.2, 3.2);
    //Associated track raw histos--------------------------------------
    TFileDirectory fTrkAssHist  = fs->mkdir("AssTracksRaw");
    int nbins_ass = (int) (pTmax_ass - pTmin_ass)/0.01;
    hEtaTrk_ass = fTrkAssHist.make<TH1F>("hEtaTrk_ass",
					 Form("%1.1f<$p_{T}<%1.1f GeV/c", pTmin_ass, pTmax_ass),
					 300, -3., 3.);
    hPtTrk_ass  = fTrkAssHist.make<TH1F>("hPtTrk_ass",
					 Form("%1.1f<$p_{T}<%1.1f GeV/c", pTmin_ass, pTmax_ass),
					 nbins_ass,  pTmin_ass, pTmax_ass);
    hPhiTrk_ass = fTrkAssHist.make<TH1F>("hPhiTrk_ass",
					 Form("%1.1f<$p_{T}<%1.1f GeV/c", pTmin_ass, pTmax_ass),
					 640, -3.2, 3.2);

    //Associated track corrected histos--------------------------------------
    TFileDirectory fTrkCorrAssHist  = fs->mkdir("AssTracksCorr");
    hEtaTrk_assCorrtd = fTrkCorrAssHist.make<TH1F>("hEtaTrk_assCorrtd",
						   Form("%1.1f<$p_{T}<%1.1f GeV/c", pTmin_ass, pTmax_ass),
						   300, -3., 3.);
    hPtTrk_assCorrtd  = fTrkCorrAssHist.make<TH1F>("hPtTrk_assCorrtd",
						   Form("%1.1f<$p_{T}<%1.1f GeV/c", pTmin_ass, pTmax_ass),
						   nbins_ass, pTmin_ass, pTmax_ass);
    hPhiTrk_assCorrtd = fTrkCorrAssHist.make<TH1F>("hPhiTrk_assCorrtd",
						   Form("%1.1f<$p_{T}<%1.1f GeV/c", pTmin_ass, pTmax_ass),
						   640, -3.2, 3.2);
  }//if(isQA)-------
  
  
  //Define eta and phi bin limits---
  double etaW = (etamax_trg - etamin_ass - etamin_trg + etamax_ass) / nEtaBins;
  double phiW = 2.0*(TMath::Pi())/nPhiBins;
  double minEta = etamin_trg - etamax_ass - etaW/2;
  double maxEta = etamax_trg - etamin_ass + etaW/2.;
  double minPhi = -(TMath::Pi() - phiW)/2.0;
  double maxPhi = (TMath::Pi() * 3.0 - phiW)/2.0;
  
  TFileDirectory fSignalHist      = fs->mkdir("Signal");
  
  if(fIsMC){
    
    hSignalPPGen = fSignalHist.make<TH2D>("signal_PP_Gen",
				       Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
					    pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
				       nEtaBins+1, minEta, maxEta, nPhiBins-1, minPhi, maxPhi);
    
    hSignalPMGen = fSignalHist.make<TH2D>("signal_PM_Gen",
				       Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
					    pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
				       nEtaBins + 1, minEta, maxEta, nPhiBins-1, minPhi, maxPhi);
    hSignalMPGen = fSignalHist.make<TH2D>("signal_MP_Gen",
				       Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
					    pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
				       nEtaBins + 1, minEta, maxEta, nPhiBins - 1, minPhi, maxPhi);
    
    hSignalMMGen = fSignalHist.make<TH2D>("signal_MM_Gen",
				       Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
					    pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
				       nEtaBins + 1, minEta, maxEta, nPhiBins - 1, minPhi, maxPhi);
 }
  
  hSignalPP = fSignalHist.make<TH2D>("signal_PP",
				     Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
					  pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
				     nEtaBins+1, minEta, maxEta, nPhiBins-1, minPhi, maxPhi);
  
  hSignalPM = fSignalHist.make<TH2D>("signal_PM",
				     Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
					  pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
				     nEtaBins + 1, minEta, maxEta, nPhiBins-1, minPhi, maxPhi);
  hSignalMP = fSignalHist.make<TH2D>("signal_MP",
				     Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
					  pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
				     nEtaBins + 1, minEta, maxEta, nPhiBins - 1, minPhi, maxPhi);
  
  hSignalMM = fSignalHist.make<TH2D>("signal_MM",
				     Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
					  pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
				     nEtaBins + 1, minEta, maxEta, nPhiBins - 1, minPhi, maxPhi);
  
  
  //For Background-------------
  TFileDirectory fBackgroundHist  = fs->mkdir("Background");
  
  if(fIsMC){
    hBackgroundPPGen = fBackgroundHist.make<TH2D>("background_PP_Gen",
						  Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
						       pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
						  nEtaBins + 1, minEta, maxEta, nPhiBins - 1, minPhi, maxPhi );
    hBackgroundPMGen = fBackgroundHist.make<TH2D>("background_PM_Gen",
						  Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
						       pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
						  nEtaBins + 1, minEta, maxEta, nPhiBins - 1, minPhi, maxPhi );
    hBackgroundMPGen = fBackgroundHist.make<TH2D>("background_MP_Gen",
						  Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
						       pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
						  nEtaBins + 1, minEta, maxEta, nPhiBins - 1, minPhi, maxPhi );
    
    hBackgroundMMGen = fBackgroundHist.make<TH2D>("background_MM_Gen",
						  Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
						       pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
						  nEtaBins + 1, minEta, maxEta, nPhiBins - 1, minPhi, maxPhi );
    
    
  }
  
  hBackgroundPP = fBackgroundHist.make<TH2D>("background_PP",
					     Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
						  pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
					     nEtaBins + 1, minEta, maxEta, nPhiBins - 1, minPhi, maxPhi );
  hBackgroundPM = fBackgroundHist.make<TH2D>("background_PM",
					     Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
						  pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
					     nEtaBins + 1, minEta, maxEta, nPhiBins - 1, minPhi, maxPhi );
  hBackgroundMP = fBackgroundHist.make<TH2D>("background_MP",
					     Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
						  pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
					     nEtaBins + 1, minEta, maxEta, nPhiBins - 1, minPhi, maxPhi );
  
  hBackgroundMM = fBackgroundHist.make<TH2D>("background_MM",
					     Form("%1.1f<p_{T}^{trg}<%1.1f GeV/c, %1.1f<p_{T}^{ass}<%1.1f GeV/c;#Delta#eta;#Delta#phi",
						  pTmin_trg, pTmax_trg, pTmin_ass, pTmax_ass),
					     nEtaBins + 1, minEta, maxEta, nPhiBins - 1, minPhi, maxPhi ); 
  
  
  
}

CmsTwoPartCorrAnalysis::~CmsTwoPartCorrAnalysis()
{
 
   delete evt;
}


//
// member functions
//

// ------------ method called for each event  ------------
void CmsTwoPartCorrAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // ----- centrality selection -----
   // Get calo centrality collection by token
   edm::Handle< reco::Centrality > centrality;
   iEvent.getByToken(centralityTags, centrality);
   // Get calo centrality bin by token
   edm::Handle< int > cbin;
   iEvent.getByToken(centralityBinTags, cbin);
   int centBin = *cbin;
   if(centBin < 0)
   {
       edm::LogWarning ("Invalid value") <<"Invalid centrality value";
   }

   // ----- Vertex selection---------------------
   // Get vertex collection by token
   edm::Handle< reco::VertexCollection > vertices;
   iEvent.getByToken(vtxTags, vertices);
   if( !vertices->size() )
   {
      edm::LogWarning ("Missing Collection") <<"Invalid or empty vertex collection!";
      return;
   }

   //---Loop over vertices--------------------------------------
   xBestVtx = -999.;
   yBestVtx = -999.;
   zBestVtx = -999.;
   rhoBestVtx = -999.;
   xBestVtxError = -999.;
   yBestVtxError = -999.;
   zBestVtxError = -999.;
   
   LoopVertices(iEvent, iSetup);
   
   if( zBestVtx < zminVtx || zBestVtx  > zmaxVtx ) return; //cut on vertex Z position
   if( rhoBestVtx > rhomaxVtx ) return; //cut on vertex XY position
   
   hZBestVtx->Fill( zBestVtx );
   
   // ----- Define event classification (either centrality or Ntrk^off) -----
   int evtclass = -1;
   //int evtclass = 0;
   switch(evtclassifier)
     {
     case 0:
       evtclass = centBin;
       if(evtclass < centmin*2 || evtclass >= centmax*2) return;
       break;
     case 1:
       evtclass = -999.0; //set by hand for now
       if(evtclass < noffmin || evtclass >= noffmax)
	 return;
     default:
       evtclass = -1;
     }

   nTrkTot_trgGen    = 0;
   nTrkTot_assGen    = 0;
   
   nTrkTot_trg       = 0;
   nTrkTot_trgCorrtd = 0;
   nTrkTot_ass       = 0;
   nTrkTot_assCorrtd = 0;
   
   // ----- Track selection -----
   if(fIsMC){ //Only Gen level info---
     LoopTracksMC(iEvent, iSetup, true,  evtclass, fIsMC); //trigger tracks
     LoopTracksMC(iEvent, iSetup, false, evtclass, fIsMC); //associated tracks
     getSpherocity(iEvent, iSetup, evtclass, fIsMC, spValue);//get spherocity value
   }
   //Either MC reco, reco+corr or data---
   LoopTracks(iEvent, iSetup, true,  evtclass); //trigger tracks
   LoopTracks(iEvent, iSetup, false, evtclass); //associated tracks
   // ----- Fill and push evt containers -----
   evt->run   = iEvent.id().run();
   evt->event = iEvent.id().event();
   evt->zvtx  = zBestVtx;

   evtVec.push_back(*evt);

   hCent->Fill(centBin);
   // ----- Reset evt container -----
   evt->reset();
}

void CmsTwoPartCorrAnalysis::getSpherocity(const edm::Event& iEvent, const edm::EventSetup& iSetup,
					int evtclass, bool isMC, double& spValue)
{

edm::Handle< reco::GenParticleCollection > tracksGen;
iEvent.getByToken(ftrackMC, tracksGen);
  if( !tracksGen->size() )
    {
      edm::LogWarning ("Missing MC Gen Collection") <<"Invalid or empty MC-Gen track collection!";
      return;
    }

  float spherocity = -10.0;
  float pFull = 0;
  const float fSizeStepESA = 1;
  float Spherocity = 2;
  float sumapt = 0;      


  for(Int_t i = 0; i < 360/(fSizeStepESA); ++i)
    {
      float numerator = 0;
      float phiparam  = 0;
      float nx = 0;
      float ny = 0;
      phiparam=( (TMath::Pi()) * i * fSizeStepESA ) / 180; // parametrization of the angle
      nx = TMath::Cos(phiparam);            // x component of an unitary vector n
      ny = TMath::Sin(phiparam);            // y component of an unitary vector n
	  
 for (reco::GenParticleCollection::const_iterator iii = tracksGen->begin();iii != tracksGen->end(); ++iii)
	{

      if( iii->status() != 1 ) continue;
      if( iii->charge() == 0 ) continue;
		
	if(abs(eta)< 0.8) 
	    {
                  pt_pp[iii]= iii->pt();
		  eta_pp[iii]= iii->eta();
		  phi_pp[iii]= iii->phi();
		  if (i==0) sumapt +=  pt_pp[iii];
		  Float_t pxA = pt_pp[iii] * TMath::Cos( phi_pp[iii] );
		  Float_t pyA = pt_pp[iii] * TMath::Sin( phi_pp[iii] );
		  numerator += TMath::Abs( ny * pxA - nx * pyA );//product between p  projection in XY plane and the unitary vector
		    
	}// eta selection
	}// track loop
   pFull = TMath::Power((numerador / sumapt),2 );
   if(pFull < Spherocity)//maximization of pFull
	{
	  Spherocity = pFull;
	}	
	    
    }// unit vector loop 
 spherocity=((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;
 spValue = spherocity;
}
// ------------ method called once each job just before starting event loop  ------------
void CmsTwoPartCorrAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void CmsTwoPartCorrAnalysis::endJob()
{
   std::cout<< "Start sorting the events!" << std::endl;
   std::sort(evtVec.begin(), evtVec.end());
   std::cout<< "Finish sorting the events!" << std::endl;

   std::cout<< "Total of " << evtVec.size() << " events are selected!" << std::endl;

   std::cout<< "Start running correlation analysis!" << std::endl;
   for( unsigned int i = 0; i < evtVec.size(); i++ )
   {
     //if( i % 100 == 0 ) std::cout << "Processing " << i << "th event" << std::endl;
     //std::cout<<"Signal event :  " << i << std::endl<<std::endl<<std::endl;
     if(fIsMC) FillHistsSignal( i, kTRUE );
     FillHistsSignal( i, kFALSE );
      unsigned int mixstart = i - bkgFactor/2;
      unsigned int mixend   = i + bkgFactor/2 + 1;
      
      if(i < bkgFactor)
	{
	  mixstart = 0;
	  mixend   = 2*bkgFactor + 1;
	}
      else if(i > evtVec.size() - bkgFactor - 1)
	{
	  mixstart = evtVec.size() - 2*bkgFactor - 1;
	  mixend   = evtVec.size();
	}
      
      if( mixend > evtVec.size() )
	mixend = evtVec.size();
      std::cout << mixstart << " : " << mixend << std::endl;
      for( unsigned int j = mixstart; j < mixend; j++ )
	{
          if(i == j) continue;
          //if(Nmix >= 10) continue;
	  //if(Nmix >= 5) continue;
	  
          //std::cout << i << " : " << j << std::endl;
          double deltazvtx = evtVec[i].zvtx - evtVec[j].zvtx;
          if(fabs(deltazvtx) > 2.0) continue;
	  //if(fabs(deltazvtx) > 0.5) continue;
	  if(fIsMC)FillHistsBackground( i, j, kTRUE );
	  FillHistsBackground( i, j, kFALSE );
	}
   }
   
   
   std::cout<< "Finish running correlation analysis!" << std::endl;
   
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void CmsTwoPartCorrAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//=========================================================================================

void CmsTwoPartCorrAnalysis::LoopVertices(const edm::Event& iEvent,
                                 const edm::EventSetup& iSetup)
{
  
  
  
   edm::Handle< reco::VertexCollection > vertices;
   iEvent.getByToken(vtxTags, vertices);
   if(!vertices->size())
   {
      std::cout<<"Invalid or empty vertex collection!"<<std::endl;
      return;
   }

   reco::VertexCollection recoVertices = *vertices;

   if( selectVtxByMult )
   {
       std::sort( recoVertices.begin(),
                  recoVertices.end(),
                  [](const reco::Vertex &a, const reco::Vertex &b)
                  {
                     if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2();
                          return a.tracksSize() > b.tracksSize();
                  }
                );
   }

   for( reco::VertexCollection::const_iterator itVtx = recoVertices.begin();
        itVtx != recoVertices.end();
        ++itVtx )
   {
        // Drop fake vertex and vertex with less than 2 tracks attached to it
        if( !itVtx->isFake() && itVtx->tracksSize() >= nTrkAssoToVtx )
        {
            // x,y,z vertex position
            double xVtx = itVtx->x();
            double yVtx = itVtx->y();
            double zVtx = itVtx->z();
            // x,y,z vertex position error
            double xVtxError = itVtx->xError();
            double yVtxError = itVtx->yError();
            double zVtxError = itVtx->zError();
            // Radial vertex position
            double rho = sqrt(xVtx*xVtx + yVtx*yVtx);
            // Increase N valid vertex in the collection
            ++nVtx;

            //Get the first vertex as the best one (greatest sum p_{T}^{2})
            if( itVtx == recoVertices.begin() )
            {
                xBestVtx = xVtx;
                yBestVtx = yVtx;
                zBestVtx = zVtx;
		rhoBestVtx = rho;
                xBestVtxError = xVtxError;
                yBestVtxError = yVtxError;
                zBestVtxError = zVtxError;
            }
        }
   }
   
}
//=========================================================================================
//=========================================================================================
void CmsTwoPartCorrAnalysis::LoopTracksMC(const edm::Event& iEvent, const edm::EventSetup& iSetup,
					  bool istrg, int evtclass, bool isMC )
{
  
  edm::Handle< reco::GenParticleCollection > tracksGen;
  iEvent.getByToken(ftrackMC, tracksGen);
  if( !tracksGen->size() )
    {
      edm::LogWarning ("Missing MC Gen Collection") <<"Invalid or empty MC-Gen track collection!";
      return;
    }
  
  // Loop over tracks
  for( reco::GenParticleCollection::const_iterator itTrk = tracksGen->begin();itTrk != tracksGen->end(); ++itTrk )
    {
      
      // Get eta, pt, phi and charge of the track
      if( itTrk->status() != 1 ) continue;
      if( itTrk->charge() == 0 ) continue;
      // Get eta, pt, phi and charge of the track
      double eta      = itTrk->eta();
      double pt       = itTrk->pt();
      double phi      = itTrk->phi();
      int charge      = itTrk->charge();
      
      float eff1 = 1., eff2 = 1.;
      if( istrg ){ //for trigger particle--
	if( pt <= pTmin_trg || pt > pTmax_trg ) continue;
	if( eta < etamin_trg || eta > etamax_trg ) continue;
	
	hPhidistGen->Fill(phi);
	hEtadistGen->Fill(eta);
	hPTdistGen->Fill(pt);
	if(charge>0){
	  hPTdistPlusGen->Fill(pt);
	  hetadistPlusGen->Fill(eta);
	}
	else {
	  hPTdistMinusGen->Fill(pt);
	  hetadistMinusGen->Fill(eta);
	}
	
      } //else asso.
      else{
	if( pt <= pTmin_ass || pt > pTmax_ass ) continue;
	if( eta < etamin_ass || eta > etamax_ass ) continue;
      }
      
      AssignpTbins(pt, eta, phi, charge, eff1, eff2,  istrg, isMC );
      
      //std::cout<<" pt is="<<pt<<std::endl;
    }
  
  //Fill trk histograms
  
  if(istrg)
    { 
      hMultTrgGen->Fill(nTrkTot_trgGen);
      (evt->nMultCorrVect_trg)[0] = nTrkTot_trgGen;
    }
  else
    {
      hMultAssoGen->Fill(nTrkTot_assGen);
      (evt->nMultCorrVect_ass)[0] = nTrkTot_assGen;
    }
  
}
//==========================================================================================
  
//=========================================================================================
void CmsTwoPartCorrAnalysis::LoopTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup,
                               bool istrg, int evtclass)
{

   edm::Handle<std::vector<float>> mvaoutput;
   iEvent.getByToken(mvaSrc, mvaoutput);

// Get track collection by token
   edm::Handle< reco::TrackCollection > tracks;
   iEvent.getByToken(trackTags, tracks);
   if( !tracks->size() )
   {
       edm::LogWarning ("Missing Collection") <<"Invalid or empty track collection!";
       return;
   }

   // Loop over tracks-----------
   int it =0;
   
   for( reco::TrackCollection::const_iterator itTrk = tracks->begin(); itTrk != tracks->end(); ++itTrk ){

     // Select tracks based on proximity to best vertex
     math::XYZPoint bestvtx( xBestVtx, yBestVtx,zBestVtx );
     double dzvtx    = itTrk->dz(bestvtx);
     double dxyvtx   = itTrk->dxy(bestvtx);
     double dzerror  = sqrt(itTrk->dzError()*itTrk->dzError() + zBestVtxError*zBestVtxError);
     double dxyerror = sqrt(itTrk->d0Error()*itTrk->d0Error() + xBestVtxError*yBestVtxError);
     double pterror  = itTrk->ptError();
     // Get eta, pt, phi and charge of the track
     double eta      = itTrk->eta();
     double pt       = itTrk->pt();
     double phi      = itTrk->phi();
     int charge      = itTrk->charge();
     // HI specific cuts
     double chi2n    = itTrk->normalizedChi2();
     double nlayers  = itTrk->hitPattern().trackerLayersWithMeasurement();
     chi2n           = chi2n/nlayers;
     int nHits       = itTrk->numberOfValidHits();
     int algo        = itTrk->originalAlgo();
     
     float trkMVA = (*mvaoutput)[it];    
     it++;  
     
     // Select track based on quality
     if( !itTrk->quality(reco::TrackBase::highPurity) ) continue;
     if( charge == 0 ) continue;
     if( nHits < 11) continue ;  
     if((pterror)/ pt >= 0.1 ) continue;
     if(chi2n >= 0.18 ) continue ;                                                                        
     if(( fabs(dxyvtx / dxyerror)) >= 3.0) continue;
     if(fabs(dzvtx / dzerror) >= 3.0 ) continue;
     if( pt <= 0.5 ) continue;
     if( eta < -2.4 || eta > 2.4 ) continue;
     if(algo==6 && trkMVA < 0.98)continue;
     
     
     float eff1 = 1., eff2 = 1.;
     if(IsEffCorrection){
     eff1  = trkEFF1->getCorrection(pt, eta, evtclass);
     eff2  = trkEFF2->getCorrection(pt, eta, evtclass);
     }
     
     if( istrg ){ //for trigger particle--
       if( pt <= pTmin_trg || pt > pTmax_trg ) continue;
       if( eta < etamin_trg || eta > etamax_trg ) continue;

       hPhidist->Fill(phi);
       hEtadist->Fill(eta);
       hPTdist->Fill(pt);
       if(charge>0){
	 hPTdistPlus->Fill(pt);
	 hetadistPlus->Fill(eta);
       }
       else {
	 hPTdistMinus->Fill(pt);
	 hetadistMinus->Fill(eta);
       }
     } //else asso.
     else{
       if( pt <= pTmin_ass || pt > pTmax_ass ) continue;
       if( eta < etamin_ass || eta > etamax_ass ) continue;
     }
     
     AssignpTbins(pt, eta, phi, charge, eff1, eff2,  istrg, kFALSE );
     //std::cout<<" reco pt is="<<pt<<std::endl;
     
   }
   
   // Fill trk histograms
   if(istrg)
     {     
       hMultTrg->Fill(nTrkTot_trg);
       hMultTrgCorrtd->Fill(nTrkTot_trgCorrtd);
       (evt->nMultCorrVect_trg)[1] = nTrkTot_trgCorrtd;
     }
   else
     {
       hMultAsso->Fill(nTrkTot_ass);
       hMultAssoCorrtd->Fill(nTrkTot_assCorrtd);
       (evt->nMultCorrVect_ass)[1] = nTrkTot_assCorrtd;
     }
   
}
//=========================================================================================

//=========================================================================================
double CmsTwoPartCorrAnalysis::GetEffWeight(double eta, double pt, int evtclass)
{
  double effweight = 1.0;
  if(evtclass == -1)
    {
      effweight = heff[0]->GetBinContent(heff[0]->FindBin(eta,pt));
    }
  else
    {
      int centIdx = 0;
      for(int icent = 0; icent < static_cast<int>(effCorrBinMin.size()); ++icent)
	{
	  if(evtclass >= effCorrBinMin[icent]*2 && evtclass < effCorrBinMax[icent]*2)
	    {
	      centIdx = icent;
	      continue;
	    }
	}
      effweight = heff[centIdx]->GetBinContent(heff[centIdx]->GetXaxis()->FindBin(eta), heff[centIdx]->GetYaxis()->FindBin(pt));
    }
  return effweight;
}

//=========================================================================================
//
//=========================================================================================
void CmsTwoPartCorrAnalysis::AssignpTbins(double pt,  double eta, double phi, int charge, double eff1, double eff2, bool istrg, bool isGen )
{

  TLorentzVector pvector;
  pvector.SetPtEtaPhiM(pt, eta, phi, 0.140);
  
  double EffValue[2] = { eff1, eff2 };
  int idx;
  charge > 0 ? idx = 0 : idx = 1; 

  if(istrg)
    {
      if(isGen){
	(evt->pVect_trg[0]).push_back(pvector);
	(evt->chgVect_trg[0]).push_back(charge);
	(evt->effVect_trg[0]).push_back(1.0);
	nTrkTot_trgGen += 1.;
      }
      else{
	if(isQAhisto){
	  hEtaTrk_trg->Fill(eta);
	  hPtTrk_trg->Fill(pt);
	  hPhiTrk_trg->Fill(phi);
	  //corrected----
	  hEtaTrk_trgCorrtd->Fill(eta, EffValue[idx]);
	  hPtTrk_trgCorrtd->Fill(pt, EffValue[idx]);
	  hPhiTrk_trgCorrtd->Fill(phi);
	}//QA--
	//std::cout<<" pt value=" << pt <<std::endl;
	(evt->pVect_trg[1]).push_back(pvector);
	(evt->chgVect_trg[1]).push_back(charge);
	(evt->effVect_trg[1]).push_back(1.0/EffValue[idx]);
	nTrkTot_trg       += 1.;
	nTrkTot_trgCorrtd += EffValue[idx];
      }
      
    }//if trg---
  else  
    {
      if(isGen){
	(evt->pVect_ass[0]).push_back(pvector);
	(evt->chgVect_ass[0]).push_back(charge);
	(evt->effVect_ass[0]).push_back(1.0);
	nTrkTot_assGen  += 1.;
      }
      else{
	if(isQAhisto){
	  hPtTrk_ass->Fill(pt);
	  hEtaTrk_ass->Fill(eta);
	  hPhiTrk_ass->Fill(phi);
	  
	  hPtTrk_assCorrtd->Fill(pt, EffValue[idx]);
	  hEtaTrk_assCorrtd->Fill(eta, EffValue[idx]);
	  hPhiTrk_assCorrtd->Fill(phi);
	}
	
	(evt->pVect_ass[1]).push_back(pvector);
	(evt->chgVect_ass[1]).push_back(charge);
	(evt->effVect_ass[1]).push_back(1.0/EffValue[idx]);
	nTrkTot_ass       += 1.;
	nTrkTot_assCorrtd += EffValue[idx];
      }
      
    }//else asso--
  
}

//--------------------------------------------------------------------------------------
double CmsTwoPartCorrAnalysis::GetDeltaEta(double eta_trg, double eta_ass)
  
{
  double deltaEta = eta_ass - eta_trg;
  //double deltaEta = eta_trg - eta_ass;
  return deltaEta;
}

double CmsTwoPartCorrAnalysis::GetDeltaPhi(double phi_trg, double phi_ass)
{
  //double deltaPhi = phi_trg - phi_ass;
  double deltaPhi = phi_ass - phi_trg;
  
  if(deltaPhi > 1.5*TMath::Pi())
    deltaPhi = deltaPhi - 2.0*TMath::Pi();
  
  else if(deltaPhi < -1.0*TMath::Pi() / 2.0)
    deltaPhi = deltaPhi + 2.0*TMath::Pi();
  
  return deltaPhi;
  
  
}

double CmsTwoPartCorrAnalysis::GetDPhiStar(double phi1, double pt1, int charge1, double phi2, double pt2, int charge2, double radius, double bSign)
{
  // calculates dphistar
  double dphistar = phi1 - phi2 - charge1*bSign*TMath::ASin(0.003*radius / pt1) + charge2*bSign*TMath::ASin(0.003*radius / pt2);
  
  static const double kPi = TMath::Pi();
  
  /*
    if (dphistar > 1.5 * kPi)
    dphistar = dphistar - 2.0*kPi;
    else if (dphistar < -0.5 * kPi)
    dphistar = dphistar + 2.0*kPi;*/
  
  if (dphistar > kPi)
    dphistar = kPi * 2.0 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2.0 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2.0 - dphistar;
  
  return dphistar;
}


void CmsTwoPartCorrAnalysis::FillHistsSignal(int ievt, bool isGen)
{
  int jj;
  isGen ? jj = 0 : jj = 1;
  
  unsigned int ntrgsize = evtVec[ievt].pVect_trg[jj].size();
  unsigned int nasssize = evtVec[ievt].pVect_ass[jj].size();
  double nMult_corr_trg = evtVec[ievt].nMultCorrVect_trg[jj];
 
  for( unsigned int itrg = 0; itrg < ntrgsize; itrg++ )
    {
      TLorentzVector pvector_trg = (evtVec[ievt].pVect_trg[jj])[itrg];
      double effweight_trg       = (evtVec[ievt].effVect_trg[jj])[itrg];
      int chg_trg                = (evtVec[ievt].chgVect_trg[jj])[itrg];
      double eta_trg             = pvector_trg.Eta();
      double phi_trg             = pvector_trg.Phi();
      double pt_trg              = pvector_trg.Pt();
      //std::cout<<" pt_trg value=" << pt_trg <<std::endl;
      for( unsigned int jass = 0; jass < nasssize; jass++ )
	{
	  TLorentzVector pvector_ass = (evtVec[ievt].pVect_ass[jj])[jass];
	  double effweight_ass       = (evtVec[ievt].effVect_ass[jj])[jass];
	  int chg_ass                = (evtVec[ievt].chgVect_ass[jj])[jass];
	  double eta_ass             = pvector_ass.Eta();
	  double phi_ass             = pvector_ass.Phi();
	  double pt_ass              = pvector_ass.Pt();
	  
	  double deltaPhi            = GetDeltaPhi(phi_trg, phi_ass);
	  double deltaPhi2           = GetDeltaPhi(phi_ass, phi_trg);
	  double deltaEta            = GetDeltaEta(eta_trg, eta_ass);
	  
	  //Skip the loop when trg, ass particles are the same
	  if( deltaEta == 0.0  &&  deltaPhi == 0.0 && deltaPhi2 == 0.0 && pt_trg == pt_ass) continue;
	  
	  double effweight = effweight_trg * effweight_ass * nMult_corr_trg;
	  double absDelEta = TMath::Abs(deltaEta);
	  //std::cout<<"charge triger=" << chg_trg <<std::endl;
	  
	  //std::cout<<"deltaEta value=" << deltaEta <<std::endl;
	  //Fill and symmetrize the distribution	  
	  if( chg_trg > 0 && chg_ass > 0){
	    if(isGen){
	      hSignalPPGen->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalPPGen->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalPPGen->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight);
	      hSignalPPGen->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight ); 
	    }
	    else{
	      hSignalPP->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalPP->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalPP->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight);
	      hSignalPP->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight );
	    }
	  }//+ve, +ve particle only
	  else if( chg_trg > 0 && chg_ass < 0){
	    if(isGen){
	      hSignalPMGen->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalPMGen->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalPMGen->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight);
	      hSignalPMGen->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight );
	    }
	    else{
	      hSignalPM->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalPM->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalPM->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight);
	      hSignalPM->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight );
	    }
	  }//+ve. -ve particle only
	  else if( chg_trg < 0 && chg_ass > 0){
	    if(isGen){
	      hSignalMPGen->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalMPGen->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalMPGen->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight);
	      hSignalMPGen->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight );
	    }
	    else
	      {
	    hSignalMP->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
	    hSignalMP->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	    hSignalMP->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight);
	    hSignalMP->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight );
	      }
	  }//-ve, +ve particle only
	  else if(chg_trg < 0 && chg_ass < 0){
	    if(isGen){
	      hSignalMMGen->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalMMGen->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalMMGen->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight);
	      hSignalMMGen->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight );
	    }
	    else{
	      hSignalMM->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalMM->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hSignalMM->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight);
	      hSignalMM->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight );
	    }
	  }//-ve, -ve particle only
	  
	}//nass----
       
    }//ntrg----
}

//======================================================================================
void CmsTwoPartCorrAnalysis::FillHistsBackground(int ievt_trg, int jevt_ass, bool isGen)
{
  
  if( evtVec[ievt_trg].run   ==  evtVec[jevt_ass].run &&
      evtVec[ievt_trg].event == evtVec[jevt_ass].event )
  {
      std::cout << "Event are the same. Skipping it" << std::endl;
      return;
  }

  int kk;
  isGen ? kk = 0 : kk = 1;
  
  unsigned int ntrgsize = evtVec[ievt_trg].pVect_trg[kk].size();
  unsigned int nasssize = evtVec[jevt_ass].pVect_ass[kk].size();
  double nMult_corr_trg = evtVec[ievt_trg].nMultCorrVect_trg[kk];
  //std::cout<<" ntrgsize value=" << ntrgsize <<std::endl;
  for( unsigned int itrg = 0; itrg < ntrgsize; itrg++ )
    {
      
      TLorentzVector pvector_trg = (evtVec[ievt_trg].pVect_trg[kk])[itrg];
      double effweight_trg       = (evtVec[ievt_trg].effVect_trg[kk])[itrg];
      int chg_trg                = (evtVec[ievt_trg].chgVect_trg[kk])[itrg];
      double eta_trg             = pvector_trg.Eta();
      double phi_trg             = pvector_trg.Phi();
      double pt_trg              = pvector_trg.Pt();
      
      for( unsigned int jass = 0; jass < nasssize; jass++ )
	{
	  TLorentzVector pvector_ass = (evtVec[jevt_ass].pVect_ass[kk])[jass];
	  double effweight_ass       = (evtVec[jevt_ass].effVect_ass[kk])[jass];
	  int chg_ass                = (evtVec[jevt_ass].chgVect_ass[kk])[jass];
	  double eta_ass             = pvector_ass.Eta();
	  double phi_ass             = pvector_ass.Phi();
	  double pt_ass              = pvector_ass.Pt();
	  
	  double deltaPhi            = GetDeltaPhi( phi_trg, phi_ass );
	  double deltaPhi2           = GetDeltaPhi( phi_ass, phi_trg);
	  double deltaEta            = GetDeltaEta( eta_trg, eta_ass );
	  
	  if( deltaEta == 0.0  &&  deltaPhi == 0.0 && deltaPhi2 == 0.0  && pt_trg == pt_ass) continue;
	
	  //Total weight
	  double effweight = effweight_trg * effweight_ass* nMult_corr_trg;
	  double absDelEta = TMath::Abs(deltaEta);
	  
	  if( chg_trg > 0 && chg_ass > 0)
	    {
	      if(isGen){
		hBackgroundPPGen->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
		hBackgroundPPGen->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
		hBackgroundPPGen->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight );
		hBackgroundPPGen->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight ); 
	      }
	      else{
		hBackgroundPP->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
		hBackgroundPP->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
		hBackgroundPP->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight );
		hBackgroundPP->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight ); 
	      }
	    }
	  else if( chg_trg > 0 && chg_ass < 0){
	    if(isGen){
	      hBackgroundPMGen->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight  );
	      hBackgroundPMGen->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hBackgroundPMGen->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight );
	      hBackgroundPMGen->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight ); 
	    }
	    else{
	      hBackgroundPM->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight  );
	      hBackgroundPM->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hBackgroundPM->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight );
	      hBackgroundPM->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight ); 
	    }
	  }
	  else if( chg_trg < 0 && chg_ass > 0){
	    if(isGen){
	      hBackgroundMPGen->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hBackgroundMPGen->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hBackgroundMPGen->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight );
	      hBackgroundMPGen->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight ); 
	    }
	    else{
	      hBackgroundMP->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hBackgroundMP->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hBackgroundMP->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight );
	      hBackgroundMP->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight ); 
	    }
	  }
	  else if( chg_trg < 0 && chg_ass < 0){
	    if(isGen){
	      hBackgroundMMGen->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hBackgroundMMGen->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hBackgroundMMGen->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight );
	      hBackgroundMMGen->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight ); 
	    }
	    else{
	      hBackgroundMM->Fill( absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hBackgroundMM->Fill( -1.*absDelEta, deltaPhi, 1.0/4.0/effweight );
	      hBackgroundMM->Fill( absDelEta, deltaPhi2, 1.0/4.0/effweight );
	      hBackgroundMM->Fill( -1.*absDelEta, deltaPhi2, 1.0/4.0/effweight ); 
	    }
	  }
	}//jass--
    }//trg
}
//-------------------------------------------
double CmsTwoPartCorrAnalysis::getEff(const TLorentzVector pvector, int evtclass)
{
   double effweight = 1.0;
   for(int i_cf=0; i_cf < static_cast<int>(effCorrBinMin.size()); i_cf++) if(evtclass >= effCorrBinMin[i_cf] && evtclass < effCorrBinMax[i_cf]) {
       effweight = heff[i_cf]->GetBinContent(
					      heff[i_cf]->GetXaxis()->FindBin(pvector.Eta()),
					      heff[i_cf]->GetYaxis()->FindBin(pvector.Pt())
					      );
     }
   return effweight;
   }

 //------------------------------------------------

  
//}
//}

//define this as a plug-in
DEFINE_FWK_MODULE(CmsTwoPartCorrAnalysis);
