#if !defined(__CINT__) && !defined(__CLING__)
#include <iostream>
#include <vector>

#include <TNamed.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>

#include <FairMCTrack.h>
#include <MpdPid.h>
#include "MpdMCEventHeader.h"
#include "MpdEvent.h"
#include "MpdTrack.h"
#include <MpdTpcKalmanTrack.h>
#include <MpdKalmanTrack.h>
#include <MpdVertex.h>
#include <MpdTofMatching.h>
#endif
char name[1000];
void spectra_from_dst_5(const Int_t events2analyze = 250) {
  

  //    TString inname          =    "urqmd-BiBi-09.0GeV-0-14fm-eos0-250-2498-3.reco.root";

    TString outEffFile      =    "testKalmanEfficiencyData_9GEV.root";
    TString outdEdxFile     =    "testKalmandEdxData_9GEV.root";
    TString outSpectraFile  =    "testKalmanSpectraData_9GEV.root";

    Double_t fProbCut = 0.95;

    Double_t sigM = 2.0, sigE = 2.0, energy = 9., koef = 1.; /// n-sigma bands for PID selection
    TString Generator = "URQMD", Tracking = "CF";
    MpdPid *pid = new MpdPid(sigM, sigE, energy, koef, Generator, Tracking, "pikapr");

    Int_t nbinsEffp = 30, nbinsEffpt = 25, nbinsSpectra = 25, nbinsdEdx = 1000;

    Float_t minpt = 0.;
    Float_t maxpt = 2.5;
    Float_t maxp = 3.;

    TH1F *mcpplus6   = new TH1F("mcpplus6",   "mcpplus dedx and m2",       nbinsEffp,minpt,maxp);
    TH1F *mcpminus6  = new TH1F("mcpminus6",  "mcpminus dedx and m2",      nbinsEffp,minpt,maxp);
    TH1F *mcpiplus6  = new TH1F("mcpiplus6",  "mcpiplus dedx and m2",      nbinsEffp,minpt,maxp);
    TH1F *mcpiminus6 = new TH1F("mcpiminus6", "mcpiminus dedx and m2",     nbinsEffp,minpt,maxp);
    TH1F *mckplus6   = new TH1F("mckplus6",   "mckplus dedx and m2",       nbinsEffp,minpt,maxp);
    TH1F *mckminus6  = new TH1F("mckminus6",  "mckminus dedx and m2",      nbinsEffp,minpt,maxp);
    TH1F *pplus6     = new TH1F("pplus6",     "pplus dedx and m2",         nbinsEffp,minpt,maxp);
    TH1F *pminus6    = new TH1F("pminus6",    "pminus dedx and m2",        nbinsEffp,minpt,maxp);
    TH1F *piplus6    = new TH1F("piplus6",    "piplus dedx and m2",        nbinsEffp,minpt,maxp);
    TH1F *piminus6   = new TH1F("piminus6",   "piminus dedx and m2",       nbinsEffp,minpt,maxp);
    TH1F *kplus6     = new TH1F("kplus6",     "kplus dedx and m2",         nbinsEffp,minpt,maxp);
    TH1F *kminus6    = new TH1F("kminus6",    "kminus dedx and m2",        nbinsEffp,minpt,maxp);
    TH1F *tpplus6    = new TH1F("tpplus6",    "tptrueplus dedx and m2",    nbinsEffp,minpt,maxp);
    TH1F *tpminus6   = new TH1F("tpminus6",   "tptrueminus dedx and m2",   nbinsEffp,minpt,maxp);
    TH1F *tpiplus6   = new TH1F("tpiplus6",   "tpitrueplus dedx and m2",   nbinsEffp,minpt,maxp);
    TH1F *tpiminus6  = new TH1F("tpiminus6",  "tpitrueminus dedx and m2",  nbinsEffp,minpt,maxp);
    TH1F *tkplus6    = new TH1F("tkplus6",    "tktrueplus dedx and m2",    nbinsEffp,minpt,maxp);
    TH1F *tkminus6   = new TH1F("tkminus6",   "tktrueminus dedx and m2",   nbinsEffp,minpt,maxp);
    TH1F *fpplus6    = new TH1F("fpplus6",    "fpfalseplus dedx and m2",   nbinsEffp,minpt,maxp);
    TH1F *fpminus6   = new TH1F("fpminus6",   "fpfalseminus dedx and m2",  nbinsEffp,minpt,maxp);
    TH1F *fpiplus6   = new TH1F("fpiplus6",   "fpifalseplus dedx and m2",  nbinsEffp,minpt,maxp);
    TH1F *fpiminus6  = new TH1F("fpiminus6",  "fpifalseminus dedx and m2", nbinsEffp,minpt,maxp);
    TH1F *fkplus6    = new TH1F("fkplus6",    "fkfalseplus dedx and m2",   nbinsEffp,minpt,maxp);
    TH1F *fkminus6   = new TH1F("fkminus6",   "fkfalseminus dedx and m2",  nbinsEffp,minpt,maxp);


    TH1F *mcpplus6pt  = new TH1F("mcpplus6pt",  "mcpplus dedx and m2",      nbinsEffpt,minpt,maxpt);
    TH1F *mcpiplus6pt = new TH1F("mcpiplus6pt", "mcpiplus dedx and m2",     nbinsEffpt,minpt,maxpt);
    TH1F *mckplus6pt  = new TH1F("mckplus6pt",  "mckplus dedx and m2",      nbinsEffpt,minpt,maxpt);
    TH1F *pplus6pt    = new TH1F("pplus6pt",    "pplus dedx and m2",        nbinsEffpt,minpt,maxpt);
    TH1F *piplus6pt   = new TH1F("piplus6pt",   "piplus dedx and m2",       nbinsEffpt,minpt,maxpt);
    TH1F *kplus6pt    = new TH1F("kplus6pt",    "kplus dedx and m2",        nbinsEffpt,minpt,maxpt);
    TH1F *tpplus6pt   = new TH1F("tpplus6pt",   "tptrueplus dedx and m2",   nbinsEffpt,minpt,maxpt);
    TH1F *tpiplus6pt  = new TH1F("tpiplus6pt",  "tpitrueplus dedx and m2",  nbinsEffpt,minpt,maxpt);
    TH1F *tkplus6pt   = new TH1F("tkplus6pt",   "tktrueplus dedx and m2",   nbinsEffpt,minpt,maxpt);
    TH1F *fpplus6pt   = new TH1F("fpplus6pt",   "fpfalseplus dedx and m2",  nbinsEffpt,minpt,maxpt);
    TH1F *fpiplus6pt  = new TH1F("fpiplus6pt",  "fpifalseplus dedx and m2", nbinsEffpt,minpt,maxpt);
    TH1F *fkplus6pt   = new TH1F("fkplus6pt",   "fkfalseplus dedx and m2",  nbinsEffpt,minpt,maxpt);

    TH1F *mcpminus6pt  = new TH1F("mcpminus6pt",   "mcpplus dedx and m2",       nbinsEffpt,minpt,maxpt);
    TH1F *mcpiminus6pt = new TH1F("mcpiminus6pt",  "mcpiplus dedx and m2",      nbinsEffpt,minpt,maxpt);
    TH1F *mckminus6pt  = new TH1F("mckminus6pt",   "mckplus dedx and m2",       nbinsEffpt,minpt,maxpt);
    TH1F *pminus6pt    = new TH1F("pminus6pt",     "pplus dedx and m2",         nbinsEffpt,minpt,maxpt);
    TH1F *piminus6pt   = new TH1F("piminus6pt",    "piplus dedx and m2",        nbinsEffpt,minpt,maxpt);
    TH1F *kminus6pt    = new TH1F("kminus6pt",     "kplus dedx and m2",         nbinsEffpt,minpt,maxpt);
    TH1F *tpminus6pt   = new TH1F("tpminus6pt",    "tptrueplus dedx and m2",    nbinsEffpt,minpt,maxpt);
    TH1F *tpiminus6pt  = new TH1F("tpiminus6pt",   "tpitrueplus dedx and m2",   nbinsEffpt,minpt,maxpt);
    TH1F *tkminus6pt   = new TH1F("tkminus6pt",    "tktrueplus dedx and m2",    nbinsEffpt,minpt,maxpt);
    TH1F *fpminus6pt   = new TH1F("fpminus6pt",    "fpfalseplus dedx and m2",   nbinsEffpt,minpt,maxpt);
    TH1F *fpiminus6pt  = new TH1F("fpiminus6pt",   "fpifalseplus dedx and m2",  nbinsEffpt,minpt,maxpt);
    TH1F *fkminus6pt   = new TH1F("fkminus6pt",    "fkfalseplus dedx and m2",   nbinsEffpt,minpt,maxpt);


    TH1F *hRefMult = new TH1F("hRefMult", "Reference multiplicity; refMult", 350, -0.5, 349.5);


    TH2F *dEdXfrommpdPosAll  = new TH2F("dEdXfrommpdPosAll",   "dEdX for all pos charge", nbinsSpectra, minpt, maxpt, nbinsdEdx, 0, 40.e04);
    TH2F *dEdXfrommpdNegAll  = new TH2F("dEdXfrommpdNegAll",   "dEdX for all neg charge", nbinsSpectra, minpt, maxpt, nbinsdEdx, 0, 40.e04);

    TH2F *dEdXfrommpdPosP    = new TH2F("dEdXfrommpdPosP",     "dEdX for all pos charge", nbinsSpectra, minpt, maxpt, nbinsdEdx, 0, 5.e04);
    TH2F *dEdXfrommpdPosPi   = new TH2F("dEdXfrommpdPosPi",    "dEdX for all pos charge", nbinsSpectra, minpt, maxpt, nbinsdEdx, 0, 5.e04);
    TH2F *dEdXfrommpdPosKa   = new TH2F("dEdXfrommpdPosKa",    "dEdX for all pos charge", nbinsSpectra, minpt, maxpt, nbinsdEdx, 0, 5.e04);
    TH2F *dEdXfrommpdNegP    = new TH2F("dEdXfrommpdNegP",     "dEdX for all pos charge", nbinsSpectra, minpt, maxpt, nbinsdEdx, 0, 5.e04);
    TH2F *dEdXfrommpdNegPi   = new TH2F("dEdXfrommpdNegPi",    "dEdX for all pos charge", nbinsSpectra, minpt, maxpt, nbinsdEdx, 0, 5.e04);
    TH2F *dEdXfrommpdNegKa   = new TH2F("dEdXfrommpdNegKa",    "dEdX for all pos charge", nbinsSpectra, minpt, maxpt, nbinsdEdx, 0, 5.e04);


     TH2F* MultPiPlus     = new TH2F("MultPiPlus",    "Pion pos momentum spectra",   350, -0.5, 349.5, nbinsSpectra, minpt, maxpt);
     TH2F* MultPiMinus    = new TH2F("MultPiMinus",   "Pion neg momentum spectra",   350, -0.5, 349.5, nbinsSpectra, minpt, maxpt);
     TH2F* MultKaPlus     = new TH2F("MultKaPlus",    "Kaon pos momentum spectra",   350, -0.5, 349.5, nbinsSpectra, minpt, maxpt);
     TH2F* MultKaMinus    = new TH2F("MultKaMinus",   "Kaon neg momentum spectra",   350, -0.5, 349.5, nbinsSpectra, minpt, maxpt);
     TH2F* MultPPlus      = new TH2F("MultPPlus",     "Proton pos momentum spectra", 350, -0.5, 349.5, nbinsSpectra, minpt, maxpt);
     TH2F* MultPMinus     = new TH2F("MultPMinus",    "Proton neg momentum spectra", 350, -0.5, 349.5, nbinsSpectra, minpt, maxpt);


    Int_t mcNtracks, mcTrackId, nTofMatch, pdg, pdgmc, mother, charge, fNtracks, ID, tofFlag, maxloca, nentries, noprobcut=0;

    Double_t px, py ,pz, mpdPt, mpdP, mcPt, dedx, m2, Eta, AbsEta, pion, kaon, proton, electron, maxprob, NumChargeTracks=0;
    Bool_t ret;


    for (Int_t nf = 2497; nf <2499; nf++){ 
      for (Int_t k=1; k<4; k++){
	sprintf(name,"urqmd-BiBi-09.0GeV-0-14fm-eos0-250-%d-%d.reco.root",nf,k);
	
	TString inname=name;
	cout <<" Openning file " << name << endl;
	TFile fileInput(inname.Data());
	if (fileInput.IsZombie()) continue;
	
	TChain chain("mpdsim");
	
	//    TChain chain("cbmsim");
	chain.Add(inname);
	
	TClonesArray *mcTracks = NULL;
	TClonesArray *mpdTracks = NULL;
	FairMCEventHeader *fmcHeader = NULL;
	TClonesArray *mpdKalmanTracks = NULL;
	TClonesArray *tofMatches = NULL;
	
	chain.SetBranchAddress("TOFMatching", &tofMatches);
	chain.SetBranchAddress("MCEventHeader.", &fmcHeader);
	chain.SetBranchAddress("MCTrack", &mcTracks);
	
	MpdEvent *event = new MpdEvent();
	chain.SetBranchAddress("MPDEvent.", &event);
	
	FairMCTrack *mctrack = new FairMCTrack;
	MpdTrack *mpdtrack = new MpdTrack;
	
	
	//    KalmanTracks = (TClonesArray*) TFile(inname).FindObjectAny("TpcKalmanTrack");
	chain.SetBranchAddress("TpcKalmanTrack", &mpdKalmanTracks);
	
	
	MpdKalmanTrack *KalmanTrack = new MpdKalmanTrack;
	
	MpdTofMatchingData *match = new MpdTofMatchingData;
	map<Int_t,Int_t> mapTof;
	
	
	if (events2analyze == 0) {
	  nentries = chain.GetEntries(); // get all entries
	  //cout<<nentries<<endl;
	}
	else {
	  nentries = events2analyze;
	}
	
	
	
	for (Int_t i=0; i<nentries; i++)  // entries loop
	  {
	    
	    chain.GetEntry(i); // current entry
	    
	    mcNtracks = mcTracks->GetEntriesFast();
	    
	    mpdTracks = (TClonesArray*) event->GetGlobalTracks();
	    fNtracks = mpdTracks->GetEntriesFast();
	    
	    nTofMatch = tofMatches->GetEntriesFast(); // making a ToF match map
	    for (Int_t itof = 0; itof < nTofMatch; ++itof) {
	      match = (MpdTofMatchingData*) tofMatches->UncheckedAt(itof);
	      mapTof[match->GetKFTrackIndex()] = itof;
	    }
	    
	    
	    // Int_t NumTracks = (Int_t) event->GetEventInfoNofGlobalTracks();
	    
	    
	    // cout<<mcNtracks<<endl;
	    // cout<<fNtracks<<endl;
	    
	    cout << i+1 << "/" << nentries << " event is processed... \r" << flush;
	    
	    //---------------------- Apply event cuts here
	    
	    //if (fmcHeader->GetB() > 14.0) continue;
	    
	    
	    //----------------------------------------
	    
	    
	    for(Int_t iTrack = 0; iTrack < fNtracks; iTrack++) //tracks loop
	      {
		
		mpdtrack = (MpdTrack*) mpdTracks->UncheckedAt(iTrack);
		KalmanTrack = (MpdKalmanTrack*) mpdKalmanTracks->UncheckedAt(iTrack);
		
		
		ID = KalmanTrack->GetTrackID();
		// cout << iTrack << " track is processed... \r";
		// cout<<endl;
		// cout<<"momentum      "<<mctrack->GetPt()<<endl;
		
		mctrack = (FairMCTrack*) mcTracks->UncheckedAt(ID);
		
		
		mcPt = mctrack->GetPt();
		px = KalmanTrack->Momentum3().Px(); py = KalmanTrack->Momentum3().Py(); pz = KalmanTrack->Momentum3().Pz();
		mpdP = TMath::Sqrt(px*px + py*py + pz*pz);
		
		//cout<<"px = "<<px<<"  py = "<<py<<"  pz = "<<pz<<endl;
		
		
		Eta = KalmanTrack->Momentum3().Eta();
		AbsEta = TMath::Abs(Eta);
		
		mpdPt = KalmanTrack->Momentum3().Pt();
		charge = KalmanTrack->Charge(); 
		//charge = 1; if (mpdPt > 0.) charge = -1;
		mpdPt = TMath::Abs(mpdPt);
		
		pdgmc = mctrack->GetPdgCode();
		mother = mctrack->GetMotherId();
		
		// dedx
		dedx = KalmanTrack->GetPartID();
		
		//tofFlag = mpdtrack->GetTofFlag();
		
		ret = kFALSE;
		//	m2 = mpdtrack->GetTofMass2();
		
		
		
		
		//---------------------- Apply track cuts here
		
		if (mother != -1) continue;
		if (KalmanTrack->GetNofHits() < 20) continue;
		if (AbsEta > 0.5) continue;
		if ((mpdPt < 0.1) || (mpdPt > 2.5)) continue;
		
		
		//----------------------------------------
		
		if (charge != 0) { NumChargeTracks++;  }

		switch(pdgmc) // Fill mc histograms
		  {
		  case 211:  mcpiplus6->Fill(mpdP); mcpiplus6pt->Fill(mpdPt); break;
		  case 321:  mckplus6->Fill(mpdP); mckplus6pt->Fill(mpdPt); break;
		  case 2212: mcpplus6->Fill(mpdP); mcpplus6pt->Fill(mpdPt); break;
		  case -211:  mcpiminus6->Fill(mpdP); mcpiminus6pt->Fill(mpdPt); break;
		  case -321:  mckminus6->Fill(mpdP); mckminus6pt->Fill(mpdPt); break;
		  case -2212: mcpminus6->Fill(mpdP); mcpminus6pt->Fill(mpdPt); break;
		  default: break;
		  }
		
		
		//if ( (tofFlag == 2) || (tofFlag == 6) ) // change to tof matching map
		if (mapTof.count(iTrack) > 0)
		  {
		    m2 = ((MpdTofMatchingData*)tofMatches->UncheckedAt(mapTof[iTrack]))->GetMass2();
		    
		    ret = pid->FillProbs(mpdP, dedx, m2, charge); /// dE/dx+TOF
		    if ((!ret) && (mpdP < 0.8)) ret = pid->FillProbs(mpdP, dedx, charge); /// dE/dx only
		  }
		else
		  {
		    if (mpdP < 0.8) ret = pid->FillProbs(mpdP, dedx, charge); /// dE/dx only
		    
		  }
		
		if (!ret) continue;
		
		pion = pid->GetProbPi();
		proton = pid->GetProbPr();
		kaon = pid->GetProbKa();
		//electron = pid->GetProbEl();   //electrons are not defined in current version of MPDPID
		maxloca = 0;
		Double_t Probs[] = { pion, kaon, proton};
		maxprob = TMath::MaxElement(3, Probs);
		maxloca = TMath::LocMax(3, Probs);
		pdg = 0;
		// maxloca = 0;
		
		if (maxprob<fProbCut) {noprobcut++; continue;}
		if (charge > 0)
		  {
		    dEdXfrommpdPosAll->Fill(mpdP, dedx);
		    //dEdXfrommpdPosAlleta->Fill(mpdP, dedx, Eta);
		    switch(maxloca)
		      {
		      case 0: pdg=211; piplus6->Fill(mpdP); piplus6pt->Fill(mpdPt); dEdXfrommpdPosPi->Fill(mpdP, dedx); break;
		      case 1: pdg=321; kplus6->Fill(mpdP); kplus6pt->Fill(mpdPt); dEdXfrommpdPosKa->Fill(mpdP, dedx); break;
		      case 2: pdg=2212; pplus6->Fill(mpdP); pplus6pt->Fill(mpdPt); dEdXfrommpdPosP->Fill(mpdP, dedx); break;
		      default: pdg=0;
		      }
		  }
		
		if (charge < 0)
		  {
		    dEdXfrommpdNegAll->Fill(mpdP, dedx);
		    //dEdXfrommpdNegAlleta->Fill(mpdP, dedx, Eta);
		    switch(maxloca)
		      {
		      case 0: pdg=-211; piminus6->Fill(mpdP); piminus6pt->Fill(mpdPt); dEdXfrommpdNegPi->Fill(mpdP, dedx); break;
		      case 1: pdg=-321; kminus6->Fill(mpdP); kminus6pt->Fill(mpdPt); dEdXfrommpdNegKa->Fill(mpdP, dedx); break;
		      case 2: pdg=-2212; pminus6->Fill(mpdP); pminus6pt->Fill(mpdPt); dEdXfrommpdNegP->Fill(mpdP, dedx); break;
		      default: pdg=0;
		      }
		  }
		
		
		
		if (charge > 0)
		  { // Fill true and false for charge > 0
		    
		    if(pdgmc == pdg) //pdgmc == pdg fill true identified
		      {
			switch(maxloca)
			  {
			  case 0: tpiplus6->Fill(mpdP); tpiplus6pt->Fill(mpdPt); break;
			  case 1: tkplus6->Fill(mpdP); tkplus6pt->Fill(mpdPt); break;
			  case 2: tpplus6->Fill(mpdP); tpplus6pt->Fill(mpdPt); break;
			  default: break;
			  }
		      }
		    else ///pdgmc != pdg fill false identified
		      {
			switch(maxloca)
			  {
			  case 0: fpiplus6->Fill(mpdP); fpiplus6pt->Fill(mpdPt); break;
			  case 1: fkplus6->Fill(mpdP); fkplus6pt->Fill(mpdPt); break;
			  case 2: fpplus6->Fill(mpdP); fpplus6pt->Fill(mpdPt); break;
			  default: break;
			  }
		      }
		  } // end of charge = 1
		
		
		
		if (charge < 0)
		  { // Fill true and false for charge < 0
		    if(pdgmc == pdg) //pdgmc == pdg fill true identified
		      {
			switch(maxloca)
			  {
			  case 0: tpiminus6->Fill(mpdP); tpiminus6pt->Fill(mpdPt); break;
			  case 1: tkminus6->Fill(mpdP); tkminus6pt->Fill(mpdPt); break;
			  case 2: tpminus6->Fill(mpdP); tpminus6pt->Fill(mpdPt); break;
			  default: break;
			  }
		      }
		    else ///pdgmc != pdg fill falsly identified
		      {
			switch(maxloca)
			  {
			  case 0: fpiminus6->Fill(mpdP); fpiminus6pt->Fill(mpdPt); break;
			  case 1: fkminus6->Fill(mpdP); fkminus6pt->Fill(mpdPt);  break;
			  case 2: fpminus6->Fill(mpdP); fpminus6pt->Fill(mpdPt);  break;
			  default: break;
			  }
		      }
		  } // end of charge = -1

		if (charge > 0)
		  {
		    switch(maxloca)
		      {
		      case 0: pdg=211;   MultPiPlus->Fill(NumChargeTracks,mpdPt); break;
		      case 1: pdg=321;   MultKaPlus->Fill(NumChargeTracks,mpdPt);break;
		      case 2: pdg=2212;  MultPPlus->Fill(NumChargeTracks,mpdPt); break;
		      default: pdg=0;
		      }
		  } // end of charge = 1
		if (charge < 0)
		  {
		    switch(maxloca)
		      {
		      case 0: pdg=-211;   MultPiMinus->Fill(NumChargeTracks,mpdPt); break;
		      case 1: pdg=-321;   MultKaMinus->Fill(NumChargeTracks,mpdPt);break;
		      case 2: pdg=-2212;  MultPMinus->Fill(NumChargeTracks,mpdPt); break; 
		      default: pdg=0;
		      }
		  } // end of charge = -1
	      } // end of globalTracks

	    hRefMult->Fill(NumChargeTracks);
	    // NumChargeTracks = 0;
	    
	    
	    // 
	    // 
	    
	    
	    /*   for(Int_t iTrack = 0; iTrack < fNtracks; iTrack++) //tracks loop
		 {
		 
		 mpdtrack = (MpdTrack*) mpdTracks->UncheckedAt(iTrack);
		 KalmanTrack = (MpdKalmanTrack*) mpdKalmanTracks->UncheckedAt(iTrack);
		 
		 ID = KalmanTrack->GetID();
		 mctrack = (FairMCTrack*) mcTracks->UncheckedAt(ID);
		 
		 
		 mcPt = mctrack->GetPt();
		 px = KalmanTrack->GetPx(); py = KalmanTrack->GetPy(); pz = KalmanTrack->GetPz();
		 mpdP = TMath::Sqrt(px*px + py*py + pz*pz);
		 
		 Eta = KalmanTrack->GetEta();
		 AbsEta = TMath::Abs(Eta);
		 
		 mpdPt = KalmanTrack->GetPt();
		 charge = 1; if (mpdPt > 0.) charge = -1;
		 mpdPt = TMath::Abs(mpdPt);
		 
		 pdgmc = mctrack->GetPdgCode();
		 mother = mctrack->GetMotherId();
		 
		 tofFlag = KalmanTrack->GetTofFlag();
		 
		 ret = kFALSE;
		 m2 = 0.;
		 
		 
		 //---------------------- Apply track cuts here
		 
		 if (mother != -1) continue;
		 if (KalmanTrack->GetNofHits() < 25) continue;
		 if (AbsEta > 0.5) continue;
		 if ((mpdPt < 0.1) || (mpdPt > 2.)) continue;
		 
		 
		 
		 //----------------------------------------
		 
		 
		 
		 if ( (tofFlag == 2) || (tofFlag == 6) )
		 {
		 
		 ret = pid->FillProbs(mpdP, dedx, m2, charge); /// dE/dx+TOF
		 if ((!ret) && (mpdP < 0.8)) ret = pid->FillProbs(mpdP, dedx, charge); /// dE/dx only
		 }
		 else
		 {
		 if (mpdP < 0.8) ret = pid->FillProbs(mpdP, dedx, charge); /// dE/dx only
		 
		 }
		 
		 if (!ret) continue;
		 
		 
		 pion = pid->GetProbPi();
		 proton = pid->GetProbPr();
		 kaon = pid->GetProbKa();
		 maxloca = 0;
		 Double_t Probs[] = { pion, kaon, proton};
		 maxprob = TMath::MaxElement(3, Probs);
		 maxloca = TMath::LocMax(3, Probs);
		 pdg = 0;
		 Double_t weight = 0;
	  */
	    NumChargeTracks = 0;
	    mapTof.clear();
	    
	  } // end of events
      } // k loop files
    } // nf loop files
    cout<<endl;
    
    
    
    
    //Save efficiency plots
    
    TFile *fout = new TFile(outEffFile, "RECREATE");
    
    mcpplus6pt->Write();
    mcpiplus6pt->Write();
    mckplus6pt->Write();
    pplus6pt->Write();
    piplus6pt->Write();
    kplus6pt->Write();
    tpplus6pt->Write();
    tpiplus6pt->Write();
    tkplus6pt->Write();
    fpplus6pt->Write();
    fpiplus6pt->Write();
    fkplus6pt->Write();
    
    mcpminus6pt->Write();
    mcpiminus6pt->Write();
    mckminus6pt->Write();
    pminus6pt->Write();
    piminus6pt->Write();
    kminus6pt->Write();
    tpminus6pt->Write();
    tpiminus6pt->Write();
    tkminus6pt->Write();
    fpminus6pt->Write();
    fpiminus6pt->Write();
    fkminus6pt->Write();
    
    mcpplus6->Write();
    mcpminus6->Write();
    mcpiplus6->Write();
    mcpiminus6->Write();
    mckplus6->Write();
    mckminus6->Write();
    pplus6->Write();
    pminus6->Write();
    piplus6->Write();
    piminus6->Write();
    kplus6->Write();
    kminus6->Write();
    tpplus6->Write();
    tpminus6->Write();
    tpiplus6->Write();
    tpiminus6->Write();
    tkplus6->Write();
    tkminus6->Write();
    fpplus6->Write();
    fpminus6->Write();
    fpiplus6->Write();
    fpiminus6->Write();
    fkplus6->Write();
    fkminus6->Write();
    
    fout->Close();
    
    
    //Save energy loss plots
    
    TFile *fFile = new TFile(outdEdxFile, "RECREATE");
    
    dEdXfrommpdPosAll->Write();
    dEdXfrommpdNegAll->Write();	
    dEdXfrommpdPosP->Write();
    dEdXfrommpdPosPi->Write();
    dEdXfrommpdPosKa->Write();
    dEdXfrommpdNegP->Write();
    dEdXfrommpdNegPi->Write();
    dEdXfrommpdNegKa->Write();
    
    fFile->Close();
    
    //Save particle spectra and refMult plots
    
    TFile *outFile = new TFile(outSpectraFile, "RECREATE");
    hRefMult->Write();
    MultPiPlus->Write();
    MultKaPlus->Write();
    MultPPlus->Write();
    MultPiMinus->Write();
    MultKaMinus->Write();
    MultPMinus->Write();
    
    outFile->Close();
    
    
    TCanvas *c1 = new TCanvas();
    c1->Divide(2,2);
    c1->cd(1);
    hRefMult->Draw();
    
    c1->cd(2);
    MultPiPlus->Draw();
    
    c1->cd(3);
    MultPPlus->Draw();
    
    c1->cd(4);
    MultKaPlus->Draw();
    
    
    cout << "End of spectra_from_dst" << endl;
    
}
