#include "EMCALChannelTiming.h"

EMCALChannelTiming::EMCALChannelTiming(int seegm, const std::string &name):
SubsysReco(name)
{
	//constructor
	int nbins = 1000; //maybe put this in the initial values???
	float Eu = 15.;
	float El = 1e-3;	
	seg=seegm;
	float log_binsize = (std::log(Eu) - std::log(El))/((float)nbins-1); //bins go linear in log(E)
	for(int i=0; i==nbins; i++)
	{
		float bin_log = std::log(El) + i*log_binsize;
		energy_bins.push_back(bin_log);
	}
	energy_bins.push_back(std::log(Eu)); 
	towerTree = new TTree("tT", "tT");
	for(int i = 0; i < 15*1536; i++)
	{
		tower* tw=new tower();
		allE->push_back(tw);
		towerTree->Branch(std::format("tower_{}", i).c_str(), &tw);
	}
	
}
int EMCALChannelTiming::Init(PHCompositeNode *topNode)
{
	int ntowers=1536*16;
	AddBins( 
		AllTowers1D, AllTowers2D, 
		highTowers1D, highTowers2D,
		lowTowers1D, lowTowers2D
	       );
	for(int i=0; i<ntowers; i++)
	{
		AddBins( 
			AllTowers1D_t->at(i), AllTowers2D_t->at(i), 
			highTowers1D_t->at(i), highTowers2D_t->at(i),
			lowTowers1D_t->at(i), lowTowers2D_t->at(i), std::to_string(i)
	       	);
	}
	auto evt = findNode::getClass<EventHeader>(topNode, "EventHeader");
	int run = evt->get_RunNumber();
	output_file_name=std::format("EMCAL_timing_run-{}_segment-{}.root", run, seg); 
	return Fun4AllReturnCodes::EVENT_OK;


}
void EMCALChannelTiming::AddBins(
		std::vector<TH1F*>* all1DHists, 
		std::vector<TH2F*>* all2DHists,
		std::vector<TH1F*>* high1DHists, 
		std::vector<TH2F*>* high2DHists,
		std::vector<TH1F*>* low1DHists, 
		std::vector<TH2F*>* low2DHists,
		std::string ntower/*=""*/ 
		)

{
	//per tower energy variation
	std::string ntowerus	= (ntower=="") ? "" : "_" + ntower;
	ntower			= (ntower=="") ? "" : " " + ntower;
	int nbins = 1000; //maybe put this in the initial values???
	TH1F* DeltaT = new TH1F(
			std::format("Delta_T{}", ntowerus).c_str(), 
			std::format("#Delta T tower{}; #Delta T [ns]", ntower).c_str(), 
			100, -20, 20);
	TH1F* DeltaT_high = new TH1F(
			std::format("Delta_T_high{}", ntowerus).c_str(), 
			std::format("#Delta T tower{}; #Delta T [ns]", ntower).c_str(), 
			100, -20, 20);
	TH1F* DeltaT_low = new TH1F(
			std::format("Delta_T_low{}", ntowerus).c_str(), 
			std::format("#Delta T tower{}; #Delta T [ns]", ntower).c_str(), 
			100, -20, 20);

	//Energy Hists
	TH1F* Energy = new TH1F(
			
			"Tow_E", "EMCAL Tower Energy; log(E) [GeV]; N_{tow}", 
			nbins, energy_bins.data()); 
	//1 MeV -15 GeV, to keep comparable between low and high energy towers 
	TH1F* Energy_high = new TH1F(
			"Tow_E_high", "EMCAL Tower Energy; log(E) [GeV]; N_{tow}", 
			nbins, energy_bins.data()); 
	//1 MeV -15 GeV, to keep comparable between low and high energy towers 
	TH1F* Energy_low = new TH1F(
			"Tow_E_low", "EMCAL Tower Energy; log(E) [GeV]; N_{tow}", 
			nbins, energy_bins.data()); 
	//1 MeV -15 GeV, to keep comparable between low and high energy towers 
	
	//Energy v Time 
	TH1F* EBarT = new TH1F(
			"E_bar_T", "E versus < #Delta T >; log(E) [GeV]; < #Delta T> [ns]",
			nbins, energy_bins.data());
	TH1F* EBarT_high = new TH1F(
			"E_bar_T_high", "E versus < #Delta T >; log(E) [GeV]; < #Delta T> [ns]",
			nbins, energy_bins.data());
	TH1F* EBarT_low = new TH1F(
			"E_bar_T_low", "E versus < #Delta T >; log(E) [GeV]; < #Delta T> [ns]",
			nbins, energy_bins.data());

	//add to the vectors 
	all1DHists->push_back(DeltaT);
	high1DHists->push_back(DeltaT_high);
	low1DHists->push_back(DeltaT_low);
	
	all1DHists->push_back(Energy);
	high1DHists->push_back(Energy_high);
	low1DHists->push_back(Energy_low);

	all1DHists->push_back(EBarT);
	high1DHists->push_back(EBarT_high);
	low1DHists->push_back(EBarT_low);


	TH2F* ETall 	= new TH2F(
		       "E_to_T_all", "Energy #Delta T correlation; log(E) [GeV]; #Delta T [ns]; N_{tow}",
			nbins, energy_bins.data(), 100, -20., 20.);      
	TH2F* EThigh 	= new TH2F(
		       "E_to_T_high", "Energy #Delta T correlation; log(E) [GeV]; #Delta T [ns]; N_{tow}",
			100, -20, 20, nbins, energy_bins.data());      
	TH2F* ETlow 	= new TH2F(
		       "E_to_T_low", "Energy #Delta T correlation; log(E) [GeV]; #Delta T [ns]; N_{tow}",
			100, -20, 20, nbins, energy_bins.data());      
	
	TH2F* EP_all 	= new TH2F(
		       "EtaPhiT_all", ";#eta; #varphi; < #Delta T >",
			96, -1.1, 1.1, 256, 0, 2*M_PI);      
	TH2F* EP_high 	= new TH2F(
		       "EtaPhiT_high", ";#eta; #varphi; < #Delta T >",
			);      
	TH2F* EP_low 	= new TH2F(
		       "EtaPhiT_all", ";#eta; #varphi; < #Delta T >",
			96, -1.1, 1.1, 256, 0, 2*M_PI);      
	all2DHists->push_back(ETall);
	high2DHists->push_back(EThigh);
	low2DHists->push_back(ETlow);
	all2DHists->push_back(EP_all);
	high2DHists->push_back(EP_high);
	low2DHists->push_back(EP_low);
}

int EMCALChannelTiming::getIndex(float phi, float eta)
{
	int index = 0;
	int neta = 4 * 24;
	int nphi = 4 * 64;
	float binphi = 2*M_PI/(float)nphi;
	float bineta = -2.2 /(float)neta;
	int iphi = phi * binphi;
	int ieta = eta * bineta;
	index = neta * iphi + ieta;
       return index;
}       

int EMCALChannelTiming::process_event(PHCompositeNode* topNode)
{
	//Need to get all towers from the emcal
	//Get the towers, subdivide and send them to a helper function 
	std::vector<tower*>* lowE  = new std::vector<tower*>{};
	std::vector<tower*>* highE = new std::vector<tower*>{};
	float avgtime = SubdivideDetector(lowE, highE, allE, topNode);
	AnaHelper(*allE, AllTowers1D, AllTowers2D, avgtime);
	AnaHelper(*highE, highTowers1D, highTowers2D, avgtime);
	AnaHelper(*lowE, lowTowers1D, lowTowers2D, avgtime);
	for(int n=0; n<(int) allE->size(); n++)
	{
		int index = getIndex(allE->at(n)->phi, allE->at(n)->eta);
		AnaHelper(allE->at(n), AllTowers1D_t->at(index), AllTowers2D_t->at(index), avgtime);
	}
	towerTree->Fill();	
	return Fun4AllReturnCodes::EVENT_OK;
}
float EMCALChannelTiming::SubdivideDetector(
		std::vector<tower*>* lowerEtowers, 
		std::vector<tower*>*  higherEtowers, 
		std::vector<tower*>* allTowers,
		PHCompositeNode* topNode
		)
{
	float avg_time = 0;
	auto emcaltowers = findNode::getClass<TowerInfoContainerv3>( topNode, emcal_tower );
	auto emcalgeom = findNode::getClass<RawTowerGeomContainer_Cylinderv1>( topNode, emcal_geom );

	emcalgeom->set_calorimeter_id(RawTowerDefs::CEMC);
	for(int n=0; n<(int)emcaltowers->size(); n++)
	{
		auto key   = emcaltowers->encode_key(n);
		int phibin = emcaltowers->getTowerPhiBin(key);
		int etabin = emcaltowers->getTowerEtaBin(key);
		float phi  = emcalgeom->get_phicenter(phibin);
		float eta  = emcalgeom->get_etacenter(etabin);
		auto tow   = emcaltowers->get_tower_at_channel(key);
		float e	   = tow->get_energy();
		float t    = tow->get_time();
		int N	   = tow->get_nsample();
		int16_t Ei = 0;
		avg_time  += t; 
		for(int i= 0; i< N; i++)
		{
			auto temp = tow->get_waveform_value(i);
			if(temp > Ei) Ei=temp;
		}
		if(Ei > 30 and E < 5000) continue;

		tower* twA = new tower(phi, eta, e, Ei, t);
		
		if( Ei <= 30 ) lowerEtowers->push_back(twA);
		else if( E >= 5000 ) higherEtowers->push_back(twA);
		
		allTowers->push_back(twA);
	}

	avg_time = avg_time/((float) emcaltowers->size());
	
	return avg_time; 
}
void EMCALChannelTiming::AnaHelper(
		std::vector<tower*> subsettowers,
		std::vector<TH1F*>* a1Doutput,
		std::vector<TH2F*>* a2Doutput, 
		float TAvg
		)
{
	std::map<int, std::pair<int, float>> tavg {};
	for(
		int i = 0; 
		i< (int)a1Doutput->at(a1DOUTPUTHISTS::E)->getNbinsX(); 
		i++
	)
	{
		tavg[i]=std::make_pair(0, 0.);
	}

	for(auto tow: subsettowers)
	{
		a1Doutput->at(a1DOUTPUTHISTS::DELTAT)
			->Fill(tow->t - TAvg);
		a1Doutput->at(a1DOUTPUTHISTS::E)->
			Fill(E);
		int binN = a1Doutput->at(a1DOUTPUTHISTS::E)->FindBin(E)
		tavg[binN].first++;
		tavg[binN].second+=tow->t - TAvg;
		a2Doutput->at(a2DOUTPUTHISTS::EtoT)->Fill(E, tow->t - TAvg);
		a2Doutput->at(a2DOUTPUTHISTS::EtaPhi)->Fill(tow->eta, tow->phi, tow->t - TAvg);
	}
	for(auto m:tavg)
	{
		m.second.second= m.second.second/(float)m.second.first;
		float Ec= a1Doutput->at(a1DOUTPUTHISTS::E)->GetBinCenter(m.first);
		a1Doutput->at(a1DOUTPUTHISTS::EBART)->Fill(Ec, m.second.second);
	}
	return;
}
void EMCALChannelTiming::AnaHelper(
		tower* tow,
		std::vector<TH1F*>* a1Doutput,
		std::vector<TH2F*>* a2Doutput, 
		float TAvg
		)
{
	std::map<int, std::pair<int, float>> tavg {};

	a1Doutput->at(a1DOUTPUTHISTS::DELTAT)
		->Fill(tow->t - TAvg);
	a1Doutput->at(a1DOUTPUTHISTS::E)->
		Fill(E);
	a1Doutput->at(a1DOUTPUTHISTS::EBART)->Fill(E, tow->t - TAvg);
	a2Doutput->at(a2DOUTPUTHISTS::EtoT)->Fill(E, tow->t - TAvg);
}

void EMCALChannelTiming::Print(const std::string &what) const
{
	TFile* f = new TFile(output_file_name.c_str(), "RECREATE");
	f->cd();
	T->Write();
	TDirectory* td = new TDirectory("EMCAL_TOWS", "EMCAL_TOWS");
	td->cd();
	for (int i = 0; i<(int)AllTowers1D_t->size(); i++);
	{
		TDirectory* td_T=new TDirectory(std::format("Tower_{}", i).c_str(), std::format("Tower_{}", i).c_str();
		td_T->cd();
		for(int j=0; j<(int)AllTowers1D_t->at(i)->size(); j++{
			AllTowers1D_t->at(i)->at(j)->Write();
			highTowers1D_t->at(i)->at(j)->Write();
			lowTowers1D_t->at(i)->at(j)->Write();
		}
		AllTowers2D_t->at(i)->at(0)->Write();
		highTowers2D_t->at(i)->at(0)->Write();
		lowTowers2D_t->at(i)->at(0)->Write();

		td->cd();
	}
	f->cd();
	for(int i=0; i<(int)AllTowers1D->size(); i++){
		AllTowers1D->at(0)->Write();
		highTowers1D->at(0)->Write();
		lowTowers1D->at(0)->Write();
	}
	for(int i=0; i<(int)AllTowers2D->size(); i++)
	{
		AllTowers2D->at(i)->Write();
		highTowers2D->at(i)->Write();
		lowTowers2D->at(i)->Write();
	}
	f->Write();
	f->Close();
}

