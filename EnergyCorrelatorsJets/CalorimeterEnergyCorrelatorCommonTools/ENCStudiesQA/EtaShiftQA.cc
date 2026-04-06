#include "EtaShiftQA.h"

EtaShiftQA::EtaShiftQA(TTree* T)
{
	EMCALQA = new PerCaloQAPlots("EMCAL");
	IHCALQA	= new PerCaloQAPlots("IHCAL");
	OHCALQA = new PerCaloQAPlots("OHCAL");
	MetaEQA = new PerCaloQAPlots("Meta Calo, EMCAL radius");
	MetaIQA = new PerCaloQAPlots("Meta Calo, IHCAL radius");
	MetaOQA = new PerCaloQAPlots("Meta Calo, OHCAL radius");
	
	caluclatedJetpt 	= new TH1F(
				"h_calc_Jet_pt", 
				"Jet p_{T} calculated from unshifted consituents; p_{T}^{jet}[GeV]; N_{jet}",
				50, -1, 99
			);
	caluclatedShiftedJetpt 	= new TH1F(
				"h_shift_calc_Jet_pt", 
				"Jet p_{T} calculated from shifted consituents; p_{T}^{jet} [GeV]; N_{jet}",
				50, -1, 99
			);
	hzVTX			= new TH1F(
				"h_zvtx", 
				"Vertex position; z_{vtx} [cm]; N_{events}"
				300, -150.5, 149.5
			);
				
	for(int i = 0; i<6; i++)
	{
		int z_range = (i+1)*10;
		EMCAL_Z_QA->at(i)	= new PerCaloQAPlots(std::format("EMCAL, |z|_{vtx} < {}", z_range));
		IHCAL_Z_QA->at(i) 	= new PerCaloQAPlots(std::format("IHCAL, |z|_{vtx} < {}", z_range));
		OHCAL_Z_QA->at(i)  	= new PerCaloQAPlots(std::format("OHCAL, |z|_{vtx} < {}", z_range));
		MetaE_Z_QA->at(i)  	= new PerCaloQAPlots(std::format("Meta Calo, EMCAL radius, |z|_{vtx} < {}", z_range));
		MetaI_Z_QA->at(i)  	= new PerCaloQAPlots(std::format("Meta Calo, IHCAL radius, |z|_{vtx} < {}", z_range));
		MetaO_Z_QA->at(i)  	= new PerCaloQAPlots(std::format("Meta Calo, OHCAL radius, |z|_{vtx} < {}", z_range));
	}
	//initialize variables from the tree
	T->SetBranchAddress("EventInfo",&eventInfo);
	T->SetBranchAddress("JetInfo_r04",&recoJets);
	T->SetBranchAddress("TowerInfo",&recoTowers);
	if(sim)
	{
		T->SetBranchAddress("TruthParticles",&truthParticles);
		T->SetBranchAddress("TruthJetInfo_r04",&truthJets);
	}
	emMetaTowerBuilder = new BuildMetaTowers();
	ihMetaTowerBuilder = new BuildMetaTowers();
	ohMetaTowerBuilder = new BuildMetaTowers();

	emMetaTowerBuilder->setCaloRadius(BuildMetaTowers::CALO::EMCAL)
	ihMetaTowerBuilder->setCaloRadius(BuildMetaTowers::CALO::IHCAL)
	ohMetaTowerBuilder->setCaloRadius(BuildMetaTowers::CALO::OHCAL)
}

float EtaShiftQA::CalculateJetPt(std::vector<stdd::array<float,2>> const_pxy )
{
	//using E scheme
	float jet_pt=0.;
	std::array<float,2> jet_pxy {0., 0.};
	for(auto c:const_pt) jet_pxy+=c;
	jet_pt=std::sqrt(std::pow(jet_pxy[0], 2) + std::pow(jet_pxy[1], 2));
	return jet_pt;
}

void EtaShiftQA::AnalyzeEvent(int n_evt)
{
	T->GetEntry(n_evt);
	float zvtx = eventInfo->get_z_vtx();
	//Load in the recoTowers, handle seperating the calos offline
	emMetaTowerBuilder->LoadVandyTowers(*recoTowers, zvtx);
	ihMetaTowerBuilder->LoadVandyTowers(*recoTowers, zvtx);
	ohMetaTowerBuilder->LoadVandyTowers(*recoTowers, zvtx);
	//Build the MetaTowers	
	emMetaTowerBuilder->RunMetaTowerBuilder(zvtx);
	ihMetaTowerBuilder->RunMetaTowerBuilder(zvtx);
	ohMetaTowerBuilder->RunMetaTowerBuilder(zvtx);

	//Get the tower info and load it in 
	
}	



