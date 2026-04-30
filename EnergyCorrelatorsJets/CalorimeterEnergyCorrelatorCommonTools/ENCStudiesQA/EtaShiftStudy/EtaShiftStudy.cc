//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in EtaShiftStudy.h.
//
// EtaShiftStudy(const std::string &name = "EtaShiftQA")
// everything is keyed to EtaShiftStudy, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// EtaShiftStudy::~EtaShiftQA()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int EtaShiftStudy::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int EtaShiftStudy::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int EtaShiftStudy::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT; 
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT; 
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int EtaShiftStudy::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int EtaShiftStudy::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int EtaShiftStudy::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int EtaShiftStudy::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void EtaShiftStudy::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
// [[maybe_unused]] suppresses compiler warnings if topNode is not used in this method
//
//____________________________________________________________________________..

#include "EtaShiftStudy.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

EtaShiftStudy::EtaShiftStudy(const std::string ofn, [[maybe_unused]] const std::string &name)
{
	output_file_name = ofn;

	EMCALQA = new PerCaloQAPlots("EMCAL");
	IHCALQA	= new PerCaloQAPlots("IHCAL");
	OHCALQA = new PerCaloQAPlots("OHCAL");
	MetaEQA = new PerCaloQAPlots("Meta Calo, EMCAL radius");
	MetaIQA = new PerCaloQAPlots("Meta Calo, IHCAL radius");
	MetaOQA = new PerCaloQAPlots("Meta Calo, OHCAL radius");
	
	calculatedJetpt 	= new TH1F(
				"h_calc_Jet_pt", 
				"Jet p_{T} calculated from unshifted consituents; p_{T}^{jet}[GeV]; N_{jet}",
				50, -1, 99
			);
	calculatedShiftedJetpt 	= new TH1F(
				"h_shift_calc_Jet_pt", 
				"Jet p_{T} calculated from shifted consituents; p_{T}^{jet} [GeV]; N_{jet}",
				50, -1, 99
			);
	hzVTX			= new TH1F(
				"h_zvtx", 
				"Vertex position; z_{vtx} [cm]; N_{events}",
				300, -150.5, 149.5
			);
				
	for(int i = 0; i<6; i++)
	{
		int z_range = (i+1)*10;
		EMCAL_Z_QA->at(i)	= new PerCaloQAPlots(std::format("EMCAL, |z|_{{vtx}} < {}", z_range));	
		IHCAL_Z_QA->at(i) 	= new PerCaloQAPlots(std::format("IHCAL, |z|_{{vtx}} < {}", z_range));
		OHCAL_Z_QA->at(i)  	= new PerCaloQAPlots(std::format("OHCAL, |z|_{{vtx}} < {}", z_range));
		MetaE_Z_QA->at(i)  	= new PerCaloQAPlots(std::format("Meta Calo, EMCAL radius, |z|_{{vtx}} < {}", z_range));
		MetaI_Z_QA->at(i)  	= new PerCaloQAPlots(std::format("Meta Calo, IHCAL radius, |z|_{{vtx}} < {}", z_range));
		MetaO_Z_QA->at(i)  	= new PerCaloQAPlots(std::format("Meta Calo, OHCAL radius, |z|_{{vtx}} < {}", z_range));
	}
	//initialize variables from the tree
	emMetaTowerBuilder = new BuildMetaTowers( BuildMetaTowers::CALO::EMCAL, "Fun4AllTowers" );
	ihMetaTowerBuilder = new BuildMetaTowers( BuildMetaTowers::CALO::IHCAL, "Fum4AllTowers" );
	ohMetaTowerBuilder = new BuildMetaTowers( BuildMetaTowers::CALO::OHCAL, "Fun4AllTowers" );

	
}

float EtaShiftStudy::CalculateJetPt(std::vector<std::array<float,2>> const_pxy )
{
	//using E scheme
	float jet_pt=0.;
	std::array<float,2> jet_pxy {0., 0.};
	for(auto c:const_pxy)
	{
		jet_pxy[0]+=c[0];
		jet_pxy[1]+=c[1];
	}
	jet_pt=std::sqrt(std::pow(jet_pxy[0], 2) + std::pow(jet_pxy[1], 2));
	return jet_pt;
}

void EtaShiftStudy::AnalyzeEvent(PHCompositeNode* topNode)
{
	std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 4>* rawTowersEM	
		= new std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 4> {}; 
	std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 4>* rawTowersIH	
		= new std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 4> {}; 
	std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 4>* rawTowersOH	
		= new std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 4> {}; 
	for(int i = 0; i<4; i++){
		rawTowersEM->at(i)=new std::array<BuildMetaTowers::TowerArrayEntry*, 1536> {};
		rawTowersIH->at(i)=new std::array<BuildMetaTowers::TowerArrayEntry*, 1536> {};
		rawTowersOH->at(i)=new std::array<BuildMetaTowers::TowerArrayEntry*, 1536> {};
	}
	std::cout<<__LINE__<<std::endl;
	float zvtx = 0.;
        try{
                GlobalVertexMap* vertexmap=findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
                if(vertexmap){
                        if(vertexmap->empty())
                                std::cout<<"Empty Vertex Map. \n Setting vertex to origin" <<std::endl;
                        else{

                                GlobalVertex* gl_vtx=nullptr;
                                for(auto vertex_iter:*vertexmap){
                                        if(vertex_iter.first == GlobalVertex::VTXTYPE::MBD || vertex_iter.first == GlobalVertex::VTXTYPE::SVTX_MBD )
                                        {
                                                gl_vtx=vertex_iter.second;
                                        }
                                }
                                if(gl_vtx){
                                        zvtx=gl_vtx->get_z();
                                }
                        }
                }
        }
        catch(std::exception& e){std::cout<<"Could not find the vertex. \n Setting to origin" <<std::endl;}
	std::cout<<__LINE__<<std::endl;	
	hzVTX->Fill(zvtx);
	std::cout<<__LINE__<<std::endl;	
	//Load in the recoTowers, handle seperating the calos offline
	rawTowersEM	= emMetaTowerBuilder->LoadFun4AllTowers(topNode);
	rawTowersIH	= ihMetaTowerBuilder->LoadFun4AllTowers(topNode);
	rawTowersOH	= ohMetaTowerBuilder->LoadFun4AllTowers(topNode);
	std::cout<<__LINE__<<std::endl;	
	//Build the MetaTowers	
	emMetaTowerBuilder->RunMetaTowerBuilder(zvtx);
	ihMetaTowerBuilder->RunMetaTowerBuilder(zvtx);
	ohMetaTowerBuilder->RunMetaTowerBuilder(zvtx);
	std::cout<<__LINE__<<std::endl;	

	//Get the tower info and load it in 
	std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 4>* emMetaTowers	= new std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 4> {};    	
	std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 4>* ihMetaTowers = new std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 4> {};    	
	std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 4>* ohMetaTowers = new std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 4> {};    	
	std::cout<<__LINE__<<std::endl;	
	
	grabTowerArray(emMetaTowerBuilder, emMetaTowers);
	grabTowerArray(ihMetaTowerBuilder, ihMetaTowers);
	grabTowerArray(ohMetaTowerBuilder, ohMetaTowers);
	std::cout<<__LINE__<<std::endl;	
	
	//need to hand it the full MetaTower set then break it down
	compareTowerValue( emMetaTowers->at(3), rawTowersEM->at(3), zvtx, 0, emMetaTowerBuilder );
	compareTowerValue( ihMetaTowers->at(3), rawTowersIH->at(3), zvtx, 1, ihMetaTowerBuilder );
	compareTowerValue( ohMetaTowers->at(3), rawTowersOH->at(3), zvtx, 2, ohMetaTowerBuilder );
	std::cout<<__LINE__<<std::endl;	
	
	//do the jet analysis
	grabJetConstituents(topNode, zvtx);	
	std::cout<<__LINE__<<std::endl;	
	
		
	return;
}	

void EtaShiftStudy::grabTowerArray(BuildMetaTowers* MTB, std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 4>* caloSet)
{
	caloSet->at(0)	= MTB->getEMReTowers();
	caloSet->at(1)	= MTB->getIHCaTowers();
	caloSet->at(2)	= MTB->getOHCaTowers();
	caloSet->at(3)	= MTB->getMetaTowers();
	return;
}
void EtaShiftStudy::compareTowerValue(std::array<BuildMetaTowers::TowerArrayEntry*, 1536>* CaloShifted, std::array<BuildMetaTowers::TowerArrayEntry*, 1536>* CaloBase, float zvtx, int meta, BuildMetaTowers* MTB)
{
	PerCaloQAPlots* CPQA {nullptr};
	std::array<PerCaloQAPlots*, 6>* ZRest {nullptr};
	if( meta == 0 ) 
	{
		CPQA 	= MetaEQA;
		ZRest	= MetaE_Z_QA;
	}
	else if( meta == 1 ) 
	{
		CPQA	= MetaIQA;
		ZRest	= MetaI_Z_QA;
	}
	else if( meta == 2 ) 
	{
		CPQA	= MetaOQA;
		ZRest	= MetaO_Z_QA;
	}

	int non_neg = 0;
	float avg_eta = 0.;
	float EtSum = 0.;
	float DeltaEtSum = 0.;
	for(int i=0; i<(int)CaloBase->size(); i++)
	{
		auto CaloShiftTower	= CaloShifted->at(i);
		auto CaloBaseTower	= CaloBase->at(i);
		if( CaloShiftTower->Energy <= 0 || CaloBaseTower->Energy <= 0 ) continue;
		non_neg++;
		float E		= CaloShiftTower->Energy;
		float eta	= CaloShiftTower->eta;
		float Et	= E / std::sinh(eta);
		float DeltaEt	= Et - (CaloBaseTower->Energy)/std::sinh(CaloBaseTower->eta);
		float Deltaeta	= CaloShiftTower->eta - CaloBaseTower->eta;
		avg_eta 	+= eta;
		EtSum		+= Et;
		DeltaEtSum	+= DeltaEt;
		CPQA->shifteta->Fill(Deltaeta);
		for(int j=0; j<(int)ZRest->size(); j++)
		{
			if(std::abs(zvtx) <  10*j+10){
				ZRest->at(j)->shifteta->Fill(Deltaeta);
			}
		}

	}
	avg_eta	= avg_eta/(float)non_neg;
	CPQA->avgeta->Fill(avg_eta);
	CPQA->deltaEt->Fill(DeltaEtSum);
	CPQA->Et->Fill(EtSum);
	CPQA->zVTXdeltaeta->Fill(zvtx, avg_eta);
	for(int i=0; i<(int)ZRest->size(); i++)
	{
		if(std::abs(zvtx) < 10*(i+1)){
			ZRest->at(i)->avgeta->Fill(avg_eta);
			ZRest->at(i)->deltaEt->Fill(DeltaEtSum);
			ZRest->at(i)->Et->Fill(EtSum);
			ZRest->at(i)->zVTXdeltaeta->Fill(zvtx, avg_eta);
		}
	}
	std::array<double, 25> shiftedEdge = MTB->getShiftedEdges();
	for(int i = 0; i<(int) shiftedEdge.size() - 1; i++)
	{
		double low 	= shiftedEdge[i];
		double high 	= shiftedEdge[i+1];
		double diff	= std::abs(high - low);
		double ddiff	= diff - 0.0917;
	        CPQA->Deltaetabin->Fill(ddiff);
                for(int j=0; j<(int)ZRest->size(); j++)
                {
                        if(std::abs(zvtx) <  10*j+10){
                                ZRest->at(j)->Deltaetabin->Fill(ddiff);
                        }
                }
	}
	delete CPQA;
	delete ZRest;
	return;	
}
void EtaShiftStudy::grabJetConstituents( PHCompositeNode* topNode, float zVTX)
{
	//grab the jet consituents and save the momentum 
	auto jetConts	= findNode::getClass<JetContainerv1>(topNode, "AntiKt_r04_calib");
	if(!jetConts) return;
	if(jetConts->size() <1) return;
	std::string ohcal_energy_towers	= "TOWERINFO_CALIB_HCALOUT";
	std::string ihcal_energy_towers	= "TOWERINFO_CALIB_HCALIN"; 
	std::string emcal_energy_towers	= "TOWERINFO_CALIB_CEMC_RETOWER";
	auto emcal_tower_energy	= findNode::getClass<TowerInfoContainer>(topNode,  emcal_energy_towers );
	auto ihcal_tower_energy	= findNode::getClass<TowerInfoContainer>(topNode, ihcal_energy_towers );
	auto ohcal_tower_energy	= findNode::getClass<TowerInfoContainer>(topNode,  ohcal_energy_towers );
	auto ohcal_geom	= findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALOUT");
	auto emcal_geom	= findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_CEMC"   );
	auto ihcal_geom	= findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALIN");
	ohcal_geom->set_calorimeter_id(RawTowerDefs::HCALOUT);
	ihcal_geom->set_calorimeter_id(RawTowerDefs::HCALIN);
	emcal_geom->set_calorimeter_id(RawTowerDefs::CEMC);
	for(auto jet: *jetConts){
		if(!jet) continue;
		std::vector< std::array< float, 2> > jetC {}; 
		std::array<float, 2> pxy {0., 0.};
		std::vector< std::array< float, 2> > jetCsft {}; 
		std::array<float, 2> pxysft {0., 0.};
		auto cmp_vec = jet->get_comp_vec();
		for(auto iter:cmp_vec)
		{
			unsigned int tower_id=iter.second;
			float px= 0.;
		       	float py= 0.;
			float e	= 0.;
			float phicenter	= 0.;
			float etacenter = 0.;
			float r	= 1.;
			Jet::SRC source=iter.first;
			if( source == Jet::SRC::HCALOUT_TOWER)
			{
				try{
					e		= ohcal_tower_energy->get_tower_at_channel(tower_id)->get_energy();
					int phibin	= ohcal_tower_energy->getTowerPhiBin(tower_id);
					int etabin	= ohcal_tower_energy->getTowerEtaBin(tower_id);
					phicenter	= ohcal_geom->get_phicenter(phibin);
					etacenter	= ohcal_geom->get_etacenter(etabin);	
					r 		= ohMetaTowerBuilder->getRadius();
				
				}
				catch(std::exception& x){
					std::cout<<"Bad tower id found for source " <<Jet::SRC::HCALOUT_TOWER 
						<<" despite actual source being " <<source 
					<<" and tower id is " <<tower_id<<std::endl;
				}
			}
			else if( source == Jet::SRC::HCALIN_TOWER)
			{
				try{
					e		= ihcal_tower_energy->get_tower_at_channel(tower_id)->get_energy();
					int phibin	= ihcal_tower_energy->getTowerPhiBin(tower_id);
					int etabin	= ihcal_tower_energy->getTowerEtaBin(tower_id);
					phicenter	= ihcal_geom->get_phicenter(phibin);
					etacenter	= ihcal_geom->get_etacenter(etabin);	
					r 		= ihMetaTowerBuilder->getRadius();

				}
				catch(std::exception& x){
					std::cout<<"Bad tower id found for source " <<Jet::SRC::HCALOUT_TOWER 
						<<" despite actual source being " <<source 
					<<" and tower id is " <<tower_id<<std::endl;
				}
			}
			else if( source == Jet::SRC::CEMC_TOWER)
			{
				try{
					e		= emcal_tower_energy->get_tower_at_channel(tower_id)->get_energy();
					int phibin	= emcal_tower_energy->getTowerPhiBin(tower_id);
					int etabin	= emcal_tower_energy->getTowerEtaBin(tower_id);
					phicenter	= emcal_geom->get_phicenter(phibin);
					etacenter	= emcal_geom->get_etacenter(etabin);	
					r 		= emMetaTowerBuilder->getRadius();
				}
				catch(std::exception& x){
					std::cout<<"Bad tower id found for source " <<Jet::SRC::HCALOUT_TOWER 
						<<" despite actual source being " <<source 
					<<" and tower id is " <<tower_id<<std::endl;
				}
			}
			else 
				std::cout<<"Couldn't deal with the source: " <<source <<"\n Skipping this particle" <<std::endl;  
			float sfteta	= asinh((r*sinh(etacenter)+zVTX)/r);
			float pt	= e / cosh(etacenter); 
			float ptsft	= e / cosh(sfteta);
			px		= pt * cos(phicenter);
			py		= pt * sin(phicenter);
			pxy[0]		= px;
			pxy[1]		= py;
			jetC.push_back(pxy);
		
			px		= ptsft * cos(phicenter);
			py		= ptsft * sin(phicenter);
			pxysft[0]	= px;
			pxysft[1]	= py;
			jetCsft.push_back(pxysft);	
		}
		calculatedJetpt->Fill(CalculateJetPt(jetC));
	}
	return;

}

int EtaShiftStudy::process_event(PHCompositeNode *topNode)
{
	//if(vebose) std::cout << "EtaShiftStudy::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
	AnalyzeEvent(topNode); 
  	return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..

int EtaShiftStudy::Reset( [[maybe_unused]] PHCompositeNode *topNode)
{
  std::cout << "EtaShiftStudy::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EtaShiftStudy::End([[maybe_unused]] PHCompositeNode *topNode)
{
 	std::cout << "EtaShiftStudy::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
 	TFile* output_file = new TFile(output_file_name.c_str(), "RECREATE");
	output_file->cd();
	calculatedJetpt->Write();
	calculatedShiftedJetpt->Write();
	hzVTX->Write();
	std::array<TDirectory*, 7>* base_dirs = new std::array<TDirectory*, 7> {};
	std::array<std::array<TDirectory*, 6>*, 7>* deep_dirs = new std::array<std::array<TDirectory*, 6>*, 7>{}; 
	for(int i = 0; i<(int)base_dirs->size(); i++)
	{
		base_dirs->at(i) = output_file->mkdir(std::format("Z_leq_{}_cm", 10*(i+1)).c_str()); 
		deep_dirs->at(i)->at(0) = base_dirs->at(i)->mkdir("EMCAL");
		if(i==0) EMCALQA->dumpThePlots(deep_dirs->at(i)->at(0));	
		else EMCAL_Z_QA->at(i-1)->dumpThePlots(deep_dirs->at(i)->at(0));	
		
		deep_dirs->at(i)->at(1)	= base_dirs->at(i)->mkdir("IHCAL");
		if(i==0) IHCALQA->dumpThePlots(deep_dirs->at(i)->at(1));	
		else IHCAL_Z_QA->at(i-1)->dumpThePlots(deep_dirs->at(i)->at(1));	
		
		deep_dirs->at(i)->at(2)	= base_dirs->at(i)->mkdir("OHCAL");
		if(i==0) OHCALQA->dumpThePlots(deep_dirs->at(i)->at(2));	
		else OHCAL_Z_QA->at(i-1)->dumpThePlots(deep_dirs->at(i)->at(2));	
		
		deep_dirs->at(i)->at(3)	= base_dirs->at(i)->mkdir("Meta EMCal");
		if(i==0) MetaEQA->dumpThePlots(deep_dirs->at(i)->at(3));	
		else MetaE_Z_QA->at(i-1)->dumpThePlots(deep_dirs->at(i)->at(3));	
		
		deep_dirs->at(i)->at(4)	= base_dirs->at(i)->mkdir("Meta IHCal");
		if(i==0) MetaIQA->dumpThePlots(deep_dirs->at(i)->at(4));	
		else MetaI_Z_QA->at(i-1)->dumpThePlots(deep_dirs->at(i)->at(4));	
		
		deep_dirs->at(i)->at(5)	= base_dirs->at(i)->mkdir("Meta OHCal");
		if(i==0) MetaOQA->dumpThePlots(deep_dirs->at(i)->at(5));	
		else MetaO_Z_QA->at(i-1)->dumpThePlots(deep_dirs->at(i)->at(5));	
		output_file->cd();

	}
	output_file->cd();
	for(int i = 0; i < (int)base_dirs->size(); i++)
	{
		base_dirs->at(i)->cd();
		for(int j = 0; j < (int)deep_dirs->at(i)->size(); j++)
		{
			deep_dirs->at(i)->at(j)->Write();
		}
		output_file->cd();
		base_dirs->at(i)->Write();
	}
	output_file->Write();
	output_file->Close();
  	return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void EtaShiftStudy::Print(const std::string &what) const
{
  std::cout << "EtaShiftStudy::Print(const std::string &what) const Printing info for " << what << std::endl;
}
