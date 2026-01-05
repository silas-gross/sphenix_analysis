#include "JetShapeEventObservables.h"

JetShapeEventObservables::JetShapeEventObservables()
{
}

int JetShapeEventObservables::process_event(PHCompositeNode* topNode)
{
}

void JetShapeEventObservables::filterTowerInput(PHCompositeNode* topNode) 
{
	//This is just getting things ready for the trimming
	//need to look for Truth, Clusters, Towers
	if(this->truth) getTruth(topNode);
	if(this->cluster) getClusters(topNode);
	if(this->tower) getTowers(topNode);
	return;
}
void JetShapeEventObservables::getTowers(PHCompositeNode* topNode)
{
	auto emcal_geom=findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_CEMC"   );
	auto ihcal_geom=findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALIN"   );
	auto ohcal_geom=findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALOUT"   );
	auto emcal_tower_energy=findNode::getClass<TowerInfoContainer>(topNode, emcal_energy_towers      );
	auto ihcal_tower_energy=findNode::getClass<TowerInfoContainer>(topNode, ihcal_energy_towers      );
	auto ohcal_tower_energy=findNode::getClass<TowerInfoContainer>(topNode, ohcal_energy_towers      );
	
	if(this->emcal_lookup_table.size() == 0  || n_evts==1)
	{
		MakeEMCALRetowerMap(emcal_geom, emcal_tower_energy, ohcal_geom, ohcal_tower_energy);
		std::cout<<"Lookup table has size: " <<this->emcal_lookup_table.size() <<std::endl;
	}
	std::map<int, std::array<float, 4>> emcal {}, ihcal {}, ohcal {}, emcal_retower {}, all_cal {};
	emcal_geom->set_calorimeter_id(RawTowerDefs::CEMC);
	ihcal_geom->set_calorimeter_id(RawTowerDefs::HCALIN);
	ohcal_geom->set_calorimeter_id(RawTowerDefs::HCALOUT);
	for(int n=0; n<(int) emcal_tower_energy->size(); n++)
	{
			
	return;	
}
void JetShapeEventObservables::getTruth(PHCompositeNode* topNode) 
{
	//if there is a truth node, grab it 
	std::map<int, std::array<float, 4>> truth_phg4 {}, truth_hepmc{};  
	try{ 
		findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
	}
	catch(std::exception& e ) {return; }
	PHG4TruthInfoContainer *truthinfo=findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
	if(!truthinfo) return;
	else if (truthinfo) 
	{
		PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
		for(PHG4TruthInfoContainer::ConstIterator iter = range.first, iter !=range.second; ++iter) 
		{
			PHG4Particle *part = iter->second;
			if(!part) continue;
			if(abs(part->get_pid()) >= 12  && abs(part->get_pid) <= 18) continue;
			float E = part->get_e();
			float px = part->get_px();
			float py = part->get_py();
			float pz = part->get_pz();
			float phi=atan2(py, px)+PI;
			float eta=atanh(pz/E);
			float r = 1.;
			int id = part->get_barcode();
			truth_phg4[id]=std::array<float, 4> {eta, phi, r, E};
		}
		filterInputTowers["PHG4_sPHENIX_PRIMARY_TRUTH"]=truth_phg4;
	}
	else return;
	try{
		findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
	}
	catch(std::exception& e) { return;}
	auto *hepmcinfo=findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
	if(!hepmcinfo)return;
	else if(hepmcinfo)
	{
		//as a backup look for the hepmc event
		std::map<int, std::array<float, 4>> truth_hepmc{};
		for(PHHepMCGenEventMap::ConstIter evtIter = hepmcinfo->begin(); evtIter != hepmcinfo->end(); ++evtIter)
		{
			PHHepMCGenEvent* hep_evt=evtIter->second;
			if(hep_evt)
			{
				HepMC::GenEvent* ev=hep_evt->getEvent();
				if(ev)
				{
					for(HepMC::GenEvent::particle_const_iterator iter = ev->particles_begin(); iter !=ev->particles_end(); ++iter)
					{
						if(!(*iter)) continue;
						auto particle = *iter;
						if(!(particle->end_vertex()) && particle->status() == 1) 
						{
							if(abs(particle->pdg_id()) >= 12 && abs(parrticle->pdg_id()) <=18) continue;
							float px = particle->momentum.px();
							float py = particle->momentum.py();
							float pz = particle->momentum.pz();
							float E  = particle->momentum.E();
							float phi=atan2(py, px)+PI;
							float eta=atanh(pz/E);
							float r = 1.;
							int id = particle->get_barcode();
							truth_hepmc[id]=std::array<float, 4> {eta, phi, r, E};
						}
					}
				}
			}
		}
		filterInputTowers["HepMCGenEvent_TRUTH"]=truth_hepmc;
	}
	else return;
	return;
}
void JetShapeEventObservables::getClusters(PHCompostieNode* topNode)
{
       	//find and grab the clusters
	auto emcal_geom=findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_CEMC"   );
	auto ihcal_geom=findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALIN"   );
	auto ohcal_geom=findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALOUT"   );
	auto topoClusters = findNode::getClass<RawClusterContainer>(topNode, "TOPOCLUSTER_ALLCALO");
	auto cm = topoClusters->getClustersMap();
	//int m_cluster_n = 0;
	//std::map<int, float> cluster_e, cluster_phi, cluster_eta; 
	//std::map<int, int> cluster_multipl;
	std::map<int, std::array<float, 4>> cluster_e {}, emcal_e {}, ihcal_e {}, ohcal_e {};
	for(auto entry: cm) 
	{
		RawCluster* cluster = entry.second;
		CLHEP::Hep3Vector origin (vtx[0], vtx[1], vtx[2]);
		float cl_e=cluster->get_energy();
		float cl_eta= RawClusterUtility::GetPseudorapidity(*cluster, origin);
		float cl_phi= RawClusterUtility::GetAzimuthAngle(*cluster, origin);
		int ohcal_cluster=0, emcal_cluster=0;
		float r = cluster->get_r();
		int i = cluster->get_id();
		cluster_e[i]=std::array<float, 4> {cl_eta, cl_phi, r, cl_e};
		for(const auto& [tower_id, tower_e] : cluster->get_towermap())
		{
			RawTowerDefs::CalorimeterId calo_id = RawTowerDefs::decode_caloid(tower_id);
			int iphi = RawTowerDefs::decode_index2(tower_id);
			int ieta = RawTowerDefs::decode_index1(tower_id);
			float eta = 0., phi=0.;
			std::array<float, 4> tow {eta, phi, r, 0.};
			if(calo_id == RawTowerDefs::CEMC){
				emcal_cluster++;
				emcal_geom->set_calorimeter_id(calo_id);
				tow[0]=emcal_geom->get_etacenter(ieta);
				tow[1]=emcal_geom->get_phicenter(iphi);
				tow[2]=emcal_geom->get_radius();
				tow[3]=tower_e;
			       	emcal_e[tower_id]=tow;
		     	}
			else if (calo_id == RawTowerDefs::HCALIN){
				ihcal_geom->set_calorimeter_id(RawTowerDefs::decode_caloid(tower_id));
				tow[0]=ihcal_geom->get_etacenter(ieta);
				tow[1]=ihcal_geom->get_phicenter(iphi);
				tow[2]=ihcal_geom->get_radius();
				tow[3]=tower_e;
			       	ihcal_e[tower_id]=tow;
			}
			else if(calo_id == RawTowerDefs::HCALOUT){
				ohcal_cluster++;
				ohcal_geom->set_calorimeter_id(RawTowerDefs::decode_caloid(tower_id));
				tow[0]=ohcal_geom->get_etacenter(ieta);
				tow[1]=ohcal_geom->get_phicenter(iphi);
				tow[2]=ohcal_geom->get_radius();
				tow[3]=tower_e;
			       	ohcal_e[tower_id]=tow;
			}
			else continue;
		}
		h_n_ohcal_clusters->Fill(ohcal_cluster);
		h_n_emcal_clusters->Fill(emcal_cluster);
	}
	this->filterInputTowers["TopoCluster"]=cluster_e;
	this->filterInputTowers["EMCALClusterTowers"]=emcal_e;
	this->filterInputTowers["OHCALClusterTowers"]=ohcal_e;
	this->filterInputTowers["IHCALClusterTowers"]=ihcal_e;
	return;
}	
void JetShapeEventObservables::MakeEMCALRetowerMap(RawTowerGeomContainer_Cylinderv1* em_geom, TowerInfoContainer* emcal, RawTowerGeomContainer_Cylinderv1* h_geom, TowerInfoContainer* hcal )
{
	em_geom->set_calorimeter_id(RawTowerDefs::CEMC);
	h_geom->set_calorimeter_id(RawTowerDefs::HCALOUT);

	for(int n=0; n<(int) emcal->size(); n++){
		auto key=emcal->encode_key(n);
		int phibin=emcal->getTowerPhiBin(key);
		int etabin=emcal->getTowerEtaBin(key);
		float phicenter=em_geom->get_phicenter(phibin);
		float etacenter=em_geom->get_etacenter(etabin);
		for(int j=0; j<(int) hcal->size(); j++)
		{
			bool goodPhi=false, goodEta=false;
			auto key_h=hcal->encode_key(j);
			int phibin_h=hcal->getTowerPhiBin(key_h);
			int etabin_h=hcal->getTowerEtaBin(key_h);
			float phicenter_h=h_geom->get_phicenter(phibin_h);
			float etacenter_h=h_geom->get_etacenter(etabin_h);	
			std::pair<double, double> phi_bounds=h_geom->get_phibounds(phibin_h);
			std::pair<double, double> eta_bounds=h_geom->get_etabounds(etabin_h);
			if(phicenter >= phi_bounds.first && phicenter < phi_bounds.second) goodPhi=true; 
			if(etacenter >= eta_bounds.first && etacenter < eta_bounds.second) goodEta=true;
			if(goodPhi && goodEta){
				this->emcal_lookup_table[n]=std::make_pair(etacenter_h, phicenter_h);
				break;
			}
			else continue;
		}
	}

	return;
}
