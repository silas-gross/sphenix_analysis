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

	
}
boid JetShapeEventObservables::getTruth(PHCompositeNode* topNode) 
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

