#include "ShowerTowerMatching.h"

ShowerTowerMatching::ShowerTowerMatching(const std::string name)
{
	//this is the initializer
	dataTowers = new std::array<BuildMetaTowers::TowerArrayEntry*, 1536> {};
	truthTowers = new std::array<BuildMetaTowers::TowerArrayEntry*, 1536> {};
	dataClusters = new std::array<BuildMetaTowers::TowerArrayEntry*, 1536> {};
	dataTowers = new std::array<BuildMetaTowers::TowerArrayEntry*, 1536> {};
}
void ShowerTowerMatching::buildTowerBins(int n_bins/*=100*/)
{
	tower_bins->clear();
	float first_bin = 1e-01;
	float last_realbin = 50; //use a 50 GeV upper bin 
	float lastbin = 1e+02; //catchall for above 50 GeV
	tower_bins->push_back(1e-06); //1 keV bottom bin
	tower_bins->push_back(firstbin); //100 MeV is the start of the actual bins
	float width = (std::log10(last_realbin) - std::log10(first_bin) )/((float)n_bins-2.);
	for(int i=1; i<n_bins; i++)
	{
		float bin_edge = first_bin * std::pow(10, width * i);
		if(bin_edge == tower_bins->at(i)) continue;
		else if (bin_edge >= last_realbin) break;
		tower_bins->push_back(bin_edge);
	}
	tower_bins->push_back(last_realbin);
	tower_bins->push_back(lastbin);
	return;
}
int ShowerTowerMatching::Init(PHCompositeNode* topNode)
{
	buildtowerbins();
	h_tow_fake = new TH1F("tow_fake", 
			"Fakes on Meta Towers from truth particle showers; E_{T}^{tow} [GeV]; Fake rate", 
			tower_bins->data());
	h_tow_miss = new TH1F("tow_miss",
			"Miss truth particle showers to Meta Towers; E_{T}^{particle} [GeV]; Miss rate",
			tower_bins->data());
	h_cls_fake = new TH1F("cls_fake", 
			"Fakes on Topo Clusterss from truth particle showers; E_{T}^{cls} [GeV]; Fake rate", 
			tower_bins->data());
	h_cls_miss = new TH1F("cls_miss",
			"Miss truth particle showers to Topo Clusterss; E_{T}^{particle} [GeV]; Miss rate",
			tower_bins->data());
	h_tow_fake_tr 	= new TH1F("tow_fake_tr", 
			"Fakes on Meta Towers from truth tower showers; E_{T}^{tow} [GeV]; Fake rate", 
			tower_bins->data());
	h_tow_miss_tr 	= new TH1F("tow_miss_tr",
			"Miss truth tower showers to Meta Towers; E_{T}^{tower} [GeV]; Miss rate",
			tower_bins->data());
	h_cls_fake_tr 	= new TH1F("cls_fake", 
			"Fakes on Topo Clusters from truth tower showers; E_{T}^{cls} [GeV]; Fake rate", 
			tower_bins->data());
	h_cls_miss_tr 	= new TH1F("cls_miss",
			"Miss truth tower showers to Topo Clusters; E_{T}^{tower} [GeV]; Miss rate",
			tower_bins->data());
	return Fun4AllReturnCodes::EVENT_OK;
}
int ShowerTowerMatching::process_event(PHCompositeNode*  topNode)
{
	if(Verbosity()) std::cout<<"Event number " <<n_evt <<std::endl;
	auto truthinfo= findNode::getClass<PHG4TruthInfoContainer*>(topNode, "G4TruthInfo");
	if(!truthinfo) return Fun4AllReturnCodes::EVENT_OK;
	auto shower_R = truthinfo->GetShowerRange();
	std::map<int, PHG4Shower*> showers {}; 
	for(auto iter = shower_R.first; iter != shower_R.second; ++iter)
	{
		if(!iter) continue;
		PHG4Shower* shower = iter->second;
		if(!shower) continue;
		int parent_id = shower->get_parent_particle_id();
		showers[parent_id]=shower;
	}	
	std::vector<PHG4Particle*> unmatched_truth_particles {};
	std::vector<PHG4Shower*> fake_shower {};
	std::map<int, bool> isShowerMatched {};
	for(auto s:showers) isShowerMatched[s.first] = false;
	std::map<PHG4Particle*, Shower*> particles_to_match; 
	for(
			auto iter = truthinfo->GetSPHENIXPrimaryParticleRange().first; 
			iter != truthinfo->GetSPHENIXPrimaryParticleRange().second; 
			++iter
	   )
	{
		if(!iter) continue;
		PHG4Particle* p = iter->second;
		if(!p) continue;
		bool isKinematicGood = KinCuts(p);
		
		if(!isKinematicGood) continue;
		int track_id = p->get_track_id();
		if(showers.find(track_id)== showers.end())
			unmatched_truth_particles.push_back(p);
		else{
			isShowerMatched[track_id]=true;
			Shower* sr = new Shower();
			getParticleShower(
					p, showers[track_id],
				       	sr, topNode);
			particles_to_match[p]=sr;
		}
	}
	for(auto afake:isShowerMatched)
	{
		if(afake.second == true) continue;
		else fakeShower.push_back(showers[afake.first]);
	}
	float truth_zvtx = 0; 
	auto hepmc_gen_event= findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
	if(hepmc_gen_event)
	{
		for( PHHepMCGenEventMap::ConstIter evtIter=hepmc_gen_event->begin(); evtIter != hepmc_gen_event->end(); ++evtIter)
		{
			PHHepMCGenEvent* hpev=evtIter->second;
			if(hpev){
				HepMC::GenEvent* ev=hpev->getEvent();	
				if(ev)
				{
					auto vtx = ev->signal_process_vertex();
					truth_zvtx = vtx->position().z(); 
				}
			}
		}
	}	
	std::map<BuildMetaTowers::TowerArrayEntry*, Shower*> matched_towers {};
	std::vector<BuildMetaTowers::TowerArrayEntry*> unmatched_towers {};
	buildTruthTowers(
			particles_to_match, 
			unmatched_truth_particles, 
			&matched_towers, 
			&unmatched_towers, 
			truth_zvtx );
	
	matchTheTowers(particles_to_match, matched_towers);
	matchTheClusters(particles_to_match, matched_towers);
}
void ShowerTowerMatching::buildTruthTowers(
		std::map < PHG4Particle*, Shower* > matched, 
		std::vector < PHG4Particle* > unmatched, 
		std::map < BuildMetaTowers::TowerArrayEntry*, Shower*>* matched_towers, 
		std::vector <TowerArrayEntry*>* unmatched_towers,
		float truth_zvtx;
		)
{
	BuildMetaTowers* bm = new BuildMetaTowers();
	bm->RunMetaTowerBuilder(truth_zvtx);
	BuildMetaTower* um = new BuildMetaTowers();
	um->RunMetaTowerBuilder(truth_zvtx);
	um->ConvertPhParticles(unmatched);
	unmatched_towers = &(um->getMetaTowers());
	std::vector<PHG4Particle*> matched_particles {};
	std::vector<int> shower_to_tower_index; 
	for(auto pm: matched)
	{	
		matched_particles.push_back(pm.first);
		shower_to_tower_index.push_back(-999);
	}
	bm->ConvertPhParticles(matched_particles, shower_to_tower_index);
	auto mt = bm->getMetaTowers();
	for(auto m:mt)
	{
		Shower* sh  = new Shower();
		*matched_towers[m] = sh;
	}

	for(auto sti: shower_to_tower_index)
	{
		auto m = mt[sti];
		matched_towers->at(m)->addtoShower(matched[matched_particles.at(sti)]);
	}
	return;
}
void ShowerTowerMatching::buildTopoTowers(
		PHCompositeNode* topNode
		)
{
	dataTowers->clear();
	BuildMetaTower* bm = new BuildMetaTowers( BuildMetaTowers::CALO::EMCAL, "Fun4AllTowers");
	bm->LoadFun4AllTowers(topNode);
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
	bm->RunMetaTowerBuilder(zvtx);	
	dataTowers = bm->getMetaTowers();	
	return;
}
void ShowerTowerMatching::matchTheTowers(
		std::map<PHG4Particle*, Shower*> particles,
		std::map<BuildMetaTowers::TowerArrayEntry*, Shower*> truth_towers
		)
{
	//matching the meta towers to the shower 
	std::array<bool, 1536> is_this_real_or_fake {};
	std::array<std::pair<PHG4Particle*, float>, 1536> TowerParticleWeight {};
	std::vector<std::vector<std::pair<TowerArrayEntry*, float>>> ParticleTowerWeight {};
	std::vector<bool> is_this_a_miss {};
	std::vector<float> particle_ET; 
	for(auto i:is_this_real_or_fake) i = false;
	for(auto p:particles)
	{
		bool is_matched = false;
		std::array<float, 2> eB = p.second->get_etaBounds();
		std::array<float, 2> pB = p.second->get_phiBounds();
		std::vector<Tower> sts = p.second->getStruck();
		float pz {p.first->get_pz()};
		float e	 {p.first->get_e()};
		float eta { std::atanh(pz / e)}; 
		float eT = e / std::cosh(eta);
		particle_ET.push_back(eT);
		for(int i=0; i<(int)dataTowers.size(); i++)
		{
			auto tower = dataTowers.at(i);
			float phi = tower->phi;
			float eta = tower->eta;
			if(phi >= pB[0] && phi <= pB[1])
				is_matched = (eta > eB[0] && eB[1] < eB[1]) ? true : false;
		       if(is_matched){
			       is_this_real_or_fake.at(i) = true;
			       std::pair<PHG4Particle*, float> pt_weight {p, 0.};
			       std::pair<TowerArrayEntry*, float> tw {tower, 0.};
			       for(auto t:sts)
				       if(phi >= t.philow && phi <= t.phihigh)
					       if( eta >= t.etalow && eta <= t.etahigh )
						       tw.second+=t.ET;
		       }

		}
 		is_this_a_miss.push_back(is_matched);		
	}
	for(int j=0; j<(int)particle_ET.size(); j++)
	{
		h_truth_all->Fill(particle_ET);
		if(is_this_a_miss.at(j)) h_truth_to_tower_match->Fill(particle_ET.at(i));
		else h_tow_miss->Fill(particle_ET.at(i));
	}
	for(int j=0; j<(int)dataTowers.size(); j++)
	{
		if(datatowers.at(j)->ET > 0) h_tower_all->Fill(dataTowers.at(j).ET);
	return;	
}
void ShowerTowerMatching::setWeight(std::array<std::pair<PHG4Particle*, float>, 1536>* TowerParticleWeight, PHG4Particle* p)
{
	//The weight here should be the portion of the particle energy going into the shower that hits a specific tower, and the contributions from a specific shower to a specific tower	
}
void ShowerTowerMatching::matchTheClusters()
{
}
bool ShowerTowerMatching::KinCuts(PHG4Particle* p)
{
	bool kingood {false};
	if(!p) return kingood;
	float px {p->get_px()};
	float py {p->get_py()};
	float pz {p->get_pz()};
	float e	 {p->get_e()};
	float eta { std::atanh(pz / e)}; 
	bool isEM {false};
	if( std::abs(eta) <= -1.1)
	{
		int pid = std::abs(p->get_pid());
		if( pid == 11 || pid == 13 || pid == 22) 
			isEM = true;
		else if( pid < 11 || pid > 16 ) 
			isEM = false;
		else return kingood;
		float threhold = isEM ? 0.2 : 0.5;
		if(e > threshold) isEM = true;
	}
	return kingood;
}
void ShowerTowerMatching::getParticleShower(
		PHG4Particle* p1, 
		PHG4Shower* s1, 
		Shower* shower, 
		PHCompositeNode* topNode
		)
{
	//this is where we turn a particle into a set of Hits at the calo level
	std::vector<PHG4Hit*> showerHits;
	auto em_hit = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_CEMC");
	auto ih_hit = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_HCALIN");
	auto oh_hit = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_HCALOUT");
	std::map<int, PHG4HitContainer*> hits {};
	if(em_hit) hits[em_hit->GetID()] = em_hit;
	if(ih_hit) hits[ih_hit->GetID()] = ih_hit;
	if(oh_hit) hits[oh_hit->GetID()] = oh_hit;
	
	for(auto iter = s1->begin_g4hit_id(); 
			iter != s1->end_g4hit_id();
			++iter)
	{
		int cont_id = iter->first;
		auto hit_ids = iter->second;
		if(hits.find(cont_id) == hits.end()) continue;
		auto layer_hits = hits[cont_id];
		for(auto ht_id : hit_ids)
		{
			PHG4Hit* hit = layer_hits.findHit(ht_id);
			if(!hit) continue;
			
			float xc = hit->get_avg_x();
			float yc = hit->get_avg_y();
			float zc = hit->get_avg_z();
			float ed = hit->get_edep();

			float xl = hit->get_x(0);
			float xh = hit->get_x(1);
			
			float yl = hit->get_y(0);
			float yh = hit->get_y(1);
			
			float zl = hit->get_z(0);
			float zh = hit->get_z(1);
			
			float rc = std::sqrt(std::pow(xc, 2) + std::pow(yc,2));
			float rl = std::sqrt(std::pow(xl, 2) + std::pow(yl,2));
			float rh = std::sqrt(std::pow(xh, 2) + std::pow(yh,2));
			
			float ec = std::atanh( zc / rc ); 
			float el = std::atanh( zl / rl ); 
			float eh = std::atanh( zh / rh ); 
			
			float pc = std::atan2( yc, xc);
			float pl = std::atan2( yl, xl);
			float ph = std::atan2( yh, xh);
			
			std::array<float, 3> eta { el, ec, eh};	
			std::array<float, 3> phi { pl, pc, ph};	
			
			tower tw { eta, phi, ed, ht_id, true};
			shower->AddTower(tw);
		}
	}
	return;
}
int ShowerTowerMatching::End(PHCompositeNode* topNode)
{
}
 
