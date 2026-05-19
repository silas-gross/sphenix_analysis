#include "ShowerTowerMatching.h"

ShowerTowerMatching::ShowerTowerMatching(const std::string name)
{
	//this is the initializer
}

int ShowerTowerMatching::Init(PHCompositeNode* topNode)
{

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
	buildTruthTowers(particles_to_match, unmatched_truth_particles);
	matchTheTowers(topNode);
	matchTheClusters(topNode);
}
void ShowerTowerMatching::buildTruthTowers(std::map<PHG4Particle*, Shower*> matched, std::vector<PHG4Particle*> unmatched)
{
	
	for(auto pm: matched)
	{
	}
	return;
}
void ShowerTowerMatching::matchTheTowers
{
	//matchign the meta towers to the shower 
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
 
