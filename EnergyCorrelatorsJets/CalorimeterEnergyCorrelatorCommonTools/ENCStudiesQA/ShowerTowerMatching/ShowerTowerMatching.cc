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
	auto truthinfoContainer = findNode::getClass<PHG4TruthInfoContainer*>(topNode, "PHG4TruthInfoContainer");
	
}
void ShowerTowerMatching::getParticleShower(
		PHG4Particle* p1, 
		PHG4Hit* h1, 
		Shower* shower, 
		PHCompositeNode* topNode
		)
{
	//this is where we turn a particle into a set of Hits at the calo level
	auto 
}
int ShowerTowerMatching::End(PHCompositeNode* topNode)
{
}
 
