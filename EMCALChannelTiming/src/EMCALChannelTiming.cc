#include "EMCALChannelTiming.h"

EMCALChannelTiming::EMCALChannelTiming(const std::string &name):
SubsysReco(name)
{
	//constructor
}
EMCALChcannelTiming::process_event(PHCompositeNode* topNode)
{
	//Need to get all towers from the emcal
	//Get the towers, subdivide and send them to a helper function 
	
}
EMCALChannelTiming::SubdivideDetector(
		std::vector<tower*>* lowerEtowers, 
		std::vector<tower*>*  higherEtowers, 
		PHCompositeNode* topNode
		)
{
	auto emcaltowers = findNode::getClass<TowerInfoContainer>(topNode, emcal_towers );
	auto emcalgeom = findNode::getClass<TowerInfoContianer>(topNode, emcal_geom);

	emcalgeom->set_calorimeter_id(RawTowerDefs::CEMC);
	for(int n=0; n<(int)emcaltowers->size(); n++)
	{
		auto key = emcaltowers->encode_key(j);
		int phibin = emcaltowers->getTowerPhiBin(key);
		int etabin = emcaltowers->getTowerEtaBin(key);
		float phi = emcalgeom->get_phicenter(phibin);
		float eta = emcal->get_etacenter(etabin);
		float tow = emcaltowers->get_tower_at_channel(key);
		float E	  = tow->get_energy();
		float t   = tower->get_timing();
		float Ei  = tower->get_waveform(6);
	}
}

