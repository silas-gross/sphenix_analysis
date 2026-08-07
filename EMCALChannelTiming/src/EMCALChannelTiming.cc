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
float EMCALChannelTiming::SubdivideDetector(
		std::vector<tower*>* lowerEtowers, 
		std::vector<tower*>*  higherEtowers, 
		PHCompositeNode* topNode
		)
{
	float avgtime = 0;
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
		float t   = tow->get_timing();
		int N	  = tow->get_nsample();
		int16_t Ei=0;
		avg_time += t; 
		for(int i= 0; i< N; i++)
		{
			temp = tower->get_waveform_value(i);
			if(temp > Ei) Ei=temp;
		}
		if(Ei > 30 and E < 5000) continue;
		tower* twA = new tower(phi, eta, E, Ei, t);
		if( Ei <= 30 ) lowerEtowers->push_back(twA);
		if( E >= 5000 ) higherEtowers->pushback(twA);
	}
	avg_time = avg_time/((float) emcaltowers->size());
	return avg_time; 
}
void EMCALChannelTiming::AnaHelper(
		std::vector<tower*> subsettowers,
		std::vector<TH1F*>* 1Doutput,
		std::vector<TH2F*>* 2Doutput, 
		float TAvg
		)
{
	std::map<int, std::pair<int, float>> tavg {};
	for(
		int i = 0; 
		i< (int)1Doutput->at(1DOUTPUTHISTS::E)->getNbins(); 
		i++
	)
	{
		tavg[i]=std::make_pair(0, 0.);
	}

	for(auto tow: subsettowers)
	{
		1Doutput->at(1DOUTPUTHISTS::DELTAT)
			->Fill(tow->t - TAvg);
		1Doutput->at(1DOUTPUTHISTS::E)->
			Fill(E);
		int binN = 1Doutput->at(1DOUTPUTHISTS::E)->findBin(E)
		tavg[binN].first++;
		tavg[binN].second+=tow->t - TAvg;
	}
	for(auto m:tavg)
	{
		m.second.second= m.second.second/(float)m.second.first;
		float Ec= 1Doutput->at(1DOUTPUTHISTS::E)->getBinCenter(m.first);
		1Doutput->at(1DOUTPUTHISTS::EBART)->Fill(Ec, m.second.second);
	}
}
