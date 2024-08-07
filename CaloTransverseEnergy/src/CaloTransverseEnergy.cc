#include "CaloTransverseEnergy.h"
std::vector <int> packets;
class PHCompositeNode;
class Fun4AllInputManager;
class Event;
const std::string &prdfnode="PRDF"; 
const std::string &DSTnode="DST";
const std::string &iHCALnode="HCALIN";
const std::string &oHCALnode="HCALOUT";
const std::string &EMCALnode="CEMC";
std::vector <int> n_evt_eta;
float phi_top, phi_bot;
int CaloTransverseEnergy::process_event(PHCompositeNode *topNode)
{
//	if(!ApplyCuts(e)) return 1;
	n_evt++;
	int zbin=0;
	std::vector<float> emcalenergy_vec, ihcalenergy_vec, ohcalenergy_vec; 
	if(isPRDF){
	//	if(n_evt>1000000) return 1;
		for (auto pn:packets)
		{
//			std::cout<<"Looking at packet " <<pn <<std::endl;
			Event* e= findNode::getClass<Event>(topNode, prdfnode);
			if(pn/1000 ==6 ) 
			{
			//this is the EMCal
				processPacket(pn, e, &emcalenergy_vec, false);
			}
			else if (pn/1000 == 7) 
			{
			//inner Hcal
				processPacket(pn, e, &ihcalenergy_vec, true);
			}
			else if (pn/1000 == 8)
			{	
			//outerhcal
				processPacket(pn, e, &ohcalenergy_vec, true);
			}
			else
			{
			//not a calorimeter 
				//packets.erase(&pn); //this is not quite the way to do it, but the idea is there
				continue;
			}
		}
	}
	else { // This is the DST Processor
		//
		energy=0; 
		hcalenergy=0;
		emcalenergy=0;
		energy_transverse=0;
		et_hcal=0;
		et_emcal=0;
		std::cout <<"Running on event " <<n_evt <<std::endl;
		TowerInfoContainerv2* ihe=NULL, *ohe=NULL, *eme=NULL, *ihek=NULL, *ohek=NULL, *emek=NULL;
		TowerInfoContainerv1* ihes=NULL, *ohes=NULL, *emes=NULL;// *iheks=NULL, *oheks=NULL, *emeks=NULL; 
		std::string ihcalgeom="TOWERGEOM_HCALIN", ohcalgeom="TOWERGEOM_HCALOUT", emcalgeom="TOWERGEOM_CEMC";
		if(!sim) {
			ihe=findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALIN");
			ohe=findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALOUT");
			ihek=findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_HCALIN");
			ohek=findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_HCALOUT");
			emek=findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_CEMC");
			eme=findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_CEMC");
		}
		else if(sim) {
			ihes=findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALIN");
			ohes=findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALOUT");
			//iheks=findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_HCALIN");
			//oheks=findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_HCALOUT");
			//emeks=findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_CEMC");
			emes=findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_CEMC");
		}
		RawTowerGeomContainer_Cylinderv1 *ihg=findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, ihcalgeom);
		RawTowerGeomContainer_Cylinderv1 *ohg=findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, ohcalgeom);
		RawTowerGeomContainer_Cylinderv1 *emg=findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, emcalgeom);
		float z_vtx=0;
		std::cout<<"Now finding vertex" <<std::endl; 
		GlobalVertexMap* mbdvtxmap=findNode::getClass<GlobalVertexMap>(topNode,"GlobalVertexMap");
		//get the z vertex of the event here
		if(mbdvtxmap->empty()){
//			 std::cout<<"empty mbdmap" <<std::endl;
		}
		GlobalVertex* mbdvtx=nullptr;
		if(mbdvtxmap ){ 
			for(GlobalVertexMap::ConstIter mbditer = mbdvtxmap->begin(); mbditer != mbdvtxmap->end(); ++mbditer)
			{
				mbdvtx=mbditer->second;
			}
			if(mbdvtx) z_vtx=mbdvtx->get_z();
		}
		std::cout<<"Getting Energies around vertex " <<z_vtx <<std::endl;
		int z_bin=rint(z_vtx);
		z_bin=2*(rint(z_bin)/2);
//		if(z_bin % 3 == 1) z_bin=z_bin-1;
//		if(z_bin % 3 == 2) z_bin=z_bin-2;
		if(z_bin<-20) z_bin=-20;
		else if(z_bin > 20) z_bin = 20;
		zbin=z_bin;
		std::cout<<"z bin is " <<z_bin <<std::endl;
		if(zPLTS.find(zbin) == zPLTS.end()){ 
			std::cout<<"Could not find the zbin, what is up with that?" <<std::endl;
			return 1;}
		if(!sim){
			try{	
		std::cout<<"Plots with z_bin " <<zbin <<" have high value of " <<zPLTS[zbin]->zh <<std::endl;
			phi_top=0; 
			phi_bot=0;
			processDST(eme,emek,&emcalenergy_vec, emg, false, false, RawTowerDefs::CEMC, z_vtx, zPLTS[z_bin]);
			float rat=(phi_top/3.)/(phi_bot/29.);
			zPLTS[z_bin]->em->top3_behavior->Fill(rat);
			}
		catch(std::exception& e){
			n_evt--;
			return 1;}
		try{
			phi_top=0; 
			phi_bot=0;
			processDST(ihe,ihek,&ihcalenergy_vec, ihg, true, true, RawTowerDefs::HCALIN, z_vtx, zPLTS[z_bin]);
			float rat=(phi_top/3.)/(phi_bot/29.);
			zPLTS[z_bin]->ihcal->top3_behavior->Fill(rat);
		}
		catch(std::exception& e){
			n_evt--;
			return 1;}
		try{
			phi_top=0; 
			phi_bot=0;
			processDST(ohe,ohek,&ohcalenergy_vec, ohg, true, false, RawTowerDefs::HCALOUT, z_vtx, zPLTS[z_bin]);
			float rat=(phi_top/3.)/(phi_bot/29.);
			zPLTS[z_bin]->ohcal->top3_behavior->Fill(rat);
		}
		catch(std::exception& e){
			n_evt--;
			return 1;}
		}
		else if(sim && !truth){
//			std::cout<<"Trying to run on the simulation data "<<std::endl;
			try{	
			std::cout<<"Plots with z_bin " <<zbin <<" have high value of " <<zPLTS[zbin]->zh <<std::endl;
			phi_top=0; 
			phi_bot=0;
			processDST(emes,emes,&emcalenergy_vec, emg, false, false, RawTowerDefs::CEMC, z_vtx, zPLTS[z_bin]);
			float rat=(phi_top/3.)/(phi_bot/29.);
			zPLTS[z_bin]->em->top3_behavior->Fill(rat);
		}
		catch(std::exception& e){
			n_evt--;
			return 1;}
		try{
			phi_top=0; 
			phi_bot=0;
			processDST(ihes,ihes,&ihcalenergy_vec, ihg, true, true, RawTowerDefs::HCALIN, z_vtx, zPLTS[z_bin]);
			float rat=(phi_top/3.)/(phi_bot/29.);
			zPLTS[z_bin]->ihcal->top3_behavior->Fill(rat);
		}
		catch(std::exception& e){
			n_evt--;
			return 1;}
		try{
			phi_top=0; 
			phi_bot=0;
			processDST(ohes,ohes,&ohcalenergy_vec, ohg, true, false, RawTowerDefs::HCALOUT, z_vtx, zPLTS[z_bin]);
			float rat=(phi_top/3.)/(phi_bot/29.);
			zPLTS[z_bin]->ohcal->top3_behavior->Fill(rat);
}
		catch(std::exception& e){
			n_evt--;
			return 1;}
		}
		z_vertex->Fill(z_vtx);	
//		std::cout<<"Found energy"<<std::endl;
		if(sim && truth){
			phi_top=0; 
			phi_bot=0;
			try{
				PHHepMCGenEventMap *phg=findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
				//EventHeaderv1 *evthead=findNode::getClass<EventHeaderv1>(topNode, "EventHeader");
				PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
				std::map<std::string, std::vector<float>> hep_map;
//				std::cout<<"Found the hep map"<<std::endl;
				for(PHHepMCGenEventMap::ConstIter eventIter=phg->begin(); eventIter != phg->end(); ++eventIter){
					PHHepMCGenEvent * hpev=eventIter->second;
					if(hpev && hpev->get_embedding_id() == 0){
						HepMC::GenEvent *truthevent =hpev->getEvent();		
						PHHepMCGenEvent *phe=phg->get(0);
						hep_map["VTX_Z"].push_back(phe->get_collision_vertex().z());
						if(!truthevent) return 1;
						for(HepMC::GenEvent::particle_const_iterator iter = truthevent->particles_begin(); iter != truthevent ->particles_end(); ++iter){
							if(!(*iter)->end_vertex() && (*iter)->status() ==1){
								hep_map["Energy"].push_back((*iter)->momentum().e());
								double px=(*iter)->momentum().px();
								double py=(*iter)->momentum().py();
								double pz=(*iter)->momentum().pz();
								hep_map["PT"].push_back(sqrt(px*px+py*py));
								hep_map["PZ"].push_back(pz);
								hep_map["phi"].push_back(atan2(py, px));
								hep_map["eta"].push_back(asinh(pz/hep_map["PT"].back()));
								hep_map["ID"].push_back((*iter)->pdg_id());
								hep_map["Mass"].push_back((float)(*iter)->generated_mass());
								part_mass->Fill(hep_map["Mass"].back());;
							}
						}
						break;
					}
				}
//				std::cout<<"Have the hep particle maps loaded" <<std::endl;
				if(!truthinfo) return 1;
				PHG4TruthInfoContainer::Range range=truthinfo->GetPrimaryParticleRange();
				int tpt=0;
				std::map<int, float> truth_event_e, evt_ET;
				for(int i=-100; i<=100; i++){
						if(i%10==0){
							truth_event_e[i]=0;
							evt_ET[i]=0;
						}
				}
				for(PHG4TruthInfoContainer::ConstIterator iter= range.first; iter !=range.second; ++iter)
				{
					float pmass = 0.938272;
					const PHG4Particle* truth=iter->second;
					if(!truthinfo->is_primary(truth)) continue;
					float pt=sqrt(pow(truth->get_px(),2) + pow(truth->get_py(),2));
					float pz=truth->get_pz()+hep_map["VTX_Z"].at(0);
					float E=0, ps=pow(pt,2)+pow(pz,2);
					std::vector<int>::iterator parit =std::find(baryons.begin(), baryons.end(), truth->get_pid());
					std::vector<int>::iterator aparit =std::find(baryons.begin(), baryons.end(), -truth->get_pid());
					if(parit !=baryons.end()){
						E=truth->get_e();// -sqrt(pow(truth->get_e(),2)-ps);
						if(sqrt(ps) > truth->get_e()){ std::cout<<"PS is too high Energy " <<truth->get_e() <<" to ps " <<ps <<" for particle id " <<truth->get_pid() <<std::endl;
						E=truth->get_e() - pmass;}
					}
					else if(aparit != baryons.end()){
						 E=truth->get_e();// -sqrt(pow(truth->get_e(),2)-ps)+2*pmass;
						if(sqrt(ps) > truth->get_e()){
							 std::cout<<"PS is too high Energy " <<truth->get_e() <<" to ps " <<sqrt(ps) <<" for particle id " <<truth->get_pid() <<std::endl;
							E=truth->get_e()+pmass; 
						}
					}
					else E=truth->get_e();
					float phi=atan2(truth->get_py(), truth->get_px());
					float eta=asinh(pz/pt);
					
					if(abs(eta) > 1.13 ){
						std::cout<<"Found an out of bound eta "<< eta <<std::endl; 
					//	continue;
					}
					tpt++;
//					std::cout<<"Found a particle with particle id " <<truth->get_pid() <<" with energy " <<E <<std::endl;
					int zbin=10*(hep_map["VTX_Z"].at(0)/10);
					if(zbin <-100) zbin=-100;
					if(zbin >100) zbin=100;
					zbin=0;
					//zbin=0;
					float etawidth=2.2/24., phiwidth=3.14159/32.;
					int etabin=(1.1+eta-etawidth/2.0)/etawidth, phibin=(phi+3.1415-phiwidth/2.0)/phiwidth;
					truth_event_e[zbin]+=E;
//					std::cout<<"added an energy of " <<E <<std::endl;
					try{std::cout<<szPLTS[zbin]->z_bin<<std::endl;}
					catch(std::exception& e){std::cout<<"found bad zbin " <<zbin <<std::endl;}
					szPLTS[zbin]->total->phi->Fill(phi);
					szPLTS[zbin]->total->eta->Fill(eta);
					szPLTS[zbin]->total->E_phi->Fill(phi, E);
//					std::cout<<"accessed bins of szplts" <<std::endl;
					float et=GetTransverseEnergy(E, eta);
//					float ET=(float)E/(float)cosh(eta);
				//	std::cout <<"ET is " <<et  <<"for eta " <<eta <<std::endl;
					if(etabin == 1 || etabin == 22){
						if(phibin>=14 && phibin <=19 ) phi_top+=et;
						else phi_bot+=et;
						szPLTS[0]->total->phi_projections->Fill(phibin, et/2);
					}
					evt_ET[zbin]+=et;
					et=et/2.2;
					et=et*24.0;
					szPLTS[zbin]->total->dET_eta->Fill(eta,et);
					if(phi < 5.*3.141593/8. && phi >= 7.*3.14159/16.) szPLTS[zbin]->total->dET_eta_top->Fill(eta, et);
					else if(phi < 3.*3.141593/16. && phi >=0.) szPLTS[zbin]->total->dET_eta_mid->Fill(eta, et);
					else if(phi < 13.*3.141593/8. && phi >= 23.*3.14159/16.) szPLTS[zbin]->total->dET_eta_top->Fill(eta, et);
					et=et/24.0;
					szPLTS[zbin]->total->ET_phi->Fill(phi, et);
					szPLTS[zbin]->total->ET_eta_phi->Fill(etabin, phibin, et);
					szPLTS[zbin]->total->E_eta_phi->Fill(etabin, phibin, E);
					szPLTS[zbin]->total->Hits2D->Fill(etabin, phibin);
				
					if(n_evt==1){
						std::map<std::pair<double, double>, float> c;
						c[std::make_pair(-10, -10)] =energy;
						jet_like_container[std::make_pair(eta, phi)]=c;
					}
					else if(n_evt==2){
						for(auto t:jet_like_container){
							std::pair<double, double> p=t.first;
							float DP2=phi-p.second;
							DP2=std::min(std::min(abs(DP2), (float)abs(DP2-6.283)), (float)abs(DP2-6.283));
							float R=sqrt(pow((eta-p.first),2)+pow(DP2,2));
							if(p == std::make_pair((double)eta, (double)phi)) t.second[std::make_pair(-10, -10)]+=energy;
							else if(R <= 0.4) t.second[std::make_pair((double)eta, (double)phi)]= energy;
						}
					}
					else{
						for(auto t:jet_like_container){
							std::pair<double, double> p=t.first;
							float DP2=phi-p.second;
							DP2=std::min(std::min(abs(DP2), (float)abs(DP2-6.283)), (float)abs(DP2-6.283));
							float R=sqrt(pow((eta-p.first),2)+pow(DP2,2));
							if(p == std::make_pair((double)eta, (double)phi)) t.second[std::make_pair(-10, -10)]+=energy;
							else if(R <= 0.4) t.second[std::make_pair((double)eta, (double)phi)]+=energy;
						}
					}
				}
			//	std::cout<<"Now loading the energy into the plts" <<std::endl;
				for(auto E:truth_event_e) szPLTS[0]->total->Energy->Fill(E.second);
				for(auto E:evt_ET) szPLTS[0]->total->ET->Fill(E.second);
				float rat=(phi_top/3.)/(phi_bot/29.);
				szPLTS[0]->total->top3_behavior->Fill(rat);		
				
			}
			catch(std::exception& e){ return 1;}
		}
	}
	float emcaltotal=GetTotalEnergy(emcalenergy_vec,1)/2.26; //not sure about the calibration factor, need to check
	float ihcalcalib=1, ohcalcalib=1;
//	if(isPRDF){
//		ihcalcalib=1966.26;
//		ohcalcalib=311.6;
//	}
	float ihcaltotal=GetTotalEnergy(ihcalenergy_vec,ihcalcalib)/2.2;
	float ohcaltotal=GetTotalEnergy(ohcalenergy_vec,ohcalcalib)/2.2;
	std::cout<<"Adding the data to the histograms" <<std::endl;
	zPLTS[zbin]->em->ET->Fill(emcaltotal*2.26);
	zPLTS[zbin]->em->dET->Fill(emcaltotal);
	zPLTS[zbin]->ihcal->ET->Fill(ihcaltotal*2.2);
	zPLTS[zbin]->ihcal->dET->Fill(ihcaltotal);
	zPLTS[zbin]->ohcal->ET->Fill(ohcaltotal*2.2);
	zPLTS[zbin]->ohcal->dET->Fill(ohcaltotal);
	evt_data.emcal_et=emcaltotal;
	evt_data.ohcal_et=ohcaltotal;
	evt_data.ihcal_et=ihcaltotal;
	//float total=emcaltotal+ihcaltotal+ohcaltotal;
//	if(!truth) std::cout<<"Base level sum is " <<total <<std::endl;
	std::map<std::string, float> data {{"EMCAL", emcaltotal}, {"OHCAL", ohcaltotal}, {"IHCAL", ihcaltotal}};
	zPLTS[zbin]->event_data.push_back(data);
//	PLTS.event_data.push_back(std::make_pair("OHCAL", ohcaltotal));
//	PLTS.event_data.push_back(std::make_pair("IHCAL", ihcaltotal));
//	ETOTAL->Fill(emcaltotal+ihcaltotal+ohcaltotal);
	datatree->Fill();
	evt_data.emcal_et=0;
	evt_data.ihcal_et=0;
	evt_data.ohcal_et=0;
	evt_data.emcal_tower_et.clear();
	evt_data.ihcal_tower_et.clear();
	evt_data.ohcal_tower_et.clear();
	
	ihcalenergy_vec.clear();
	ohcalenergy_vec.clear();
	emcalenergy_vec.clear();
	return 1;
}
void CaloTransverseEnergy::processDST(TowerInfoContainerv1* calo_event, TowerInfoContainerv1* calo_key, std::vector<float>* energies, RawTowerGeomContainer_Cylinderv1* geom, bool hcalorem, bool inner, RawTowerDefs::CalorimeterId caloid, float z_vtx, plots* PLTS)
{
	int maxphi=0, maxeta=0;
	geom->set_calorimeter_id(caloid);
	std::cout<<"Have now entered the DST processor and set the caloid on event " <<n_evt <<std::endl;
	//This processes all events in the DST in the processor 
	int ntowers=calo_event->size();
	int n_live_towers=0, n_time=0;
//	std::cout<<"N Towers: " <<ntowers <<std::endl;
	try{
	for(int i =0;i<ntowers; i++ )
	{
		try{calo_event->get_tower_at_channel(i);}
		catch(std::exception& e){
			std::cout<<"could not find the tower "<<std::endl;
			continue;}
		auto tower=calo_event->get_tower_at_channel(i);
		//TowerInfov1* tower=calo_event->get_tower_at_channel(i);
		float energy1=tower->get_energy();
		if(hcalorem && inner) PLTS->ihcal->E_f->Fill(energy1);
		else if(hcalorem) PLTS->ohcal->E_f->Fill(energy1);
		else PLTS->em->E_f->Fill(energy1);
		PLTS->total->E_f->Fill(energy1);
		if(energy1<0) continue;
		
		if(hcalorem && inner) PLTS->ihcal->E_m->Fill(energy1);
		else if(hcalorem) PLTS->ohcal->E_m->Fill(energy1);
		else PLTS->em->E_m->Fill(energy1);
		PLTS->total->E_m->Fill(energy1);
//		if(hcalorem) std::cout<<"energy is " <<energy1 <<std::endl;
		float time=tower->get_time();
		if(inner){
		//	 if(time<-2 || time > 2) energy1=-10;
			PLTS->ihcal->timing->Fill(time);
		} //using a specific timing for the relevant runs I'm using, should do the systematic from cdb when available 
		if (hcalorem && ! inner){ 
		//		if (time<-2 || time >2) energy1=-10;
				PLTS->ohcal->timing->Fill(time);
		}
		else{
		//	 if(time < -1 || time > 1 ) energy1=-10; 
			PLTS->em->timing->Fill(time);
		}
		PLTS->total->timing->Fill(time);
		try{calo_key->encode_key(i);}
		catch(std::exception& e){
			std::cout<<"could not find the calokey"<<std::endl;
			continue;}
		int key=calo_key->encode_key(i);
		//if(hcalorem) std::cout<<"encoding key " <<key <<std::endl;
		int phibin=calo_key->getTowerPhiBin(key);
		int etabin=calo_key->getTowerEtaBin(key);
		float radius, z;
		radius=geom->get_radius();
		n_time++;
		//std::cout<<"Tower is at radius " <<radius <<std::endl;
		//RawTowerGeom *towergeom=geom->get_tower_geometry(key);
		//assert(towergeom);
		if(energy1 <= -5){continue;}
		if(hcalorem && inner) PLTS->ihcal->E_s->Fill(energy1);
		else if(hcalorem) PLTS->ohcal->E_s->Fill(energy1);
		else PLTS->em->E_s->Fill(energy1);
		PLTS->total->E_s->Fill(energy1);
		n_live_towers++;
		if(etabin>maxeta) maxeta=etabin;
		if(phibin>maxphi) maxphi=phibin;
	//	if(inner)std::cout<<"Have figures stuff out" <<std::endl;
		double eta=geom->get_etacenter(etabin);
		double phi=geom->get_phicenter(phibin);
		std::pair<double, double> etabounds=geom->get_etabounds(etabin);
		std::pair <double, double> phibounds=geom->get_phibounds(phibin); 
	//	double phi_width=abs(phibounds.first - phibounds.second);
		if(!hcalorem){
			 eta=1.134/96*(1+2*etabin)-1.134;
			etabounds.first=eta-1.134/96.;
			etabounds.second=eta+1.134/96.;
		}
		z=radius*sinh(eta);	
		z+=z_vtx;
//		if(!hcalorem) std::cout<<"Distance of " <<z<<std::endl;
		eta=asinh(z/radius);
		etabounds.first= asinh(z_vtx/radius+sinh(etabounds.first));
		etabounds.second=asinh(radius*z_vtx+sinh(etabounds.second));
		if(hcalorem){
			if(inner && PLTS->ihcal->counter==false)
			{
				std::pair<int, int> iedge {etabin, phibin}; 
				PLTS->ihcal->tower_edges[iedge]=std::make_pair(etabounds, phibounds);
				
			}
			else if (!inner && PLTS->ohcal->counter==false){
				std::pair<int, int> iedge {etabin, phibin};
				PLTS->ohcal->tower_edges[iedge]=std::make_pair(etabounds, phibounds);
			}
		}
		if(!hcalorem && PLTS->ihcal->counter && PLTS->ohcal->counter){
			if(PLTS->ihcal->filled_em_cal==false || PLTS->ohcal->filled_em_cal==false) countTheTowers(PLTS, etabounds.first,etabounds.second, phibounds.first, phibounds.second);
		}
		if(hcalorem && !inner){
		if(n_evt==1){
			std::map<std::pair<double, double>, float> c;
			c[std::make_pair(-10, -10)] =energy1;
			jet_like_container[std::make_pair(eta, phi)]=c;
		}
		else if(n_evt==2){
			for(auto t:jet_like_container){
				std::pair<double, double> p=t.first;
				float DP2=phi-p.second;
				DP2=std::min(std::min(abs(DP2), (float)abs(DP2-6.283)), (float)abs(DP2-6.283));
				float R=sqrt(pow((eta-p.first),2)+pow(DP2,2));
				if(p == std::make_pair(eta, phi)) t.second[std::make_pair(-10, -10)]+=energy1;
				else if(R <= 0.4) t.second[std::make_pair(eta, phi)]=energy1;
			}
		}
		else{
			for(auto t:jet_like_container){
				std::pair<double, double> p=t.first;
				float DP2=phi-p.second;
				DP2=std::min(std::min(abs(DP2), (float)abs(DP2-6.283)), (float)abs(DP2-6.283));
				float R=sqrt(pow((eta-p.first),2)+pow(DP2,2));
				if(p == std::make_pair(eta, phi)) t.second[std::make_pair(-10, -10)]+=energy1;
				else if(R <= 0.4) t.second[std::make_pair(eta, phi)]+=energy1;
				}
			}
		}
		double eta_width=abs(etabounds.first - etabounds.second);
		if(!hcalorem){
			 //if( n_evt == 600 && etabin_em->GetBinContent(etabin) <= 0.001 )  etabin_em->SetBinContent(etabin, eta);
			 //if(n_evt == 600 && phibin_em->GetBinContent(phibin) <= 0.001) phibin_em->SetBinContent(phibin, phi_width);
//			std::cout<<"Trying to put in the z vertex for the emcal" <<z_vtx <<std::endl;
			try{
				PLTS->em->z_val->Fill(z_vtx);
//				std::cout<<"Does this work?" <<std::endl;
			}
			catch(std::exception& e){ std::cout<<"cant find zvtx " <<std::endl; 
			}
		} 
		else{
	//		std::cout<<"Hi anyone here?" <<std::endl;
			if(inner) PLTS->ihcal->z_val->Fill(z_vtx);
			else PLTS->ohcal->z_val->Fill(z_vtx);
	//		std::cout<<"putting in the vertex for hcal" <<std::endl;
		}
		float et=GetTransverseEnergy((float)energy1, eta);
		float et_div=energy1/cosh(eta);
		if(eta_width>0) et_div=et/(float)eta_width;
		energies->push_back(GetTransverseEnergy(energy1, eta));
		if(!hcalorem){
			//float rat=PLTS->em->dET_eta->GetNbinsX();
			float rat=4.0;
			rat=et_div*96.0/rat;
			emcalenergy+=energy1;
			//float etamin=PLTS->em->dET_eta->GetBinCenter(1), etamax=PLTS->em->dET_eta->GetBinCenter(96); 
			//if(etamin < eta || etamax > eta) std::cout<<"Out of bounds with eta : " <<eta <<" bounds are " <<etamin << " to " <<etamax<<std::endl;
//			if(et_div > 0 ) std::cout<<"Energy " <<et_div <<" with eta " <<eta <<std::endl;
			PLTS->em->phi->Fill(phibin);
			PLTS->em->eta->Fill(etabin);
			PLTS->em->E_phi->Fill(phi, energy1);
			PLTS->em->ET_phi->Fill(phi, et);
			PLTS->em->dET_eta->Fill(eta,(float)et_div);
			if(phi < 5.*3.141593/8. && phi >= 7.*3.14159/16.) PLTS->em->dET_eta_top->Fill(eta, et_div);
			else if(phi < 3.*3.141593/16. && phi >=0.) PLTS->em->dET_eta_mid->Fill(eta, et_div);
			else if(phi < 13.*3.141593/8. && phi >= 23.*3.14159/16.) PLTS->em->dET_eta_top->Fill(eta, et_div);
		
//			std::cout<<"eta " <<eta <<std::endl;
			//if(et_div >0 ) std::cout<<"filled data in to the det/eta of" <<et_div <<" and readsback as " <<PLTS->em->dET_eta->GetBinContent(PLTS->em->dET_eta->FindBin(eta)) << "\n with sample value " <<PLTS->em->dET_eta->FindBin(eta) <<" not in " <<etamin <<" , to " <<etamax <<std::endl;
			PLTS->em->ET_eta_phi->Fill(eta, phi, et_div);
			PLTS->em->Hits2D->Fill(etabin, phibin);
			PLTS->em->ET_z_eta->Fill(z_vtx, eta, et_div);
			PLTS->em->ET_z->Fill(z, et);
			PLTS->em->Eta_width->Fill(etabin, eta_width);
			EMCALE->Fill(energy1);
			if(etabin==1 || etabin==254){
				 PLTS->em->phi_projections->Fill(phibin, et/2);
				if(phi>= 7*3.14159/16. && phi< 31.4159/16.) phi_top+=et;
				else phi_bot+=et;
			}
			try{evt_data.emcal_tower_et[eta]+=et_div;}
			catch(std::exception& e){evt_data.emcal_tower_et[eta]=et_div;}
				eteeta[eta].push_back(et);
		
				etephi[phi].push_back(et);
				EtaPhi->Fill(eta, phi,et, et);	
				
		}
		if(hcalorem){
			hcalenergy+=energy1;
//			std::cout<<"Eta val of " <<etabin*12./1.1 <<" is now " <<eta <<" in the hcal" <<std::endl;
			if(inner){
				PLTS->ihcal->phi->Fill(phibin);
				PLTS->ihcal->eta->Fill(etabin);
				PLTS->ihcal->E_phi->Fill(phi, energy1);
				PLTS->ihcal->ET_phi->Fill(phi, et);
			//	if(eta< PLTS->ihcal->etamin || eta > PLTS->ihcal->etamax) std::cout<<"Non-good eta " <<eta <<std::endl;
				PLTS->ihcal->dET_eta->Fill(eta, (float)et_div);
				if(phi < 5.*3.141593/8. && phi >= 7.*3.14159/16.) PLTS->ihcal->dET_eta_top->Fill(eta, et_div);
				else if(phi < 3.*3.141593/16. && phi >=0.) PLTS->ihcal->dET_eta_mid->Fill(eta, et_div);
				else if(phi < 13.*3.141593/8. && phi >= 23.*3.14159/16.) PLTS->ihcal->dET_eta_top->Fill(eta, et_div);
				PLTS->ihcal->ET_eta_phi->Fill(eta, phi, et_div);
				PLTS->ihcal->Hits2D->Fill(etabin, phibin);
				PLTS->ihcal->Eta_width->Fill(etabin, eta_width);
				IHCALE->Fill(energy1);
				if(etabin==1 || etabin==22){ PLTS->ihcal->phi_projections->Fill(phibin, et/2);
				if(phi>= 7*3.14159/16. && phi< 31.4159/16.) phi_top+=et;
				else phi_bot+=et;
			}
				try{evt_data.ihcal_tower_et[eta]+=et_div;}
				catch(std::exception& e){evt_data.ihcal_tower_et[eta]=et_div;}
			}
					
			else{
				try{evt_data.ohcal_tower_et[eta]+=et_div;}
				catch(std::exception& e){evt_data.ohcal_tower_et[eta]=et_div;}
				PLTS->ohcal->phi->Fill(phibin);
				PLTS->ohcal->eta->Fill(etabin);
				PLTS->ohcal->E_phi->Fill(phi, et);
				PLTS->ohcal->dET_eta->Fill(eta, (float)et_div);
				if(phi < 5.*3.141593/8. && phi >= 7.*3.14159/16.) PLTS->ohcal->dET_eta_top->Fill(eta, et_div);
				else if(phi < 3.*3.141593/16. && phi >=0.) PLTS->ohcal->dET_eta_mid->Fill(eta, et_div);
				else if(phi < 13.*3.141593/8. && phi >= 23.*3.14159/16.) PLTS->ohcal->dET_eta_top->Fill(eta, et_div);
//				std::cout<<"eta " <<eta <<std::endl;
				PLTS->ohcal->ET_eta_phi->Fill(eta, phi, et_div);
				PLTS->ohcal->Hits2D->Fill(etabin, phibin);
				PLTS->ohcal->Eta_width->Fill(etabin, eta_width);
				if(etabin==1 || etabin==22){ PLTS->ohcal->phi_projections->Fill(phibin, et/2);
				if(phi>= 7*3.14159/16. && phi< 31.4159/16.) phi_top+=et;
				else phi_bot+=et;
			}
				OHCALE->Fill(energy1);
			}
			etheta[eta].push_back(et);
			EtaPhi->Fill(eta, phi, et, 0, et);
			ethphi[phi].push_back(et);
			}
		
		//std::cout<<"Have added the energy " <<energy1 <<std::endl;
		
	}
	}
	catch (std::exception& e){std::cout<<"Exception found " <<e.what() <<std::endl;}
	//if(n_evt==10) std::cout<<"Calorimeter " <<caloid <<" max: etabin: " <<maxeta <<", phibin: " <<maxphi <<std::endl;
//	std::cout<<"Had " <<ntowers <<" total towers, but only " <<n_live_towers <<" passed the energy cut and " <<n_time <<" even got to the  cut" <<std::endl;
	if (!hcalorem && PLTS->ohcal->filled_em_cal && PLTS->ohcal->counter){
		for(auto a:PLTS->ohcal->n_emcaltowers){
			PLTS->ohcal->EMCAL_towering->Fill(a.second);
		}
	}
	if (!hcalorem && PLTS->ihcal->filled_em_cal && PLTS->ihcal->counter){
		for (auto a:PLTS->ihcal->n_emcaltowers){
				PLTS->ihcal->EMCAL_towering->Fill(a.second);
			}
		}
	
	if(hcalorem){
		if(inner && PLTS->ihcal->tower_edges.size() > 0) PLTS->ihcal->counter=true;
		if(!inner && PLTS->ohcal->tower_edges.size() > 0) PLTS->ohcal->counter=true;
	}
	if(!hcalorem){
		if(PLTS->ihcal->counter && PLTS->ihcal->n_emcaltowers.size() > 0 ) PLTS->ihcal->filled_em_cal=true;
		if(PLTS->ohcal->counter && PLTS->ohcal->n_emcaltowers.size() > 0 ) PLTS->ohcal->filled_em_cal=true;
		std::cout<<"Have done the first filling of the thing" <<std::endl;
	}
}
void CaloTransverseEnergy::processDST(TowerInfoContainerv2* calo_event, TowerInfoContainerv2* calo_key, std::vector<float>* energies, RawTowerGeomContainer_Cylinderv1* geom, bool hcalorem, bool inner, RawTowerDefs::CalorimeterId caloid, float z_vtx, plots* PLTS)
{
	int maxphi=0, maxeta=0;
	geom->set_calorimeter_id(caloid);
	std::cout<<"Have now entered the DST processor and set the caloid on event " <<n_evt <<std::endl;
	//This processes all events in the DST in the processor 
	int ntowers=calo_event->size();
	int n_live_towers=0, n_time=0;
	std::cout<<"N Towers: " <<ntowers <<std::endl;
	try{
	for(int i =0;i<ntowers; i++ )
	{
		try{calo_event->get_tower_at_channel(i);}
		catch(std::exception& e){
//			std::cout<<"could not find the tower "<<std::endl;
			continue;}
		auto tower=calo_event->get_tower_at_channel(i);
		//TowerInfov1* tower=calo_event->get_tower_at_channel(i);
		float energy1=tower->get_energy();
		if(hcalorem && inner) PLTS->ihcal->E_f->Fill(energy1);
		else if(hcalorem) PLTS->ohcal->E_f->Fill(energy1);
		else PLTS->em->E_f->Fill(energy1);
		PLTS->total->E_f->Fill(energy1);
		if(energy1<0) continue;
		
		if(hcalorem && inner) PLTS->ihcal->E_m->Fill(energy1);
		else if(hcalorem) PLTS->ohcal->E_m->Fill(energy1);
		else PLTS->em->E_m->Fill(energy1);
		PLTS->total->E_m->Fill(energy1);
//		if(!hcalorem) std::cout<<"energy is " <<energy1 <<std::endl;
		float time=tower->get_time();
		if(inner){
//			 if(time<-2 || time > 2) energy1=-10;
			PLTS->ihcal->timing->Fill(time);
		} //using a specific timing for the relevant runs I'm using, should do the systematic from cdb when available 
		if (hcalorem && ! inner){ 
//				if (time<-2 || time >2) energy1=-10;
				PLTS->ohcal->timing->Fill(time);
		}
		else{
//			 if(time < -1 || time > 1 ) energy1=-10; 
			PLTS->em->timing->Fill(time);
		}
		PLTS->total->timing->Fill(time);
		try{calo_key->encode_key(i);}
		catch(std::exception& e){
//			std::cout<<"could not find the calokey"<<std::endl;
			continue;}
		int key=calo_key->encode_key(i);
		//if(hcalorem) std::cout<<"encoding key " <<key <<std::endl;
		int phibin=calo_key->getTowerPhiBin(key);
		int etabin=calo_key->getTowerEtaBin(key);
		bool status=tower->get_isGood();
		//if(!status) continue;
		float radius, z;
		radius=geom->get_radius(); 
		n_time++;
//		if(!hcalorem) std::cout<<"Tower is at radius " <<radius <<std::endl;
		//RawTowerGeom *towergeom=geom->get_tower_geometry(key);
		//assert(towergeom);
		//if(energy1 <= -5){continue;}
		if(hcalorem && inner) PLTS->ihcal->E_s->Fill(energy1);
		else if(hcalorem) PLTS->ohcal->E_s->Fill(energy1);
		else PLTS->em->E_s->Fill(energy1);
		PLTS->total->E_s->Fill(energy1);
		n_live_towers++;
		if(etabin>maxeta) maxeta=etabin;
		if(phibin>maxphi) maxphi=phibin;
	//	if(inner)std::cout<<"Have figures stuff out" <<std::endl;
		double eta=geom->get_etacenter(etabin);
		double phi=geom->get_phicenter(phibin);
		std::pair<double, double> etabounds=geom->get_etabounds(etabin);
		std::pair <double, double> phibounds=geom->get_phibounds(phibin); 
	//	double phi_width=abs(phibounds.first - phibounds.second);
		if(!hcalorem){
			 eta=1.134/96*(1+2*etabin)-1.134;
			etabounds.first=eta-1.134/96.;
			etabounds.second=eta+1.134/96.;
		}
		z=radius*sinh(eta);	
		z+=z_vtx;
//		if(!hcalorem) std::cout<<"Distance of " <<z<<std::endl;
		eta=asinh(z/radius);
		etabounds.first= asinh(z_vtx/radius+sinh(etabounds.first));
		etabounds.second=asinh(radius*z_vtx+sinh(etabounds.second));
		if(hcalorem ){
		/*	try{
				std::cout<<"The calorimeter load in has size " <<oh_em.empty() <<std::endl;
			}
			catch(const std::exception& e1){ std::throw_with_nested(std::runtime_error("calo map failed")); 
			std::cerr<<"Found an exception when trying to access the container "<<e1.what() <<std::endl;}
			std::cout<<"ok moving forwards" <<std::endl;*/
	//	std::cout<<"Trying to fill map info. OHCAL map has size " <<oh_em.size() <<" IHmap has size " <<ih_em.size() <<std::endl;	
			if(inner && PLTS->ihcal->counter==false)
			{
				std::pair<int, int> iedge {etabin, phibin}; 
				PLTS->ihcal->tower_edges[iedge]=std::make_pair(etabounds, phibounds);
				
			}
			else if (!inner && PLTS->ohcal->counter==false){
				std::pair<int, int> iedge {etabin, phibin};
				PLTS->ohcal->tower_edges[iedge]=std::make_pair(etabounds, phibounds);
			}
		}
		if(!hcalorem && PLTS->ihcal->counter && PLTS->ohcal->counter){
			if(PLTS->ihcal->filled_em_cal==false || PLTS->ohcal->filled_em_cal==false) countTheTowers(PLTS, etabounds.first,etabounds.second, phibounds.first, phibounds.second);
		}
		if(!status) continue;
		if(hcalorem && !inner){
		if(jet_like_container.size() < 10 ){
			std::cout<<"trying the seeds"<<std::endl;
			std::map<std::pair<double, double>, float> c;
			c[std::make_pair(-10.0, -10.0)] =energy1;
			jet_like_container[std::make_pair(eta, phi)]=c;
			std::cout<<"have a jet cont with n seeds " <<jet_like_container.size()<<std::endl;
			}
		else if( jet_like_container[std::make_pair(eta,phi)].size() == 1 ){
			for(auto t:jet_like_container){
				std::pair<double, double> p=t.first;
				float DP2=phi-p.second;
				DP2=std::min(std::min(abs(DP2), (float)abs(DP2-6.283)), (float)abs(DP2-6.283));
				float R=sqrt(pow((eta-p.first),2.0)+pow(DP2,2.0));
				if(R < 0.8) std::cout<<"have the pair " <<p.first <<" eta, phi " <<p.second <<" with R " <<R <<std::endl;
				if(R<0.8) std::cout<<"The Key pair is " <<p.first <<p.second <<" compared to " <<eta <<phi<<std::endl; 
				if(p == std::make_pair(eta, phi)) t.second[std::make_pair(-10.0, -10.0)]+=energy1;
				else if(R <= 0.4){ t.second[std::make_pair(eta, phi)]= energy1;
				std::cout<<"added a tower in the jet seperated by R="<<R<<std::endl;}
				}
			}
		else{
			for(auto t:jet_like_container){
				std::pair<double, double> p=t.first;
				if(p == std::make_pair(eta, phi)) t.second[std::make_pair(-10, -10)]+=energy1;
				else if(t.second.find(std::make_pair(eta,phi)) != t.second.end() )  t.second[std::make_pair(eta, phi)]+=energy1;
				}
			}
		}
		double eta_width=abs(etabounds.first - etabounds.second);
		if(!hcalorem){
			 //if( n_evt == 600 && etabin_em->GetBinContent(etabin) <= 0.001 )  etabin_em->SetBinContent(etabin, eta);
			 //if(n_evt == 600 && phibin_em->GetBinContent(phibin) <= 0.001) phibin_em->SetBinContent(phibin, phi_width);
//			std::cout<<"Trying to put in the z vertex for the emcal" <<z_vtx <<std::endl;
			try{
				PLTS->em->z_val->Fill(z_vtx);
//				std::cout<<"Does this work?" <<std::endl;
			}
			catch(std::exception& e){ std::cout<<"cant find zvtx " <<std::endl; 
			}
		} 
		else{
	//		std::cout<<"Hi anyone here?" <<std::endl;
			if(inner) PLTS->ihcal->z_val->Fill(z_vtx);
			else PLTS->ohcal->z_val->Fill(z_vtx);
	//		std::cout<<"putting in the vertex for hcal" <<std::endl;
		}
		float et=GetTransverseEnergy(energy1, eta);
		float et_div=et;
		if(eta_width>0) et_div=et/(float)eta_width;
		//std::cout<<"Eta value is " <<eta <<std::endl;
		energies->push_back(GetTransverseEnergy(energy1, eta));
		if(!hcalorem){
			//float rat=PLTS->em->dET_eta->GetNbinsX();
			float rat=4.0;
			rat=et_div*96.0/rat;
			emcalenergy+=energy1;
	//		float etamin=PLTS->em->etamin, etamax=PLTS->em->etamax; 
	//		if(etamin < eta || etamax > eta) std::cout<<"Out of bounds with eta : " <<eta <<" bounds are " <<etamin << " to " <<etamax<<std::endl;
			PLTS->em->phi->Fill(phibin);
			PLTS->em->eta->Fill(etabin);
			PLTS->em->E_phi->Fill(phi, energy1);
			if(phi < 5.*3.141593/8. && phi >= 7.*3.14159/16.) PLTS->em->dET_eta_top->Fill(eta, et_div);
			else if(phi < 3.*3.141593/16. && phi >=0.) PLTS->em->dET_eta_mid->Fill(eta, et_div);
			else if(phi < 13.*3.141593/8. && phi >= 23.*3.14159/16.) PLTS->em->dET_eta_top->Fill(eta, et_div);
			PLTS->em->ET_phi->Fill(phi, et);
			PLTS->em->dET_eta->Fill(eta,rat);
			PLTS->em->ET_eta_phi->Fill(eta, phi, et_div);
			PLTS->em->Hits2D->Fill(etabin, phibin);
			PLTS->em->ET_z_eta->Fill(z_vtx, eta, et_div);
			PLTS->em->ET_z->Fill(z, et);
			PLTS->em->Eta_width->Fill(etabin, eta_width);
			EMCALE->Fill(energy1);
			if(etabin==1 || etabin==254){ PLTS->em->phi_projections->Fill(phibin, et/2);
				if(phi>= 7*3.14159/16. && phi< 31.4159/16.) phi_top+=et;
				else phi_bot+=et;
			}
			try{evt_data.emcal_tower_et[eta]+=et_div;}
			catch(std::exception& e){evt_data.emcal_tower_et[eta]=et_div;}
				eteeta[eta].push_back(et);
		
				etephi[phi].push_back(et);
				EtaPhi->Fill(eta, phi,et, et);	
				
		}
		if(hcalorem){
			hcalenergy+=energy1;
//			std::cout<<"Eta val of " <<etabin*12./1.1 <<" is now " <<eta <<" in the hcal" <<std::endl;
			if(inner){
				PLTS->ihcal->phi->Fill(phibin);
				PLTS->ihcal->eta->Fill(etabin);
				PLTS->ihcal->E_phi->Fill(phi, energy1);
				PLTS->ihcal->ET_phi->Fill(phi, et);
		//		if(eta< PLTS->ihcal->etamin || eta > PLTS->ihcal->etamax) std::cout<<"Non-good eta " <<eta <<std::endl;
				PLTS->ihcal->dET_eta->Fill(eta, et_div);
				if(phi < 5.*3.141593/8. && phi >= 7.*3.14159/16.) PLTS->ihcal->dET_eta_top->Fill(eta, et_div);
				else if(phi < 3.*3.141593/16. && phi >=0.) PLTS->ihcal->dET_eta_mid->Fill(eta, et_div);
				else if(phi < 13.*3.141593/8. && phi >= 23.*3.14159/16.) PLTS->ihcal->dET_eta_top->Fill(eta, et_div);
				PLTS->ihcal->ET_eta_phi->Fill(eta, phi, et_div);
				PLTS->ihcal->Hits2D->Fill(etabin, phibin);
				PLTS->ihcal->Eta_width->Fill(etabin, eta_width);
				IHCALE->Fill(energy1);
				if(etabin==1 || etabin==22){ PLTS->ihcal->phi_projections->Fill(phibin, et/2);
				if(phi>= 7*3.14159/16. && phi< 31.4159/16.) phi_top+=et;
				else phi_bot+=et;
			}
				try{evt_data.ihcal_tower_et[eta]+=et_div;}
				catch(std::exception& e){evt_data.ihcal_tower_et[eta]=et_div;}
			}
					
			else{
				try{evt_data.ohcal_tower_et[eta]+=et_div;}
				catch(std::exception& e){evt_data.ohcal_tower_et[eta]=et_div;}
				PLTS->ohcal->phi->Fill(phibin);
				PLTS->ohcal->eta->Fill(etabin);
				PLTS->ohcal->E_phi->Fill(phi, et);
				PLTS->ohcal->dET_eta->Fill(eta, et_div);
				if(phi < 5.*3.141593/8. && phi >= 7.*3.14159/16.) PLTS->ohcal->dET_eta_top->Fill(eta, et_div);
				else if(phi < 3.*3.141593/16. && phi >=0.) PLTS->ohcal->dET_eta_mid->Fill(eta, et_div);
				else if(phi < 13.*3.141593/8. && phi >= 23.*3.14159/16.) PLTS->ohcal->dET_eta_top->Fill(eta, et_div);
				PLTS->ohcal->ET_eta_phi->Fill(eta, phi, et_div);
				PLTS->ohcal->Hits2D->Fill(etabin, phibin);
				PLTS->ohcal->Eta_width->Fill(etabin, eta_width);
				OHCALE->Fill(energy1);
				if(etabin==1 || etabin==22){ PLTS->ohcal->phi_projections->Fill(phibin, et/2);
				if(phi>= 7*3.14159/16. && phi< 31.4159/16.) phi_top+=et;
				else phi_bot+=et;
			}
			}
			etheta[eta].push_back(et);
			EtaPhi->Fill(eta, phi, et, 0, et);
			ethphi[phi].push_back(et);
			}
		
		//std::cout<<"Have added the energy " <<energy1 <<std::endl;
		
	}
	}
	catch (std::exception& e){std::cerr<<"Exception found " <<e.what() <<std::endl;}
	if (jet_like_container.size() > 10 &&  hcalorem && !inner ) std::cout<<"The number of towers in the jet like thing is " <<jet_like_container.begin()->second.size() <<std::endl;
	//if(n_evt==10) std::cout<<"Calorimeter " <<caloid <<" max: etabin: " <<maxeta <<", phibin: " <<maxphi <<std::endl;
//	std::cout<<"Had " <<ntowers <<" total towers, but only " <<n_live_towers <<" passed the energy cut and " <<n_time <<" even got to the  cut" <<std::endl;
	if (!hcalorem && PLTS->ohcal->filled_em_cal && PLTS->ohcal->counter){
		for(auto a:PLTS->ohcal->n_emcaltowers){
			PLTS->ohcal->EMCAL_towering->Fill(a.second);
		}
	}
	if (!hcalorem && PLTS->ihcal->filled_em_cal && PLTS->ihcal->counter){
		for (auto a:PLTS->ihcal->n_emcaltowers){
				PLTS->ihcal->EMCAL_towering->Fill(a.second);
			}
		}
	
	if(hcalorem){
		if(inner && PLTS->ihcal->tower_edges.size() > 0) PLTS->ihcal->counter=true;
		if(!inner && PLTS->ohcal->tower_edges.size() > 0) PLTS->ohcal->counter=true;
	}
	if(!hcalorem){
		if(PLTS->ihcal->counter && PLTS->ihcal->n_emcaltowers.size() > 0 ) PLTS->ihcal->filled_em_cal=true;
		if(PLTS->ohcal->counter && PLTS->ohcal->n_emcaltowers.size() > 0 ) PLTS->ohcal->filled_em_cal=true;
		std::cout<<"Have done the first filling of the thing" <<std::endl;
	}
	
}
void CaloTransverseEnergy::processPacket(int packet, Event * e, std::vector<float>* energy, bool HorE)
{
	try{
		e->getPacket(packet);
	}
	catch(std::exception* e) {return;} 
	Packet *p= e->getPacket(packet);
	if(!p) return;
//	std::cout<<"Have the packet with " <<p->iValue(0, "CHANNELS") <<" many channels"<<std::endl; 
	for(int c=0; c<p->iValue(0, "CHANNELS"); c++)
	{
		int etabin=0;
		float eta=0, baseline=0, en=0;
		if(HorE){
		 	 int cable_channel=(c%16), cable_number=c/64;
			 cable_channel=cable_channel/2; //doubling of eta to get extra phi resolution
			 etabin=cable_channel+8*cable_number; //24 eta bins
			 eta= (float) etabin/12.0-1.1; //HCal
		}
		else{
			//int cable_channel=(c%16), cable_number=c/64; 
			eta=(((c%16)/2+c/4)-48)/48+1/96; //EMCal
		}
	//	std::cout<<"Eta value is " <<eta <<std::endl;
		for(int s=0; s<p->iValue(c, "SAMPLES"); s++)
		{
			if(s<3) baseline+=p->iValue(s,c); 
			if(s==3) baseline=baseline/3; //first 3 for baseline
			float val=p->iValue(s,c);
			if(val>en) en=val;
		}
		en=en-baseline;
		if(packet/1000 == 7 ) en=en/1966.26;
		if(packet/1000 == 8 ) en=en/311.6;
		if(energy->size() >= 2) if(GetTransverseEnergy(en,eta)>2*energy->at(energy->size()-1)) continue; 
		if(en <=0 ) continue;
		//control for high channels and 0 energy conditions
		energy->push_back(GetTransverseEnergy(en,eta));
	//	std::cout<<"Energy is " <<en <<std::endl;
		//For the hcal get the phi distribution
		//take packet#, each packet has 8 phi bins
		//
		int phibin=(c%64)/16;
		phibin=phibin*2+c%2;
		phibin+=(packet%10)*8;
		float phi=phibin*PI/16;	
		float phval=GetTransverseEnergy(en,eta);
		//PhiD->Fill(phi, phval);	
	//	EtaD->Fill(etabin, phval*24.0/2.2);
		tep->Fill(etabin, phibin, phval);
		teps->Fill(etabin, phibin, 1);
		if(packet/1000 == 7) 
		{
			ihPhiD->Fill(phi, phval);
			ihEtaD->Fill(etabin,phval*24/2.2);
			ihep->Fill(etabin, phibin, phval);
			iheps->Fill(etabin, phibin, 1);
		}
		if(packet/1000 == 8) 
		{
			ohPhiD->Fill(phi, phval);
			ohEtaD->Fill(etabin,phval*24/2.2);
			ohep->Fill(etabin, phibin, phval);
			oheps->Fill(etabin, phibin, 1);
		}
			
		//n_evt_eta.at(etabin)++;
	}

}
float CaloTransverseEnergy::GetTransverseEnergy(float energy, float eta)
{
	float et=energy/cosh(eta);
	if(et <0) std::cout<<"transverse energy is less than zero " <<et <<" energy is " <<energy <<" cosh eta is " <<cosh(eta) <<std::endl;
	//if (eta<0.55 && eta >0.45 && et==0) std::cout<<"Zero et with energy: " <<energy <<" eta: " <<eta <<" theta: "<<sin(2*atan(exp(-eta)))<<std::endl;
	return et;
}
float CaloTransverseEnergy::GetTotalEnergy(std::vector<float> caloenergy, float calib)
{
	//This is really just a fancy sumation method
//	std::cout<<"The total energy has " <<caloenergy.size() <<" entries to sum" <<std::endl;	
	float total_energy=0;
	for(auto e:caloenergy)
	{
		total_energy+=e/calib; //right now this is summing with the calibration factor
	}
//	total_energy=total_energy/caloenergy.size(); //normalize for number of towers, should apply a "percent of towers" but that can be a next step
	//std::cout<<"There are " <<caloenergy.size() <<"active towers in this run with energy " <<total_energy <<std::endl;
	return total_energy; 
}
int CaloTransverseEnergy::countTheTowers(plots* PLTS, double etamin, double etamax, double phimin, double phimax)
{
	int n_tow=0; 
	if(!PLTS->ohcal->filled_em_cal){
		if(PLTS->ohcal->counter){
			bool goodphi=false, goodeta=false;
			for(auto a:PLTS->ohcal->tower_edges){
				std::pair<double,double> etab=a.second.first;
				std::pair<double, double> phib=a.second.second;
				double emi=std::min(etab.first, etab.second);
				double emx=std::max(etab.first, etab.second);
				double pmi=std::min(phib.first, phib.second);
				double pmx=std::max(phib.first, phib.second);
				if(etamin >= emi && etamin <=emx) goodeta=true;
				else if (etamax >= emi && etamax <= emx) goodeta=true;
				if(phimin >= pmi && phimin <= pmx) goodphi=true;
				else if (phimax >= pmi && phimax <=pmx ) goodphi=true;
				if(goodeta && goodphi){
					int ieta=a.first.first, iphi=a.first.second;
					int chann=ieta*100+iphi;
					PLTS->ohcal->EMCAL_tower_proj->Fill(ieta, iphi);
					PLTS->ohcal->n_emcaltowers[chann]++;
				}
			}
		}
	}
	if(!PLTS->ihcal->filled_em_cal){
		if(PLTS->ihcal->counter){
			bool goodphi=false, goodeta=false;
			for(auto a:PLTS->ihcal->tower_edges){
				std::pair<double,double> etab=a.second.first;
				std::pair<double, double> phib=a.second.second;
				double emi=std::min(etab.first, etab.second);
				double emx=std::max(etab.first, etab.second);
				double pmi=std::min(phib.first, phib.second);
				double pmx=std::max(phib.first, phib.second);
				if(etamin >= emi && etamin <=emx) goodeta=true;
				else if (etamax >= emi && etamax <= emx) goodeta=true;
				if(phimin >= pmi && phimin <= pmx) goodphi=true;
				else if (phimax >= pmi && phimax <=pmx ) goodphi=true;
				if(goodeta && goodphi){
					int ieta=a.first.first, iphi=a.first.second;
					int chann=ieta*100+iphi;
					PLTS->ihcal->EMCAL_tower_proj->Fill(ieta, iphi);
					PLTS->ihcal->n_emcaltowers[chann]++;
					n_tow++;
				}
			}
		}
	}
	return n_tow;
}
float CaloTransverseEnergy::EMCaltoHCalRatio(float em, float hcal)
{
	//this is to project over to phenix data
	float ratio=em/hcal;
	return ratio; //again, first order correction, probably more to it than this
	
}
void CaloTransverseEnergy::FitToData(std::vector<float> data,int which_model) 
{
	//Apply fits from CGC, Glauber, others, meant to be generalizable
}
float CaloTransverseEnergy::Heuristic(std::vector<float> model_fit)
{
	float chindf=0;
	return chindf;
}
bool CaloTransverseEnergy::ValidateDistro()
{
	//Just make sure that the distribution looks correct
	double eavg=0, ee=0;
	std::vector<double> et;
	int bad=0;
	for(auto e:etphi){
		ee=0;
		for(auto v:e.second){
			 eavg+=v;
			 ee+=v;
		}
		et.push_back(ee);
	}  
	eavg=eavg/etphi.size();
	for(auto e:et) if(e > eavg*1.25 || e < eavg*0.75 ) bad++; 
	if(bad > etphi.size()*0.25) return false; //if over half of the phi values are outside of the range someting is off
	else return true;
}	
void CaloTransverseEnergy::ProduceOutput()
{
	
	TFile* outfile=NULL;
	if(truth) outfile=new TFile(Form("../data_output/MC/Transverse_Energy_truth_run_%d_segment_%d.root",run_number, DST_Segment), "RECREATE");
	else if(sim) outfile=new TFile(Form("../data_output/MC/Transverse_Energy_run_%d_segment_%d.root",run_number, DST_Segment), "RECREATE");
	else outfile=new TFile(Form("../data_output/Transverse_Energy_run_%d_segment_%d.root",run_number, DST_Segment), "RECREATE");
	for(auto p:jet_like_container)std::cout<<"The jet type thing has tower number " <<p.second.size() <<std::endl;
	findEEC();		
	for(auto p:zPLTS){
	plots* PLTS=p.second;
		PLTS->em->getAcceptance();
		PLTS->ihcal->getAcceptance();
		PLTS->ohcal->getAcceptance();

		PLTS->em->scaleThePlots(n_evt, &PLTS->event_data);
		PLTS->ihcal->scaleThePlots(n_evt, &PLTS->event_data);
		PLTS->ohcal->scaleThePlots(1, &PLTS->event_data);
		if(PLTS->em->z==0) std::cout<<"The plots have " <<PLTS->em->dET_eta->GetEntries() <<" entries with mean" <<PLTS->em->dET_eta->GetBinContent(PLTS->em->dET_eta->GetMaximumBin()) <<std::endl;
	
	for(auto data_point:PLTS->event_data){
		float total_val=0;
		total_val+=data_point["EMCAL"];
		total_val+=data_point["IHCAL"];
		total_val+=data_point["OHCAL"];
		PLTS->total->dET->Fill(total_val);
		PLTS->total->ET->Fill(total_val*2.2);
	}
/*	for(unsigned int i=0; i<PLTS->total->hists_1.size(); i++){
		 PLTS->total->hists_1.at(i)->Add(PLTS->em->hists_1.at(i));
		 PLTS->total->hists_1.at(i)->Add(PLTS->ihcal->hists_1.at(i));
		 PLTS->total->hists_1.at(i)->Add(PLTS->ohcal->hists_1.at(i));
	}

	for(unsigned int i=0; i<PLTS->total->hists_2.size(); i++){
		 PLTS->total->hists_2.at(i)->Add(PLTS->em->hists_2.at(i));
		 PLTS->total->hists_2.at(i)->Add(PLTS->ihcal->hists_2.at(i));
		 PLTS->total->hists_2.at(i)->Add(PLTS->ohcal->hists_2.at(i));
	}*/
//	PLTS->total->BuildTotal(PLTS);
	for(auto h:PLTS->em->hists_1) h->Write();
	for(auto h:PLTS->em->hists_2)h->Write();
	for(auto h:PLTS->ihcal->hists_1)h->Write();
	for(auto h:PLTS->ihcal->hists_2)h->Write();
	for(auto h:PLTS->ohcal->hists_1)h->Write();
	for(auto h:PLTS->ohcal->hists_2)h->Write();
	for(auto h:PLTS->total->hists_1)h->Write();
	for(auto h:PLTS->total->hists_2)h->Write();
	}
	for(auto p:szPLTS){
		plots* pts=p.second;
		pts->total->scaleThePlots(n_evt, &pts->event_data);
		for(auto h:pts->total->hists_1) h->Write();
		for(auto h:pts->total->hists_2) h->Write();
	}
	he2->Write();
	heec->Write();
	he3->Write();
	h3eec->Write();
	part_mass->Write();
	EMCALE->Write();
	IHCALE->Write();
	OHCALE->Write();
	/*TH1F* acceptance_eta_em_h=new TH1F("h", "H", 96, eEtaD->GetXaxis()->GetXmin(), eEtaD->GetXaxis()->GetXmax());
	
	int n_emtow=0, n_ihtow=0, n_ohtow=0; 
	seg_acceptance.em_eta_acceptance=std::vector<float> (96,0);
	seg_acceptance.ihcal_eta_acceptance=std::vector<float> (96,0);
	seg_acceptance.ohcal_eta_acceptance=std::vector<float> (96,0);
	for(int i=0; i<eeps->GetNbinsX(); i++){
		 for(int j=0; j<eeps->GetNbinsY(); j++){
				std::cout<<"The number of towers with events is " <<n_emtow <<std::endl; 
				if(eeps->GetBinContent(i,j) > 0){
					 n_emtow++;
					 seg_acceptance.em_eta_acceptance.at(i)+=1;
					 acceptance_eta_em_h->Fill(eEtaD->GetBinCenter(i));
				}
		}
	}
	for (int i=0; i<iheps->GetNbinsX(); i++) for(int j=0; j<iheps->GetNbinsY(); j++){
			 if(iheps->GetBinContent(i,j) > 0){
					 n_ihtow++;
					 seg_acceptance.ihcal_eta_acceptance.at(i)+=1;
			}
	}
	for (int i=0; i<oheps->GetNbinsX(); i++) for(int k=0; k<oheps->GetNbinsY(); k++){
			 if(oheps->GetBinContent(i,k) > 0)
			{
	 			n_ohtow++;
				seg_acceptance.ohcal_eta_acceptance.at(i)+=1;
			}
	}
	seg_acceptance.emcal_acceptance=256*96/float(n_emtow);
	seg_acceptance.ihcal_acceptance=64*24/float(n_ihtow);
	seg_acceptance.ohcal_acceptance=64*24/float(n_ohtow);
	//std::cout<<"acceptances " <<accept_e <<" , " <<accept_ih <<" , " <<accept_oh <<std::endl;
//	eEtaD->Scale(1/n_evt);
	acceptance_eta_em_h->Scale(1/256.);
	//eEtaD->Divide(acceptance_eta_em_h);
//	for(int i=0; i<eEtaD->GetNbinsX(); i++) eEtaD->SetBinContent(i+1, eEtaD->GetBinContent(i+1)*acceptance_eta_em.at(i)/256.);
	ePhiD->Scale(accept_e/n_evt);
	ohEtaD->Scale(accept_oh/n_evt);
	ohPhiD->Scale(accept_oh/n_evt);
	ihEtaD->Scale(accept_ih/n_evt);
	ihPhiD->Scale(accept_ih/n_evt);
//	for(int i=1; i<25; i++) EtaD->Fill(i, ihEtaD->GetBinContent(i)+ohEtaD->GetBinContent(i) + eEtaD->GetBinContent(4*i)+eEtaD->GetBinContent(4+i+1)+eEtaD->GetBinContent(4*i+2)+eEtaD->GetBinContent(4*i+3));
//	for(int i=1; i<33; i++) PhiD->Fill(i, ihPhiD->GetBinContent(i)+ohPhiD->GetBinContent(i) + ePhiD->GetBinContent(4*i)+ePhiD->GetBinContent(4*i+1)+ePhiD->GetBinContent(4*i+2)+ePhiD->GetBinContent(4*i+3));
	//EtaD->Scale(1/n_evt);
	//PhiD->Scale(1/n_evt);
	EtaD->Scale(1/(float)n_evt);
	//just make a ton of histos
	EMCALE->Scale(seg_acceptance.emcal_acceptance);
	EMCALE->Write();
	IHCALE->Scale(seg_acceptance.ihcal_acceptance);
	IHCALE->Write();
	OHCALE->Scale(seg_acceptance.ohcal_acceptance);
	OHCALE->Write();
	ETOTAL=(TH1F*)EMCALE->Clone();
	ETOTAL->SetName("total");
	ETOTAL->SetTitle("Total #frac{dE^{event}_{T}}{d#eta}; #sum_{tow} #frac{d E^{tow}_{T}}{d#eta}; N_{Evt}");
	ETOTAL->Add(IHCALE);
	ETOTAL->Add(OHCALE);
	PhiD->Write();
	EtaD->Write();
	phis->Write();
	etas->Write();
	ePhiD->Write();
	eEtaD->Write();
	ephis->Write();
	eetas->Write();
	ohPhiD->Write();
	ohEtaD->Write();
	hphis->Write();
	ihPhiD->Write();
	ihEtaD->Write();
	hetas->Write();
	eep->Write();
	eeps->Write();
	tep->Write();
	teps->Write();
	ihep->Write();
	iheps->Write();
	ohep->Write();
	oheps->Write();
	etabin_hc->Write();
	phibin_hc->Write();
	etabin_em->Write();
	phibin_em->Write();
//	EtaPhi->Write();*/
//
	
	//datatree->Write();
	headertree->Fill();
	headertree->Write();
	z_vertex->Write();
	outfile->Write();
	outfile->Close();
}
bool CaloTransverseEnergy::ApplyCuts(Event* e)
{
	//pass through a set of cuts
	return true;
}
void CaloTransverseEnergy::GetNodes(PHCompositeNode *topNode)
{
	_topNode=topNode;
	_IHCALNode=findNode::getClass<PHCompositeNode>(topNode, iHCALnode);
	_OHCALNode=findNode::getClass<PHCompositeNode>(topNode, oHCALnode); 
	try{
		_EMCALNode=findNode::getClass<PHCompositeNode>(topNode, EMCALnode);
		led=false;
	}
	catch(std::exception& e){led=true;}
}
int CaloTransverseEnergy::Init(PHCompositeNode *topNode)
{
//	std::cout<<"initializing analysis on run " <<run_number <<" segement " <<DST_Segment <<std::endl;
	if(!isPRDF) GetNodes(topNode); //load the composite nodes into the global variables if using DST approach
	else{
		for(int i=0; i<24; i++) n_evt_eta.push_back(0);//std::string outname=("Transverse_Energy_%i_segment_%i.root", run_number, DST_Segment);
	}
	EtaPhi=new TNtuple("etaphi", "Eta, Phi Transverse Energy, EmCal ET, iHCal ET, oHCal ET", "eta:phi:et:eet:het");
	EMEtaPhi=new TNtuple("emetaphi", "EMCal #eta, #varphi E_{T}", "eta:phi:et");
	HEtaPhi=new TNtuple("hetaphi", "HCal #eta #varphi E_{T}", "eta, phi, et");
	datatree=new TTree("CALOET", "CALOET");
	datatree->Branch("total_energy", &energy);
	datatree->Branch("hcal_energy", &hcalenergy);
	datatree->Branch("emcal_energy", &emcalenergy);
	datatree->Branch("energy_transverse", &energy_transverse);
	datatree->Branch("angular_data_em", &evt_data);
	headertree=new TTree("HEADER", "RUN DATA");
	headertree->Branch("run number", &run_number);
	headertree->Branch("segment number", &DST_Segment);
//	headertree->Branch("acceptances", &seg_acceptance);
	headertree->Branch("number of events", &n_evt);
//	headertree->Branch("cuts", &kinematiccuts);
	//datatree->Branch("EtaPhi", EtaPhi);
	//datatree->Branch("EM_Eta_Phi", EMEtaPhi);
	//datatree->Branch("HC_Eta_Phi", HEtaPhi);
	for(int i =0; i<8; i++)
	{
		int ih=7001+i;
		int oh=8001+i;
		packets.push_back(ih);
		packets.push_back(oh);
	}
	//for(int i=0; i<128; i++){
	//	int em=6001+i;
	//	packets.push_back(em);
	//}
	return 0;	 
}
void CaloTransverseEnergy::findEEC()
{
	float ET=0;	
	for(auto q:jet_like_container){
		ET=0;
		for(auto m:q.second){
				ET+=m.second;	
			 }
		float eec=0, e3c=0;
		for(auto m:q.second){
			for(auto p:q.second){
				float DE=((float)p.first.first - (float)m.first.first), DP=(float)std::min((float)std::min((float)abs((float)p.first.second-(float)m.first.second), (float)abs((float)p.first.second-(float)m.first.second-6.283)), (float)abs((float)p.first.second-(float)m.first.second+6.283));
				float R=sqrt(pow(DE, 2.0) + pow(DP, 2.0));
	//			std::cout<<"R val is : " <<R <<"with " <<DE <<" and " <<DP <<std::endl;
				float EEC=p.second*m.second/(float)pow(ET, 2.0);
				eec+=EEC;
				he2->Fill(R, EEC);
			for(auto n:q.second){
				if(m.first==n.first || m.first == p.first || n.first == p.first) continue;
				float DE2=((float)n.first.first - (float)m.first.first), DP2=((float)n.first.second-(float)m.first.second);
				DP2=std::min(std::min(abs(DP2), (float)abs(DP2-6.283)), (float)abs(DP2-6.283));
				float DE3=((float)p.first.first - (float)n.first.first), DP3=((float)p.first.second-(float)n.first.second);
				DP3=std::min(std::min(abs(DP3), (float)abs(DP3-6.283)), (float)abs(DP3-6.283));
				float R2=sqrt(pow(DE2,2.0)+pow(DP2,2.0)), R3=sqrt(pow(DE3, 2.0)+pow(DP3, 2.0));
				float Rval[]={R, R2, R3};
				float RM=*std::max_element(Rval, Rval+3);
				if(RM>0.8) continue;
				//std::cout<<"Good tripplet with " <<p <<m <<n <<std::endl;
				//	std::cout<<"Have the EEC as " <<EEC <<std::endl; 
				//	std::cout<<"Filled that into the function" <<std::endl;
				
				float E3C=p.second*m.second*n.second/(float)(pow(ET, 3.0));
				//std::cout<<"The 3 pt gives " <<E3C<<std::endl;
				he3->Fill(RM, E3C);
			//	std::cout<<"Filled 3 pt in for " <<RM <<std::endl;
				e3c+=E3C;
			//	std::cout<<"issue is not here, so where is it?" <<std::endl;
			}
		}
	}
		
		heec->Fill(eec);
		h3eec->Fill(e3c);
		std::cout<<"one of the jet energies " <<eec <<std::endl;
	}
}		
