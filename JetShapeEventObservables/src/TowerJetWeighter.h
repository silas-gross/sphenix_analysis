#ifndef __TOWERJETWEIGHTER_H__
#define __TOWERJETWEIGHTER_H__

#include <TH1.h>
#include <TTree.h>

#include <array>
#include <vector>
#include <map>
#include <string>
#include <thread>
#include <cmath>

struct WeightedTower
{
	WeightedTower(int toweridi, float ri=-1, float phii=-100, float etai=-100, float ETi=-1, float pt_mini {-1}, std::string CaloLabeli=""):
	       	towerid(toweridi),
		r(ri),
	       	phi(phii),
	       	eta(etai),
	       	ET(ETi),
	       	pt_min(p_mini),
		CaloLabel(CaloLabel)
	{}
	int towerid {};
	float r {-1};
	float phi {-100}; 
	float eta {-100};
	float ET {-1};
	float pt_min {-1};
	std::string CaloLabel {""};
	std::map<float, float> weight {};
	std::map<float, float> ET_cut {};
};
class ShapeTrim
{
	public:
		ShapeTrim(float pt_min_i = 0.1, std::vector<float> R_vals_i ={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}){
			this->pt_min = pt_min;
			this->R_vals = R_vals_i;
		};
		~ShapeTrim(){};
		void FillWeightedTowers(std::map<int, std::array<float, 4>> inputTowers, std::string CaloLabel="")
		{
			PrepareNewEvent();
			
			for(auto i:inputTowers) 
			{
				WeightedTower* tw = new WeightedTower(i.first, i.second[0], i.second[1], i.second[2], i.second[3], this->Et_min, CaloLabel);
				this->Towers->insert_or_assign(i.first, tw);
			}
			std::vector<std::thread> TowerCollectorThreads;
			std::map< int, std::map<float, std::vector<int>>> tempHoldingSol {};
			for(auto i:*Towers)
			{
				std::map<float, std::vector<int>> holdingSolTw {};
				for(auto R:this->R_vals){
					holdingSolTw[R]=std::vector<int> {};
				}
				tempHoldingSol[i.first]=holdingSolTW;
				TowerCollectorThreads.push_back(std::thread(&ShapeTrim::CollectRelevantTowers, this, i.second, inputTowers, &(tempHoldingSol[i.first])));
			}i
			for(int i = 0; i<(int)TowerCollectorThreads.size(); i++) TowerCollectorThreads.at(i).join();
			TowerCollectorThreads.clear();
			for(auto tw:tempHoldingSol)
			{
				TowerCones[tw.first]=tw.second;
			}
			CalculateAllETiR();
			return;
		}
		std::map<int, WeightedTower*>* getTowers(){return this->Towers};
	private:
		void PrepareNewEvent()
		{
			Towers->clear(); 
			TowerCones.clear();
			Efficiency_R.clear();
			return;
		};
		void CollectRelevantTowers(std::array<float, 4> basetw, std::map<int, std::array<float, 4>> tows, std::pair<int, std::map<float, std::vector<int>>>* storageval )
		{
			for(auto R: storageval->second)
			{
				for(auto tw: tows) 
				{
					float Rij = calcR(basetw[1], basetw[2], tw.second[1], tw.second[2]);
					if(Rij <= R.first) R.second.push_back(tw.first);
					else continue;
				}
			}
			return;
		}; //just gets a relevant circle around each tower for each value of R to calculate the cone 
		void CaluclateAllETiR()
		{
			std::map<int, float> ref_copy {};
			for(auto i:*Towers)
			{
				ref_copy[i.first]=i.second->ET;
			}
			std::vector<std::thread> CalcTowerEnergyThread {};
			std::map<int, std::map<float, float>> tower_energygroups {}; //eventually what is fed into the ref 
			for(auto TwC:TowerCones) tower_energygroups[TwC.first]=std::map<float, std::vector<float>> {};
			for(auto TwC:TowerCones)
			{
				std::map<float, float> teg {};
				tower_energygroups[TwC.first].push_back(teg);
				CalcTowerEnergyThread.push_back(&ShapeTrim::DoTowerConeCalc, this, TwC.second, ref_copy, &(tower_energygroups[TwC.first].back()));
			}
			for(int i=0; i<(int)CalcTowerEnergyThread.size(); i++) CalcTowerEnergyThread.at(i).join();
			for(auto ETTw:tower_energygroups)
			{
				Towers->at(ETTw.first)->ET_cut =  ETTw.second;
			}
			return;
		}
		void DoTowerConeCalc(std::map<float, std::vector<int>> inputTowerList, std::map<int, float> ref, std::map<float, std::vector<float>>* output_energy )
		{
			for(auto inp:inputTowerList)
			{
				std::vector<float> otE;
				CollectTowerEnergies(inp.second, ref, &otE);
				float E = CalculateETiROne(otE);
				output_energy->at(inp.first) = E;
			}
			return;
		}
		void CollectTowerEnergies(std::vector<int> ref, std::map<int, float> lookup, std::vector<float>* outputTowerEnergy)
		{
			for(auto iv:ref) 
			{
				try{
					outputTowerEnergy->push_back(lookup[iv]);
				}
				catch(std::exception& e) {continue;}
			}
			return;
		}
		float CaluclateETiROne(std::vector<float> relevantTowerE)
		{
			float ET_iR=0;
			for(auto i:relevantTowerE) ET_iR+=i;
			return ET_iR;
		} //get the value of E_{T, iR} needs to just take in a processed value of the Energies

		float calcR(float phi1, float eta1, float phi2, float eta2)
		{
			float dphi = calcdPHI(phi1, phi2);
			float deta = eta1 - eta2;
			float R = std::sqrt(std::pow(dphi, 2) + std::pow(deta, 2));
			return R;
		}
		float calcdPHI(float phi1, float phi2)
		{
			float dphi = phi1-phi2;
			if(dphi > M_PI) dphi = 2*M_PI - dphi;
			else if( dphi < -M_PI ) dphi = -2*M_PI - dphi;
			return dphi;
		}
		float cutTowerETR(std::vector<float> relevant_E, float pt_min=1.);//this function runs the skiming to send ET to ET_iR 
		std::map<int, WeightedTower*>* Towers {nullptr};
		std::map<int, std::map<float, std::vector<int>>> TowerCones {};
		float Et_min {-1}; //the minimum cut --default for EMCAL = 100 MeV, default for HCAL = 650 MeV
		std::vector<float> R_vals {}; //which values of R are we using (default 0.1-0.8 in 0.1 bins) 
		std::map<float, float> Efficiency_R {}; // this is just the percentage of incoming towers that have non-zero weight at each R

};
#endif
