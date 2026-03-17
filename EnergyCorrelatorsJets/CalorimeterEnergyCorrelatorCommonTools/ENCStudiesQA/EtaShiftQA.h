#ifndef __ETASHIFTQA_H__
#define __ETASHIFTQA_H__

#include <vandyclasses/Tower.h>
#include <vandyclasses/EventInfo.h>
#include <vandyclasses/Jet.h>
#include <caloenctools/BuildMetaTowers.h>

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>

#include <math.h>
#include <array>
#include <string>
#include <vector>
#include <format>

class PerCaloQAPlots
{
	public:
		PerCaloQAPlots(std::string calorimeter)
		{
			this->calo = calorimeter;
			InitializePlots();
		}
		~PerCaloQAPlots(){};
		TH1F* Deltaetabin;
		TH1F* deltaEt;
		TH1F* shifteta;
		TH2F* etaDeltaetabin;
		TH2F* zVTXdeltaeta;
		std::string calo {"Meta"};
	private:
		void InitializePlots()
		{
			//this is to initialize the plots for each calorimeter
			Deltaetabin 	= new TH1F(
						std::format("h_{}_Detabin", calo).c_str(), 
						std::format("{} bin shift; #Delta #eta_{bin}; N_{bins}", calo).c_str(),
						100, -2, 2
					);
			deltaEt		= new TH1F(
						std::format("h_{}_dEt", calo).c_str(),
						std::format("Change in calculated E_{T} for {} ; #Sum_{{}} #Delta E_{T}; N_{events}", calo, calo).c_str(),
						2000, -1000.5, 999.5
					);
			shifteta	= new TH1F(
						std::format("h_{}_shift", calo).c_str(),
						std::format("New #eta of center {}; #eta; N_{hits}", calo).c_str(),
						100, -4, 4
					);


		}
};
class EtaShiftQA
{
	public:
		EtaShiftQA();
		~EtaShiftQA(){};

	private:
		TH1F* etaShiftEMCAL;
		TH1F* etaShiftOHCAL;
		TH1F* etaShiftIHCAL;
		TH1F* etaShiftEMCALRMeta;
		TH1F* etaShiftOHCALRMeta;
		TH1F* etaShiftIHCALRMeta;
		TH1F* calculatedJetpt;
		TH1F* calculatedShiftedJetpt;
		TH1F* zVTX;
		TH2F* etaZdeta;
		TH2F* etabineta;


};
#endif

