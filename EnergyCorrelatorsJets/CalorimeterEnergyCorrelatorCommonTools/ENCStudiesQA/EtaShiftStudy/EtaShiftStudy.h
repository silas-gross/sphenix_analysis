// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ETASHIFTSTUDY_H
#define ETASHIFTSTUDY_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

//fun4all basic stuff
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllReturnCodes.h>

//phool
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

//vertex stuff
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

//root
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

//personal
#include "../../BuildMetaTowers.h" //need to fix this later 

//c++
#include <math.h>
#include <array>
#include <string>
#include <vector>
#include <format>


class PHCompositeNode;

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
		TH1F* avgeta;
		TH1F* Et;
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
			avgeta	 	= new TH1F(
						std::format("h_{}_eta", calo).c_str(), 
						std::format("{} average #eta; <#eta >; N_{bins}", calo).c_str(),
						100, -2, 2
			deltaEt		= new TH1F(
						std::format("h_{}_dEt", calo).c_str(),
						std::format("Change in calculated E_{T} for {} ; #Sum_{{}} #Delta E_{T}; N_{events}", calo, calo).c_str(),
						2000, -100.5, 99.5
					);
			Et		= new TH1F(
						std::format("h_{}_Et", calo).c_str(),
						std::format("Transverse Energy for {}; #Sum{{}} E_{T}; N_{evts}", calo, calo).c_str(),
						1000, 0, 200
					);
			shifteta	= new TH1F(
						std::format("h_{}_shift", calo).c_str(),
						std::format("New #eta of center {}; #eta; N_{hits}", calo).c_str(),
						100, -4, 4
					);
			etaDeltaEtabin	= new TH2F(
						std::format("h_{}_eta_deta", calo).c_str(),
						std::format("{} bin shift; #eta_{phyiscal}; #Delta #eta_{bin}; N_{evts}", calo).c_str(),
					       	24, -1.1, 1.1,
						100, -2, 2,
					);
			zVTXdeltaeta	= new TH2F(
						std::format("h_{}_z_deta", calo).c_str(),
						std::format("Event average #Delta #eta in {}; z_{vtx}; < #eta >; N_{evts}", calo).c_str(),
						120, -60.5, 59.5,
						100, -4, 4
					);
			

		}
};

class EtaShiftStudy : public SubsysReco
{
	public:

		EtaShiftStudy(const std::string &name = "EtaShiftStudy");

		~EtaShiftStudy() override;

		int Init(PHCompositeNode *topNode) override;


		int process_event(PHCompositeNode *topNode) override;


		/// Called at the end of each run.
		int EndRun(const int runnumber) override;

		/// Called at the end of all processing.
		int End(PHCompositeNode *topNode) override;

		/// Reset
		int Reset(PHCompositeNode * /*topNode*/) override;

		void Print(const std::string &what = "ALL") const override;
		double caluclatedAvgEta(
				std::array<BuildMetaTowers::TowerArrayEntry*, 1536>* calovals
		)
		{
			double avgeta=0.;
			int n=0;
			for(int i = 0; i<(int)calovals->size(); i++)
			{
				if(calovals->at(i)->Energy <= 0.) continue; //don't bother with empty bins, we want to not include dead areas in a systemaic analysis
				avgeta+=calovals->at(i)->eta;
				n++;
			}
			avgeta = avgeta / (double)n;
			return avgeta;
		};
		float CalculateJetPt(
				std::vector<std::array<float, 2>>
				);
		void AnalyzeEvent(PHCompositeNode*);
		void grabTowerArray(
				BuildMetaTowers*, 
				std::array<std::array<BuildMetaTowers::TowerArrayEntry*, 1536>* 4>*
				);
		void compareTowerValue( 
				std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*, 
				std::array<BuildMetaTowers::TowerArrayEntry*, 1536>*,
				float,
				int,
				BuildMetaTowers*
				);
		void grabJetConsitents(
			PHCompositeNode*,
			float			
			);
	private:
		//over all z 
		PerCaloQAPlots* EMCALQA;
		PerCaloQAPlots* IHCALQA;
		PerCaloQAPlots* OHCALQA;
		//using EMCAL radius, IHCAL radius and OHCAL radius respectively
		PerCaloQAPlots* MetaEQA;
		PerCaloQAPlots* MetaIQA;
		PerCaloQAPlots* MetaOQA;

		TH1F* calculatedJetpt;
		TH1F* calculatedShiftedJetpt;
		
		TH1F* hzVTX;
		//meta tower builder and shifter

		//restricted z ranges
		std::array<PerCaloQAPlots*, 6>* EMCAL_Z_QA;
		std::array<PerCaloQAPlots*, 6>* IHCAL_Z_QA;
		std::array<PerCaloQAPlots*, 6>* OHCAL_Z_QA;
		std::array<PerCaloQAPlots*, 6>* MetaE_Z_QA;
		std::array<PerCaloQAPlots*, 6>* MetaI_Z_QA;
		std::array<PerCaloQAPlots*, 6>* MetaO_Z_QA;

		BuildMetaTowers* emMetaTowerBuilder 	{ nullptr };
		BuildMetaTowers* ihMetaTowerBuilder	{ nullptr };
		BuildMetaTowers* ohMetaTowerBuilder	{ nullptr };


};

#endif // ETASHIFTQA_H
