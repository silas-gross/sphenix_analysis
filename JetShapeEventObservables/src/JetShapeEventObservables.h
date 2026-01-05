// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef JETSHAPEEVENTOBSERVABLES_H
#define JETSHAPEEVENTOBSERVABLES_H

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include <string>
#include <vector>
#include <array>

#include "TowerJetWeighter.h"
#include "SoftDropComp.h"

class PHCompositeNode;
class ShapeTrim;
struct WeightedTower;

class JetShapeEventObservables : public SubsysReco
{
	public:

		JetShapeEventObservables(const std::string &name = "JetShapeEventObservables");

		~JetShapeEventObservables() override;

		/** Called during initialization.
		Typically this is where you can book histograms, and e.g.
		register them to Fun4AllServer (so they can be output to file
		using Fun4AllServer::dumpHistos() method).
		*/
		int Init(PHCompositeNode *topNode) override;

		/** Called for first event when run number is known.
		Typically this is where you may want to fetch data from
		database, because you know the run number. A place
		to book histograms which have to know the run number.
		*/
		int InitRun(PHCompositeNode *topNode) overridei {return Fun4AllReturnCodes::EVENTOK;}

		/** Called for each event.
		This is where you do the real work.
		*/
		int process_event(PHCompositeNode *topNode) override;

		/// Clean up internals after each event.
		int ResetEvent(PHCompositeNode *topNode) override {return Fun4AllReturnCodes::EVENTOK;}

		/// Called at the end of each run.
		int EndRun(const int runnumber) override {return Fun4AllReturnCodes::EVENTOK;}

		/// Called at the end of all processing.
		int End(PHCompositeNode *topNode) override {return Fun4AllReturnCodes::EVENTOK;}

		/// Reset
		int Reset(PHCompositeNode * /*topNode*/) override {return Fun4AllReturnCodes::EVENTOK;}

		void Print(const std::string &what = "ALL") const override;

	private:
		void filterTowerInput(PHCompositeNode* topNode);
		void getTowerss(PHCompositeNode* topNode);
		void getTruth(PHCompositeNode* topNode);
		void getClusters(PHCompositeNode* topNode);
		void MakeEMCALRetowerMap(RawTowerGeomContainer_Cylinderv1*, TowerInfoContainer*, RawTowerGeomContainer_Cylinderv1*, TowerInfoContainer*);
		
		ShapeTrim* ShapeTrimmer {nullptr};
		std::vector<TH1F*> efficiencies {};
		std::map<int, WeightedTowers*> trimmedTowers;	
		std::map<std::string, std::map<int, std::array<float, 4>>> filterInputTowers;
		//structure "Type Tag":{particle/cluster/tower id : {eta, phi, r, E}}
		TH1F *h_n_ohcal_clusters, *h_n_emcal_clusters;
		std::string emcal_energy_towers {"TOWERINFO_CALIB_CEMC"}, ihcal_energy_towers {"TOWERINFO_CALIB_HCALIN"}, ohcal_energy_towers {"TOWERINFO_CALIB_HCALOUT"};
		std::map<int, std::pair<float, float>> emcal_lookup_table {};
};

#endif // JETSHAPEEVENTOBSERVABLES_H
