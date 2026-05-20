// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SHOWERTOWERMATCHING_H
#define SHOWERTOWERMATCHING_H
//fun4all basic stuff
#include <fun4all/SubsysReco.h>

#include "BuildMetaTowers.h"
//root 
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>

//c++
#include <string>
#include <vector>
#include <math.h>


class PHCompositeNode;
class tower
{
	public:
		tower(
			float etaC = -999.,
			float phiC = -999. , 
			float e=-999., 
			int tn=-999, 
			bool isT=false
		){
			etaCenter = etaC;
			phiCenter = phiC;
			E 	= e;
			ET	= E / std::cosh(etaCenter);
			tower_N = tN;
			isTruth = isT;
		};
		tower(
			std::array<float, 3> etaC { -999., -999., -999.},
			std::array<float, 3> phiC { -999., -999., -999.},
			float e=-999., 
			int tn=-999, 
			bool isT=false
		){
			etalow		= etaC[0];
			etaCenter 	= etaC[1];
			etahigh		= etaC[2];

			philow		= phiC[0];
			phiCenter 	= phiC[1];
			phihigh		= phiC[2];

			E 	= e;
			ET 	= E / std::cosh(etaCenter);
			tower_N = tN;
			isTruth = isT;
		};
		~tower(){};
		int tower_N {-999};
		float etaCenter {-999.};
		float etalow {-999.};
		float etahigh {-999.};
		float phiCenter {-999.};
		float philow {-999.};
		float phihigh {-999.};
		float E {0};
		float ET {0};
		bool isTruth {false};
};
class Shower
{
	public:
		Shower(){};
		~Shower(){};
		void AddTower( tower tw)
		{
			StruckTowers.push_back(tw);
			if( tw.etalow != -999 && tw.etahigh != -999) 
			{
				if(tw.etalow < etaDown) etaDown = tw.etalow;
				if(tw.etahigh > etaUp) etaUp = tw.etahigh;
			}
			else 
			{
				if(tw.etaCenter < etaDown) etaDown = tw.etaCenter;
				if(tw.etaCenter > etaUp) etaUp = tw.etaCenter
			}
			if( tw.philow != -999 && tw.phihigh != -999) 
			{
				if(tw.philow < phiDown) phiDown = tw.philow;
				if(tw.phihigh > phiUp) phiUp = tw.phihigh;
			}
			else 
			{
				if(tw.phiCenter < phiDown) phiDown = tw.phiCenter;
				if(tw.phiCenter > phiUp) phiUp = tw.phiCenter
			}
			return;
		}
		bool IsInHits(tower tw)
		{
			bool geomMatch {false};
			if(tw.etaCenter > etaDown && tw.etaCenter < etaUp )
				if(tw.phiCenter > phiDown && tw.phiCenter < phiUp) 
					geomMatch=true;
			return geomMatch;
		}

		std::vector<Tower> getStruck() {	
			return StruckTowers;
		};
		std::array<float, 2> get_etaBounds() { 
			return std::array<float, 2> {etaDown, etaUp};
		};
		std::array<float, 2> get_phiBounds() { 
			return std::array<float, 2> {etaDown, etaUp};
		};
	private:
		std::vector<Tower> StruckTowers{};
		float etaDown	{-999.};
		float etaUp 	{-999.};
		float phiDown	{-999.};
		float phiUp 	{-999.};

};
class ShowerTowerMatching : public SubsysReco
{
	public:

		ShowerTowerMatching(const std::string &name = "ShowerTowerMatching");

		~ShowerTowerMatching() override;

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

		/** Called for each event.
		This is where you do the real work.
		*/
		int process_event(PHCompositeNode *topNode) override;

		/// Clean up internals after each event.

		/// Called at the end of each run.

		/// Called at the end of all processing.
		int End(PHCompositeNode *topNode) override;


	private:
		void buildTowerBins(int n_bins = 10);
		void matchTheTower();
		void KinCuts(PHG4Particle*);
		void getParticleShower(
			PHG4Particle*,
			PHG4Shower*,
			Shower*,
			PHCompositeNode* 
			);	
		void buildTruthTowers(
				std::map<PHG4Particle*, Shower*>, 
				std::vector<PHG4Particle*> );

		std::vector<tower> shower_having_truth_towers {};
		std::vector<tower> unshowerd_truth_towers {};

		std::vector<float>* tower_bins = new std::vector<float> {};

		TH1F* h_tow_fake {nullptr};
		TH1F* h_tow_miss {nullptr};
		TH1F* h_cls_fake {nullptr};
		TH1F* h_cls_miss {nullptr};
		
		TH1F* h_tow_fake_tr {nullptr};
		TH1F* h_tow_miss_tr {nullptr};
		TH1F* h_cls_fake_tr {nullptr};
		TH1F* h_cls_miss_tr {nullptr};
};

#endif // SHOWERTOWERMATCHING_H
