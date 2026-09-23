// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TOPOCLUSTERMATCHING_H
#define TOPOCLUSTERMATCHING_H

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

//vertex stuff
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

//G4 objects

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>


#include <phhepmc/PHHepMCGenEvent.h>  
#include <phhepmc/PHHepMCGenEventMap.h>
#include <HepMC/GenEvent.h>

//jetbase objects 

#include <jetbase/JetContainer.h>
#include <jetbase/JetContainerv1.h>
#include <jetbase/Jetv1.h>

//calo tower and cluster stuff 
//Calo towers 
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfov2.h>
#include <calobase/TowerInfov1.h>
#include <calobase/TowerInfo.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>

#include "DijetEventCuts.h"

#include <string>
#include <math.h>
#include <vector>
#include <utility>
#include <format>
#include <array>

//root
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

class PHCompositeNode;

class CandidateObj
{
 public:
	CandidateObj(PHG4Particle* truth) 
	{
		//puts a truth particle in this stripped down form
		E = truth->get_e();
		float px = truth->get_px();
		float py = truth->get_py();
		float pz = truth->get_pz();
		pt = std::sqrt(std::pow(px, 2) + std::pow(py, 2));
		phi = std::atan2(py, px);
		eta = std::atanh(pz / E);
		Et = E / std::sinh(eta);
	}
	CandidateObj(RawCluster* cl)
	{
		//puts a cluster in this stripped down form 
		E = cl->get_energy();
		phi = cl->get_phi();
		float r = cl->get_r();
		float z = cl->get_z();
		eta = std::asinh(z/r);
		Et = E / std::cosh(eta);
		pt = Et;
	}
	~CandidateObj(){};
	float E {-999.};
	float Et {-999};
	float pt {-999.};
	float phi {-999.};
	float eta {-999.};
	bool isMatch {false};

};

class TopoClusterMatching : public SubsysReco
{
 public:

  TopoClusterMatching(const float mpT=0.0, const std::string &name = "TopoClusterMatching");

  ~TopoClusterMatching(){};

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
//  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  //int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
 // int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
 // int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
//  int End([[maybe_unused]] PHCompositeNode *topNode) override;

  /// Reset
//  int Reset(PHCompositeNode * /*topNode*/) override;

  //void Print(const std::string &what = "ALL") const override;

 private:
  	int n_evt	{0};
	int n_good	{0};
	float minpt	{0.};
	float Q2	{-999.};
	float eratmin	{0.5};
	float dRmax	{0.2};
	
	std::string jet_node_name {"Truth_AntikT_r04"};
	std::pair<int, int> jet_bin_index {0,0};
	void getBins	();
	bool runDijetCut(PHCompositeNode*);
	bool isAMatch	(CandidateObj*, CandidateObj*);
	float getR 	(CandidateObj*, CandidateObj*);	
	void getEEC(std::vector<CandidateObj*>, bool );

	bool KinCuts(PHG4Particle* );
	bool KinCuts(RawCluster* );	
	void MatchAllTruthToClusters	(PHCompositeNode*);
	std::pair<int, int> findBucket	(float, float);
	
		
	std::vector<float> clusterPTBin {};
	std::vector<float> lead_bucket {};
	std::vector<std::vector<float>> sub_bucket {};	

	DijetEventCuts* event_cut {nullptr};
	
	TH1F* h_ClusterPairET {nullptr};
	TH1F* h_TruthPairET {nullptr};
	
	std::array< std::array<TH1F*, 10>, 10> h_ClusterPairET_Div{};  
	std::array< std::array<TH1F*, 10>, 10> h_TruthPairET_Div{};  
	
	TH1F* h_ClusterPairETAll {nullptr};
	TH1F* h_TruthPairETAll {nullptr};
	
	std::array< std::array<TH1F*, 10>, 10> h_ClusterPairETAll_Div{};  
	std::array< std::array<TH1F*, 10>, 10> h_TruthPairETAll_Div{};  
	
	TH1F* h_MatchedTruth {nullptr};
	TH1F* h_RealCluster {nullptr};
};

#endif // TOPOCLUSTERMATCHING_H
