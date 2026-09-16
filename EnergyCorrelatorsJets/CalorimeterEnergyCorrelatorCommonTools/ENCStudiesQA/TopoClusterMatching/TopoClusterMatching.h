// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TOPOCLUSTERMATCHING_H
#define TOPOCLUSTERMATCHING_H

#include <fun4all/SubsysReco.h>

#include <string>

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
		eta = cl->get_eta();
		Et = E / std::sinh(eta);
		pt = Et;
	}
	~CandidateObj(){};
	float E {-999.};
	float Et {-999};
	float pt {-999.};
	float phi {-999.};
	float eta {-999.};
	bool isMatch {false};

}

class TopoClusterMatching : public SubsysReco
{
 public:

  TopoClusterMatching(float mpT=0.0, [[maybe_unused]] const std::string &name = "TopoClusterMatching");

  ~TopoClusterMatching() override;

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
  int End([[maybe_unused]] PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  //void Print(const std::string &what = "ALL") const override;

 private:
  	int n_evt	{0};
	int n_good	{0};
	float minpt	{0.};
	float Q2	{-999.};
	float erat	{0.5};
	float dRmax	{0.2};

	std::pair<int, int> jet_bin_index {0,0};
	void getBins	(float);
	bool runDijetCut(PHCompositeNode*);
	bool isAMatch	(CandidateObj*, CandidateObj*);
	float getR 	(CandidateObj*, CandidateObj*);	

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
};

#endif // TOPOCLUSTERMATCHING_H
