// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EMCALCHANNELTIMING_H
#define EMCALCHANNELTIMING_H

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllReturnCodes.h>

//phool 
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <ffaobjects/EventHeader.h>

#include <calobase/TowerInfov4.h>
#include <calobase/TowerInfoContainerv4.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TDirectory.h> 

#include <string>
#include <map>
#include <vector>
#include <math.h>
#include <utility>
#include <format>

class PHCompositeNode;

class tower
{
 public:
	tower(
		float ph = -999., 
		float et = -999., 
		float e  = -999., 
		float T  = -999., 
		int ei   = -999.
	)
	{
		phi = ph;
		eta = et;
		E   = e;
		t   = T;
		Ei  = ei;
	}
	~tower(){};
	float phi {-999};
	float eta {-999};
	float E	  {-999};
	float t   {-999};
	int Ei    {-999};
};

class EMCALChannelTiming : public SubsysReco
{
 public:

  EMCALChannelTiming(int seegm = 0, int runm=0,  const std::string &name = "EMCALChannelTiming");

  ~EMCALChannelTiming(){};

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init( [[maybe_unused]] PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
 // int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
//  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
//  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End( [[maybe_unused]] PHCompositeNode* topNode) override;

  /// Reset
//  int Reset(PHCompositeNode * /*topNode*/) override;

// void Print(const std::string &what = "ALL") const override;

 private:
  float SubdivideDetector( std::vector<tower*>*, std::vector<tower*>*, std::vector<tower*>*, PHCompositeNode*);
  void AnaHelper( std::vector<tower*>, std::vector<TH1F*>*, std::vector<TH2F*>*, float);
  void AnaHelper( tower*, std::vector<TH1F*>*, std::vector<TH2F*>*, float);
  int getIndex(float, float);
  void AddBins(std::vector<TH1F*>*, std::vector<TH2F*>*,
		  std::vector<TH1F*>*, std::vector<TH2F*>*,
		  std::vector<TH1F*>*, std::vector<TH2F*>*,
		  std::string ntower="");
  std::string emcal_tower {"TOWERINFO_CALIB_CEMC"};
  std::string emcal_geom {"TOWERGEOM_CEMC"};
  std::string output_file_name {"A.root"};
  enum a1DOUTPUTHISTS
  {
	 DELTAT, 
	 E,
	 EBART,
  };
  enum a2DOUTPUTHISTS
  {
	  EtoT,
	  EtaPhi,
  };
  int seg	{0};
  int run  	{0};
  int nevt 	{0};
  std::vector<Double_t> energy_bins {}; //log bins for energy
  
  //histograms across all towers
  //1D
  std::vector<TH1F*>* AllTowers1D=new std::vector<TH1F*> {};
  std::vector<TH1F*>* highTowers1D=new std::vector<TH1F*> {};
  std::vector<TH1F*>* lowTowers1D=new std::vector<TH1F*> {};

  //2D
  std::vector<TH2F*>* AllTowers2D=new std::vector<TH2F*> {};
  std::vector<TH2F*>* highTowers2D=new std::vector<TH2F*> {};
  std::vector<TH2F*>* lowTowers2D=new std::vector<TH2F*> {};

  //histograms for individual towers 
  //1D
/*  std::vector< std::vector<TH1F*>*>* AllTowers1D_t=new std::vector< std::vector<TH1F*>*> {};
  std::vector< std::vector<TH1F*>*>* highTowers1D_t=new std::vector< std::vector<TH1F*>*> {};
  std::vector< std::vector<TH1F*>*>* lowTowers1D_t=new std::vector< std::vector<TH1F*>*> {};

  //2D
  std::vector< std::vector<TH2F*>*>* AllTowers2D_t=new std::vector< std::vector<TH2F*>*> {};
  std::vector< std::vector<TH2F*>*>* highTowers2D_t=new std::vector< std::vector<TH2F*>*> {};
  std::vector< std::vector<TH2F*>*>* lowTowers2D_t=new std::vector< std::vector<TH2F*>*> {};
*/
  std::vector<tower*>* allE = new std::vector<tower*> {}; 
  TTree* towerTree {nullptr};
};

#endif // EMCALCHANNELTIMING_H
