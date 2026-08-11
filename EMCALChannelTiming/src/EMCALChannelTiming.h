// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EMCALCHANNELTIMING_H
#define EMCALCHANNELTIMING_H

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <string>

class PHCompositeNode;

class tower
{
	tower(
		float ph, 
		float et, 
		float e, 
		float T, 
		int16_t ei)
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
	float e	  {-999};
	float t   {-999};
	int16_t ei{-999};
};

class EMCALChannelTiming : public SubsysReco
{
 public:

  EMCALChannelTiming(const std::string &name = "EMCALChannelTiming");

  ~EMCALChannelTiming() override;

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
//  int End(PHCompositeNode *topNode) override;

  /// Reset
//  int Reset(PHCompositeNode * /*topNode*/) override;

 // void Print(const std::string &what = "ALL") const override;

 private:
  float SubdivideDetector( std::vector<tower*>*, std::vector<tower*>*, PHCompositeNode*);
  void AnaHelper( std::vector<tower*>, std::vector<TH1F*>, std::vector<TH2F*>);
  std::string emcal_tower {""};
  std::string emcalgeom {""};

  enum 1DOUTPUTHISTS
  {
	 DELTAT, 
	 E,
	 EBART,
  };
  enum 2DOUTPUTHISTS
  {
	  EtoT,
	  EtaPhi,
  }
  std::vector<float> energy_bins {}; //log bins for energy
  
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
  std::vector< std::vector<TH1F*>*>* AllTowers1D_t=new std::vector< std::vector<TH1F*>*> {};
  std::vector< std::vector<TH1F*>*>* highTowers1D_t=new std::vector< std::vector<TH1F*>*> {};
  std::vector< std::vector<TH1F*>*>* lowTowers1D_t=new std::vector< std::vector<TH1F*>*> {};

  //2D
  std::vector< std::vector<TH2F*>*>* AllTowers2D_t=new std::vector< std::vector<TH2F*>*> {};
  std::vector< std::vector<TH2F*>*>* highTowers2D_t=new std::vector< std::vector<TH2F*>*> {};
  std::vector< std::vector<TH2F*>*>* lowTowers2D_t=new std::vector< std::vector<TH2F*>*> {};

};

#endif // EMCALCHANNELTIMING_H
