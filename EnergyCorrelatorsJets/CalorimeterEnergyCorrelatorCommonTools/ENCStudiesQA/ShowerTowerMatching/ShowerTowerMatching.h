// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SHOWERTOWERMATCHING_H
#define SHOWERTOWERMATCHING_H
//fun4all basic stuff
#include <fun4all/SubsysReco.h>

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
};

#endif // SHOWERTOWERMATCHING_H
