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
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

 private:
};

#endif // JETSHAPEEVENTOBSERVABLES_H
