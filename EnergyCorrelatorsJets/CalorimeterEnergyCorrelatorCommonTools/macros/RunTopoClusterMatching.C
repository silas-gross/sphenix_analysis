#pragma once

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 00, 0)

#include <fstream>
#include <iostream>
#include <string>
#include <utility>

#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/SubsysReco.h>
#include <ffamodules/CDBInterface.h>
#include "Calo_CalibLocal.C"


#include <caloreco/RawClusterBuilderTopo.h>

#include <jetbackground/RetowerCEMC.h>

#include <topoclustermatching/TopoClusterMatching.h>

R__LOAD_LIBRARY(libfun4all.so);
R__LOAD_LIBRARY(libfun4allraw.so);
R__LOAD_LIBRARY(libcalo_io.so);
R__LOAD_LIBRARY(libffamodules.so);
R__LOAD_LIBRARY(libjetbase.so);
R__LOAD_LIBRARY(libjetbackground.so);
R__LOAD_LIBRARY(libg4dst.so);
R__LOAD_LIBRARY(libTopoClusterMatching.so);

void RunTopoClusterMatching(
    std::string caloFile,
    std::string globalFile,
    std::string jetFile,
    std::string truthFile,
    int n_files = 1)
{
  Fun4AllServer *se = Fun4AllServer::instance();

  auto *inCalo =
      new Fun4AllDstInputManager("InputManagerCalo");
  auto *inGlobal =
      new Fun4AllDstInputManager("InputManagerGlobal");
  auto *inJet =
      new Fun4AllDstInputManager("InputManagerJet");
  auto *inTruth =
      new Fun4AllDstInputManager("InputManagerTruth");

  if (n_files > 1)
  {
    inCalo->AddListFile(caloFile);
    inGlobal->AddListFile(globalFile);
    inJet->AddListFile(jetFile);
    inTruth->AddListFile(truthFile);
  }
  else
  {
    inCalo->AddFile(caloFile);
    inGlobal->AddFile(globalFile);
    inJet->AddFile(jetFile);
    inTruth->AddFile(truthFile);
  }

  se->registerInputManager(inCalo);
  se->registerInputManager(inGlobal);
  se->registerInputManager(inJet);
  se->registerInputManager(inTruth);

  int runnumber = 0;
  int seg = 0;

  if (n_files > 1)
  {
    std::ifstream inputList(caloFile);
    std::string firstFile;
    std::getline(inputList, firstFile);

    const auto runSegment = Fun4AllUtils::GetRunSegment(firstFile);
    runnumber = runSegment.first;
    seg = runSegment.second;
  }
  else
  {
    const auto runSegment = Fun4AllUtils::GetRunSegment(caloFile);
    runnumber = runSegment.first;
    seg = runSegment.second;
  }

  const std::string outfile = "topocluster_match_run_"+std::to_string(runnumber)+"_seg_"+std::to_string(seg)+".root";
  std::cout<<"Writes output to file : " <<outfile<<std::endl;
  auto *rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  CDBInterface::instance()->Verbosity(0);

  RawClusterBuilderTopo *clusterBuilder =
      new RawClusterBuilderTopo("HcalRawClusterBuilderTopo");

  clusterBuilder->Verbosity(0);
  clusterBuilder->set_nodename("TOPOCLUSTER_ALLCALO");
  clusterBuilder->set_enable_HCal(true);
  clusterBuilder->set_enable_EMCal(true);
  clusterBuilder->set_noise(0.0053, 0.0351, 0.0684);
  clusterBuilder->set_significance(4.0, 2.0, 1.0);
  clusterBuilder->allow_corner_neighbor(true);
  clusterBuilder->set_do_split(true);
  clusterBuilder->set_minE_local_max(1.0, 2.0, 0.5);
  clusterBuilder->set_R_shower(0.025);
  clusterBuilder->set_use_only_good_towers(true);
  clusterBuilder->set_absE(true);

  se->registerSubsystem(clusterBuilder);

  RetowerCEMC *retower = new RetowerCEMC("RetowerCEMV");
  retower->Verbosity(0);
  retower->set_towerinfo(true);
  retower->set_frac_cut(0.5);

  se->registerSubsystem(retower);

  /*
   * Use a positive minimum cluster pT.  The TopoClusterMatching
   * constructor uses this value when constructing logarithmic bins.
   */
  TopoClusterMatching *matching =
      new TopoClusterMatching(0.5, outfile, "TopoClusterMatching");

  se->registerSubsystem(matching);

  se->run(10);
  se->End();
  se->PrintTimer();

  delete se;

  gSystem->Exit(0);
}

#endif
