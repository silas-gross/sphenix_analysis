#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <string>
#include <sstream>
#include <iostream>
#include <format>

#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/SubsysReco.h>
#include <ffamodules/CDBInterface.h>
#include "Calo_CalibLocal.C"

#include <caloreco/RawClusterBuilderTopo.h>

#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <jetbase/JetCalib.h>
#include <jetbackground/FastJetAlgoSub.h>

#include <g4jets/TruthJetInput.h>
//#include <globalvertex/GlobalVertexReco.h>
//#include <GlobalVertex.h>                                                            
//#include <MbdDigitization.h>
//#include <MbdReco.h>

#include <jetbackground/TimingCut.h>

#include <jetbackground/RetowerCEMC.h>
#include <etashiftstudy/EtaShiftStudy.h>

R__LOAD_LIBRARY(libfun4all.so);
R__LOAD_LIBRARY(libfun4allraw.so);
R__LOAD_LIBRARY(libcalo_io.so);
R__LOAD_LIBRARY(libffamodules.so);
R__LOAD_LIBRARY(libjetbase.so);
R__LOAD_LIBRARY(libjetbackgroud.so);
R__LOAD_LIBRARY(libg4dst.so);
R__LOAD_LIBRARY(libEtaShiftStudy.so);

void RunEtaShiftStudy(std::string caloFile, std::string globalFile, std::string jetFile, std::string truthFile, std::string outDir, int n_files) 
{
	std::string generatorLabel {""};
	Fun4AllServer* se = Fun4AllServer::instance();
	int runnumber = 0;
  	int seg = 0;

	bool doSim = ( truthFile.find("none") != std::string::npos ) ? true : false ;
       	if (doSim)
       	{
		generatorLabel = ( truthFile.find("Herwig") != std::string::npos ) ? "Herwig" : "Pythia" ;
       	}
       	Fun4AllInputManager *inCalo = new Fun4AllDstInputManager("InputManagerCalo");
       	Fun4AllInputManager *inGlobal = new Fun4AllDstInputManager("InputManagerGlobal");
       	Fun4AllInputManager *inJet = new Fun4AllDstInputManager("InputManagerJet");
       	Fun4AllInputManager *inTruth = new Fun4AllDstInputManager("InputManagerTruth");
       	if(n_files > 1) 
       	{
	       	std::ifstream ifs(caloDSTlist);
	       	std::string filepath;
		std::getline(ifs,filepath);
		std::pair<int,int> runseg = Fun4AllUtils::GetRunSegment(filepath);
		runnumber = runseg.first;
		seg = runseg.second;

		inCalo->AddListFile(caloFile);
		inGlobal->AddListFile(globalFile);
		inJet->AddListFile(jetFile);
		if( doSim )	
			inTruth->addListFile(truthFile);


	}
	else
	{
	}
	se->registerInputManager(inCalo);
	se->registerInputManager(inGlobal);
	se->registerInputManager(inJet);
	if( doSim )
		se->registerInputManager(inTruth);
	
	auto rc = recoConsts::instance();
	rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
	if(doSim) 
		rc->set_uint64Flag("TIMESTAMP", 28);
	else 
		rc->set_uint64Flag("TIMESTAMP", runnumber);

	CDBInterface::instance()->Verbosity(0);
	Process_Calo_Calib();

	RawClusterBuilderTopo* ClusterBuilder = new RawClusterBuilderTopo("HcalRawClusterBuilderTopo");
	ClusterBuilder->Verbosity(0);
	ClusterBuilder->set_nodename("TOPOCLUSTER_ALLCALO");
	ClusterBuilder->set_enable_HCal(true);
	ClusterBuilder->set_enable_EMCal(true);
	ClusterBuilder->set_noise(0.0053, 0.0351, 0.0684); // 3sigma of pedestal noise
	ClusterBuilder->set_significance(4.0, 2.0, 1.0);
	ClusterBuilder->allow_corner_neighbor(true);
	ClusterBuilder->set_do_split(true);
	ClusterBuilder->set_minE_local_max(1.0, 2.0, 0.5);
	ClusterBuilder->set_R_shower(0.025);
	ClusterBuilder->set_use_only_good_towers(true);
	ClusterBuilder->set_absE(true);
	se->registerSubsystem(ClusterBuilder);

	RetowerCEMC* rtcemc = new RetowerCEMC("RetowerCEMV");
	rtcemc->Verbosity(0);
	rtcemc->set_towerinfo(true);
	rtcemc->set_frac_cut(0.5);
	se->registerSubsystem(rtcemc);
	
	std::string outfile=outDir+"/EtaShift"+generatorLabel+"-"+std::to_string(runnumber)+"_"+std::to_string(seg)+".root";
	EtaShiftStudy* ess = new EtaShiftStudy(outfile);
	se->registerSubsystem(ess);

	se->run(0);
	se->End();
	se->PrintTimer();
	delete se;
	gSystem->Exit(0);
}
#endif
