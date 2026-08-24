#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <string>
#include <iostream>
#include <format>

#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/SubsysReco.h>
#include <ffamodules/CDBInterface.h>
#include <Calo_Calib.C>

#include <caloreco/RawClusterBuilderTopo.h>


#include <emcalchanneltiming/EMCALChannelTiming.h>

//#include <jetbackground/RetowerCEMC.h>


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libcalo_io.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libEMCALChannelTiming.so)


int RunEMCALChannelTiming(std::string CALOTowersfile="", int nevt=0)
{
	Fun4AllServer* se = Fun4AllServer::instance();
	Fun4AllDstInputManager* in = new Fun4AllDstInputManager("in");
	in->AddFile(CALOTowersfile);
	se->registerInputManager(in);
	
	std::pair<int, int> runseg= Fun4AllUtils::GetRunSegment(CALOTowersfile);


	auto rc = recoConsts::instance();
	rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
	rc->set_uint64Flag("TIMESTAMP", runseg.first);
	CDBInterface::instance()->Verbosity(0);
	Process_Calo_Calib();

	EMCALChannelTiming* ti = new EMCALChannelTiming(runseg.second, runseg.first, "EMCT");
	se->registerSubsystem(ti);
	se->run(nevt);
	se->End();
	delete se;
	return 0;
}
#endif
