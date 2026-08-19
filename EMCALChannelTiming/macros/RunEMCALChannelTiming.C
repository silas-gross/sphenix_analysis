
#include "EMCALChannelTiming.h"

int RunEMCALChannelTiming(std::string CALOTowersfile="")
{
	Fun4AllServer* se = Fun4AllServer::instance();
	Fun4AllDstInputManager* in = new Fun4AllDstInputManager("in");
	in->AddFile(CaloTowersfile);
	se->RegisterInputManager(in);
	int seg = ;
	EMCALChannelTiming* ti = new EMCALChannelTiming(seg, "EMCT");
	se->RegisterSubsystem(ti);
	se->run();
	ti->Print();
	return 0;
};

