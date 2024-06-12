//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in CalorimeterCorrelationsPP.h.
//
// CalorimeterCorrelationsPP(const std::string &name = "CalorimeterCorrelationsPP")
// everything is keyed to CalorimeterCorrelationsPP, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// CalorimeterCorrelationsPP::~CalorimeterCorrelationsPP()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int CalorimeterCorrelationsPP::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int CalorimeterCorrelationsPP::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int CalorimeterCorrelationsPP::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT; 
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT; 
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int CalorimeterCorrelationsPP::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int CalorimeterCorrelationsPP::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int CalorimeterCorrelationsPP::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int CalorimeterCorrelationsPP::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void CalorimeterCorrelationsPP::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

#include "CalorimeterCorrelationsPP.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

//____________________________________________________________________________..
CalorimeterCorrelationsPP::CalorimeterCorrelationsPP(const std::string &name):
 SubsysReco(name)
{
  std::cout << "CalorimeterCorrelationsPP::CalorimeterCorrelationsPP(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
CalorimeterCorrelationsPP::~CalorimeterCorrelationsPP()
{
  std::cout << "CalorimeterCorrelationsPP::~CalorimeterCorrelationsPP() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int CalorimeterCorrelationsPP::Init(PHCompositeNode *topNode)
{
  std::cout << "CalorimeterCorrelationsPP::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CalorimeterCorrelationsPP::InitRun(PHCompositeNode *topNode)
{
  std::cout << "CalorimeterCorrelationsPP::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CalorimeterCorrelationsPP::process_event(PHCompositeNode *topNode)
{
	n_evt++;
  	std::cout << "CalorimeterCorrelationsPP::process_event(PHCompositeNode *topNode) Processing Event " <<n_evt << std::endl;
	//For the time being, need to do the calo tower wbuilder myself
	
	TowerInfoContainerv2* IHtowers=NULL, *OHtowers=NULL, *EMtowers=NULL;
  	std::string IHCALgeom="TOWERGEOM_HCALIN", OHCALgeom="TOWERGEOM_HCALOUT", EMCALgeom="TOWERGEOM_CEMC";
  	IHtowers=findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALIN");
  	OHtowers=findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALOUT");
  	EMtowers=findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_CEMC");
	EMclusters=findNode::getClass<>(topNode, "CLUSTERINFO_POS_COR_CEMC");
	RawTowerGeomContainer_Cylinderv1 *ihg=findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, ihcalgeom);
	RawTowerGeomContainer_Cylinderv1 *ohg=findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, ohcalgeom);
	RawTowerGeomContainer_Cylinderv1 *emg=findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, emcalgeom);
	GlobalVertexMapv1* vtxmap=findNode::getClass<GlobalVertexMapv1>(topNode, "GlobalVertexMap");
	if(vtxmap->empty()) std::cout<<"Vertex map is empty" <<std::endl;
	GlobalVertex* vtx=nullptr;
	if(vtxmap) for(GlobalVertexMap::ConstIter vtx_iter = vtxmap->begin(); vtx_iter != vtxmap->end(); ++vtd_iter) vtx=vtx_iter->second();
	h_vtx_z->Fill(vtx->get_z());
	float rvtx=sqrt(pow(vtx->get_x(), 2) + pow(vtx->get_y(),2));
	float phivtx=atan2(vtx->get_y(), vtx->get_x());
	h_vtx_2d->Fill(vtx->get_x(), vtx->get_y());
	h_vtx_rz->Fill(rvtx, vtx->get_z());
	h_vtx_phi->Fill(phivtx);

	
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CalorimeterCorrelationsPP::ResetEvent(PHCompositeNode *topNode)
{
  std::cout << "CalorimeterCorrelationsPP::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CalorimeterCorrelationsPP::EndRun(const int runnumber)
{
  std::cout << "CalorimeterCorrelationsPP::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CalorimeterCorrelationsPP::End(PHCompositeNode *topNode)
{
  std::cout << "CalorimeterCorrelationsPP::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CalorimeterCorrelationsPP::Reset(PHCompositeNode *topNode)
{
 std::cout << "CalorimeterCorrelationsPP::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void CalorimeterCorrelationsPP::Print(const std::string &what) const
{
  std::cout << "CalorimeterCorrelationsPP::Print(const std::string &what) const Printing info for " << what << std::endl;
}
