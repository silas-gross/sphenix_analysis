// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SILICONCALOMATCHING_H
#define SILICONCALOMATCHING_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxVertexMap;
class ActsGeometry;
class TrkrClusterContainer;
class RawClusterContainer;

class SvtxTrack;
class SvtxTrackState;
class SvtxCluster;

class RawCluster;
class RawClusterContainer;

class SiliconCaloTrack;
class SiliconCaloTrackMap;


class SiliconCaloMatching : public SubsysReco
{
public:
  SiliconCaloMatching(const std::string &name = "SiliconCaloMatching");

  ~SiliconCaloMatching() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode *topNode) override;

  void setMC(const bool input) { isMC = input; }

  void setTrackMapName( const std::string &name) { _trackMapName  = name; }
  void setVertexMapName(const std::string &name) { _vertexMapName = name; }
  void setClusterContainerName(     const std::string &name) { _clusterContainerName = name; }
  void setEMCalClusterContainerName(const std::string &name) { _emcalClusName = name; }

  void setEmcalRadius(float radius) { _caloRadiusEMCal = radius; }
  void setEmcalLowEcut(float ecut)  { _emcal_low_cut = ecut; }

/*
  void setTopoCluster(bool topo)
  {
    if (topo)
    {
      setClusterContainerName("TOPOCLUSTER_EMCAL");
    }
  }
*/

protected:
  std::string _clusterContainerName = "TRKR_CLUSTER";
  std::string _actsgeometryName     = "ActsGeometry";
  std::string _trackMapName         = "SvtxTrackMap";
  std::string _vertexMapName        = "SvtxVertexMap";
  std::string _emcalClusName        = "CLUSTER_CEMC";

  double _emcal_low_cut      = 0.18;
  float  _caloRadiusEMCal    = 93.5;
  float  _caloThicknessEMCal = 20.4997;

  bool isMC = false;

private:
  bool  getNodes(PHCompositeNode *topNode);
  void  getCaloPosition(RawCluster* calo, float& x, float& y, float& z);

  float calculatePt(SvtxTrack* track, RawCluster* emc);

private:
  SvtxTrackMap*         _trackmap   = nullptr;
  SvtxVertexMap*        _vertexmap  = nullptr;
  ActsGeometry*         _geometry   = nullptr;
  TrkrClusterContainer* _clustermap = nullptr;
  RawClusterContainer*  _emcClusmap = nullptr;

  SiliconCaloTrackMap*  _sicalomap  = nullptr;

};

#endif  // SILICONSEEDSANA_H

