// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: Prototype4DSTReader.h,v 1.7 2015/02/27 23:42:23 jinhuang Exp $

/*!
 * \file Prototype4DSTReader.h
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.7 $
 * \date $Date: 2015/02/27 23:42:23 $
 */

#ifndef PROTOTYPE4_PROTOTYPE4DSTREADER_H
#define PROTOTYPE4_PROTOTYPE4DSTREADER_H

#include "RawTower_Prototype4.h"
#include "RawTower_Temperature.h"

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>
#include <vector>

class PHCompositeNode;
class TClonesArray;
class TTree;

/*!
 * \brief Prototype4DSTReader save information from DST to an evaluator, which
 * could include hit. particle, vertex, towers and jet (to be activated)
 */
class Prototype4DSTReader : public SubsysReco
{
 public:
  Prototype4DSTReader(const std::string &filename);
  virtual ~Prototype4DSTReader();

  //! full initialization
  int Init(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  void AddTower(const std::string &name) { _tower_postfix.push_back(name); }

  void AddTowerTemperature(const std::string &name)
  {
    _towertemp_postfix.push_back(name);
  }

  void AddRunInfo(const std::string &name) { _runinfo_list.push_back(name); }
  void AddEventInfo(const std::string &name)
  {
    _eventinfo_list.push_back(name);
  }

  //! zero suppression for all calorimeters
  double get_tower_zero_sup() { return _tower_zero_sup; }

  //! zero suppression for all calorimeters
  void set_tower_zero_sup(double b) { _tower_zero_sup = b; }

 protected:
  //  std::vector<std::string> _node_postfix;
  std::vector<std::string> _tower_postfix;
  //! tower temperature
  std::vector<std::string> _towertemp_postfix;
  //  std::vector<std::string> _jet_postfix;
  //  std::vector<std::string> _node_name;
  std::vector<std::string> _runinfo_list;
  std::vector<std::string> _eventinfo_list;

  int nblocks;

  typedef std::shared_ptr<TClonesArray> arr_ptr;

  struct record
  {
    unsigned int _cnt;
    std::string _name;
    arr_ptr _arr;
    TClonesArray *_arr_ptr;
    double _dvalue;

    enum enu_type
    {
      typ_hit,
      typ_part,
      typ_vertex,
      typ_tower,
      typ_jets,
      typ_runinfo,
      typ_eventinfo,
      typ_towertemp
    };
    enu_type _type;
  };
  typedef std::vector<record> records_t;
  records_t _records;

  typedef RawTower_Prototype4 RawTower_type;

  typedef RawTower_Temperature RawTowerT_type;

  int _event;

  std::string _out_file_name;

  //  TFile * _file;
  TTree *_T;

  //! zero suppression for all calorimeters
  double _tower_zero_sup;

  void build_tree();
};

#endif
