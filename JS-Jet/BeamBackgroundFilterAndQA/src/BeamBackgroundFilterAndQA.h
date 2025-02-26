/// ===========================================================================
/*! \file    BeamBackgroundFilterAndQA.h
 *  \authors Hanpu Jiang, Derek Anderson
 *  \date    10.21.2024
 *
 *  A F4A module to filter out events with significant
 *  beam background (the so-called "streaky events")
 *  and produce some relevant QA histograms.
 */
/// ===========================================================================

#ifndef BEAMBACKGROUNDFILTERANDQA_H
#define BEAMBACKGROUNDFILTERANDQA_H

// module components
#include "BaseBeamBackgroundFilter.h"
#include "NullFilter.h"
#include "StreakSidebandFilter.h"

// f4a libraries
#include <fun4all/SubsysReco.h>

// phparameters libraries
#include <phparameter/PHParameters.h>

// c++ utilities
#include <map>
#include <memory>
#include <string>
#include <vector>

// forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class QAHistManagerHistDef;
class TowerInfoContainer;



// ============================================================================
//! Filter beam background events and create QA
// ============================================================================
/*! A F4A module to filter out events with significant
 *  beam background and produce some relevant QA
 *  histograms. 
 */
class BeamBackgroundFilterAndQA : public SubsysReco {

  public:

    // ========================================================================
    //! User options for module
    // =======================================================================
    struct Config
    {

      // turn modes on/off
      bool debug      = true;
      bool doQA       = true;
      bool doEvtAbort = false;

      ///! module name
      std::string moduleName = "BeamBackgroundFilterAndQA";

      ///! flag prefix
      std::string flagPrefix = "HasBeamBackground";

      ///! histogram tags
      std::string histTag = "";

      ///! which filters to apply
      std::vector<std::string> filtersToApply = {"Null", "StreakSideband"};

      ///! filter configurations
      NullFilter::Config null;
      StreakSidebandFilter::Config sideband;
      //... add other configurations here ...//

    };

    // ctor/dtor
    BeamBackgroundFilterAndQA(const std::string& name = "BeamBackgroundFilterAndQA", const bool debug = false);
    BeamBackgroundFilterAndQA(const Config& config); 
    ~BeamBackgroundFilterAndQA() override;

    // setters
    void SetConfig(const Config& config) {m_config = config;}

    // getters
    Config GetConfig() const {return m_config;}

    // f4a methods
    int Init(PHCompositeNode* topNode) override;
    int process_event(PHCompositeNode* topNode) override;
    int End(PHCompositeNode* /*topNode*/) override;

  private:

    // private methods
    void InitFilters();
    void InitFlags(PHCompositeNode* topNode);
    void InitHistManager();
    void BuildHistograms();
    void RegisterHistograms();
    void SetDefaultFlags();
    void UpdateFlags(PHCompositeNode* topNode);
    bool ApplyFilters(PHCompositeNode* topNode);
    std::string MakeFlagName(const std::string& filter = "");

    ///! histogram manager
    Fun4AllHistoManager* m_manager;

    ///! background flags
    PHParameters m_flags;

    ///! module-wide histograms
    std::map<std::string, TH1*> m_hists;

    ///! module configuration
    Config m_config;

    ///! filters
    std::map<std::string, std::unique_ptr<BaseBeamBackgroundFilter>> m_filters;

};  // end BeamBackgroundFilterAndQA

#endif

// end ========================================================================
