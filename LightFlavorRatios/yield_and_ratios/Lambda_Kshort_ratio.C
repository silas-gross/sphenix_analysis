#include "../corrections/EfficiencyCorrection.h"
#include "../corrections/LambdaFeedDownCorrection.h"
#include "../corrections/GeoAcceptanceCorrection.h"
#include "../corrections/TrivialEfficiencyCorrection.h"

#include "ResonanceRatio.h"
//#include "calculate_ratios.C"
#include "LambdaModel.h"
#include "KshortModel.h"

void Lambda_Kshort_ratio()
{
  TFile* Ks_file = TFile::Open("/sphenix/tg/tg01/hf/mjpeters/LightFlavorResults/Kshort_3runs.root");
  TFile* lambda_file = TFile::Open("/sphenix/tg/tg01/hf/mjpeters/LightFlavorResults/Lambda_3runs.root");

  TTree* Ks_tree = (TTree*)Ks_file->Get("DecayTree");
  TTree* lambda_tree = (TTree*)lambda_file->Get("DecayTree");

  std::vector<HistogramInfo> variables =
  {
    BinInfo::final_pt_bins,
    BinInfo::final_eta_bins,
    BinInfo::final_rapidity_bins,
    BinInfo::final_phi_bins,
  };

  RooArgList Ks_args;
  RooArgList lambda_args;

  RooRealVar m_ks("K_S0_mass","K_S0_mass",0.45,0.55);
  RooRealVar m_lambda("Lambda0_mass","Lambda0_mass",1.1,1.14);

  Ks_args.add(m_ks);
  lambda_args.add(m_lambda);

  std::vector<RooRealVar> Ks_vars;
  std::vector<RooRealVar> lambda_vars;

  for(HistogramInfo& hinfo : variables)
  {
    std::string Ks_branchname = "K_S0_"+hinfo.name;
    std::string lambda_branchname = "Lambda0_"+hinfo.name;
    std::cout << Ks_branchname << " " << lambda_branchname << std::endl;
    RooRealVar Ks_var(Ks_branchname.c_str(),Ks_branchname.c_str(),hinfo.bins.front(),hinfo.bins.back());
    RooRealVar lambda_var(lambda_branchname.c_str(),lambda_branchname.c_str(),hinfo.bins.front(),hinfo.bins.back());
    Ks_vars.push_back(Ks_var);
    lambda_vars.push_back(lambda_var);
  }

  for(int i=0; i<Ks_vars.size(); i++)
  {
    Ks_args.add(Ks_vars[i]);
    lambda_args.add(lambda_vars[i]);
  }

  Ks_args.Print();
  lambda_args.Print();

  RooDataSet Ks_ds("K_S0","K_S0",Ks_args,RooFit::Import(*Ks_tree));
  RooDataSet lambda_ds("Lambda0","Lambda0",lambda_args,RooFit::Import(*lambda_tree));

  KshortModel kshort_model;
  LambdaModel lambda_model;

  std::string fd_filename = "/sphenix/tg/tg01/hf/hjheng/HF-analysis/simulation/Pythia_ppMinBias/cascade_feeddown/Cascade_feeddown_fraction.root";
  std::vector<std::vector<std::shared_ptr<CorrectionHistogram1D>>> corrections(variables.size());
  // pT
  corrections[0].push_back(std::make_shared<LambdaFeedDownCorrection>(fd_filename,"h_feeddown_frac_xi_all"));
  corrections[0].push_back(std::make_shared<EfficiencyCorrection>());
  corrections[0].push_back(std::make_shared<GeoAcceptanceCorrection>("/sphenix/u/cdean/analysis/LightFlavorRatios/geometric_acceptance/analysis/plots/Lambda0_to_KS0_geometric_acceptance_ratio.root","pT"));

  // eta
  corrections[1].push_back(std::make_shared<LambdaFeedDownCorrection>(fd_filename,"h_feeddown_frac_xi_eta_all"));
  corrections[1].push_back(std::make_shared<TrivialEfficiencyCorrection>("awef"));
  corrections[1].push_back(std::make_shared<GeoAcceptanceCorrection>("/sphenix/u/cdean/analysis/LightFlavorRatios/geometric_acceptance/analysis/plots/Lambda0_to_KS0_geometric_acceptance_ratio_eta.root","Lambda0_inGeo_#eta"));

  // rapidity
  corrections[2].push_back(std::make_shared<LambdaFeedDownCorrection>(fd_filename,"h_feeddown_frac_xi_rapidity_all"));
  corrections[2].push_back(std::make_shared<TrivialEfficiencyCorrection>("aweg"));
  corrections[2].push_back(std::make_shared<GeoAcceptanceCorrection>("/sphenix/u/cdean/analysis/LightFlavorRatios/geometric_acceptance/analysis/plots/Lambda0_to_KS0_geometric_acceptance_ratio_rap.root","Lambda0_inGeo_y"));

  // phi
  corrections[3].push_back(std::make_shared<LambdaFeedDownCorrection>(fd_filename,"h_feeddown_frac_xi_phi_all"));
  corrections[3].push_back(std::make_shared<TrivialEfficiencyCorrection>("aweh"));
  corrections[3].push_back(std::make_shared<GeoAcceptanceCorrection>("/sphenix/u/cdean/analysis/LightFlavorRatios/geometric_acceptance/analysis/plots/Lambda0_to_KS0_geometric_acceptance_ratio_phi.root","Lambda0_inGeo_#phi"));

  TFile* fout = new TFile("fits.root","RECREATE");

  ResonanceRatio analyzer(lambda_model,kshort_model,
                          fout,"lambdaKsratio","#Lambda/2K_{S}^{0} ratio",1./2.,true,
                          variables,corrections);

  analyzer.calculate_ratios_unbinned(lambda_ds,Ks_ds);
}
