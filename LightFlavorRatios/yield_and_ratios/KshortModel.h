#ifndef KSHORTMODEL_H
#define KSHORTMODEL_H

#include "ParticleModel.h"

struct KshortModel : ParticleModel
{
  KshortModel()
  {
    name = "K_S0";

    setup_default_mass_nsignal_nbackground("K_S0","K_{S}^{0}",0.45,0.55);

    add_signal_parameter("ks_mean","mean",0.49,0.48,0.5);
    add_signal_parameter("ks_mean2","mean2",0.49,0.48,0.5);
    add_signal_parameter("ks_width","width",0.005,0.0001,0.2);
    add_signal_parameter("ks_width2","width2",0.005,0.0001,0.2);

    signal_gaus1 = std::make_shared<RooGaussian>("ks_gaus1","signal_gaus1",*mass,signal_parameters[0],signal_parameters[2]);
    signal_gaus2 = std::make_shared<RooGaussian>("ks_gaus2","signal_gaus2",*mass,signal_parameters[0],signal_parameters[3]);
    n_gaus1 = std::make_shared<RooRealVar>("ngaus1","gaus 1",100.,0.,1e9);
    n_gaus2 = std::make_shared<RooRealVar>("ngaus2","gaus 2",100.,0.,1e9);

    //signal_function = std::make_shared<RooGaussian>("ks_signal","signal",*mass,signal_parameters[0],signal_parameters[2]);
    signal_function = std::make_shared<RooAddPdf>("ks_signal","double gaussian signal",RooArgList(*signal_gaus1,*signal_gaus2),RooArgList(*n_gaus1,*n_gaus2));
/*
    background_parameters.emplace_back("ks_lin_coef","lin_coef",0.,0.,1.);
    background_parameters.emplace_back("ks_quad_coef","quad_coef",0.,0.,1.);
    background_parameters.emplace_back("ks_cubic_coef","cubic_coef",0.,0.,1.);
    background_parameters.emplace_back("ks_quartic_coef","quartic_coef",0.,0.,1.);
    background_parameters.emplace_back("ks_quintic_coef","quintic_coef",0.,0.,1.);

    background_function = std::make_shared<RooChebychev>("ks_bkg","background",*mass,RooArgList(background_parameters.begin(),background_parameters.end()));
*/

    add_background_parameter("ks_k","k",-1.,-10.,10.);

    background_function = std::make_shared<RooExponential>("ks_bkg","background",*mass,background_parameters[0]);

    fit_function = std::make_shared<RooAddPdf>("ks_model","combined signal and background",RooArgList(*signal_function,*background_function),RooArgList(*n_signal,*n_background));
  }

  std::shared_ptr<RooGaussian> signal_gaus1;
  std::shared_ptr<RooGaussian> signal_gaus2;
  std::shared_ptr<RooRealVar> n_gaus1;
  std::shared_ptr<RooRealVar> n_gaus2;
};

#endif
