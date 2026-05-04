#ifndef LAMBDAMODEL_H
#define LAMBDAMODEL_H

#include "ParticleModel.h"

struct LambdaModel : ParticleModel
{
  LambdaModel()
  {
    name = "Lambda0";

    setup_default_mass_nsignal_nbackground("Lambda0","#Lambda^{0}",1.1,1.14);
/*
    signal_parameters.emplace_back("lambda_mean","mean",1.113,1.11,1.12);
    signal_parameters.emplace_back("lambda_width","width",0.01,0.0001,0.2);

    signal_function = std::make_shared<RooGaussian>("lambda_signal","signal",*mass,signal_parameters[0],signal_parameters[1]);

    background_parameters.emplace_back("lambda_lin_coef","lin_coef",0.,0.,1.);
    background_parameters.emplace_back("lambda_quad_coef","quad_coef",0.,0.,1.);
    background_parameters.emplace_back("lambda_cubic_coef","cubic_coef",0.,0.,1.);
    background_parameters.emplace_back("lambda_quartic_coef","quartic_coef",0.,0.,1.);
    background_parameters.emplace_back("lambda_quintic_coef","quintic_coef",0.,0.,1.);

    background_function = std::make_shared<RooChebychev>("lambda_bkg","background",*mass,RooArgList(background_parameters.begin(),background_parameters.end()));
*/
    add_signal_parameter("lambda_x0","x0",1.113,1.11,1.12);
    add_signal_parameter("lambda_sigmaL","sigma_L",0.01,0.0001,0.2);
    add_signal_parameter("lambda_sigmaR","sigma_R",0.01,0.0001,0.2);
    add_signal_parameter("lambda_alphaL","alpha_L",0.01,0.0001,100.);
    add_signal_parameter("lambda_nL","n_L",1.,0.0001,5.);
    add_signal_parameter("lambda_alphaR","alpha_R",0.01,0.0001,100.);
    add_signal_parameter("lambda_nR","n_R",1.,0.0001,5.);

    signal_function = std::make_shared<RooCrystalBall>("lambda_signal","signal",*mass,signal_parameters[0],signal_parameters[1],signal_parameters[2],signal_parameters[3],signal_parameters[4],signal_parameters[5],signal_parameters[6]);

    //fThresholdBackground = std::make_shared<TF1>("fThresholdBackground","TMath::Sqrt(x)*TMath::Exp([0]*x)",1.1,1.14);

    add_background_parameter("lambda_k","k",-1.,-10.,10.);

    //background_function = std::shared_ptr<RooAbsPdf>(RooFit::bindPdf(fThresholdBackground.get(),*mass));
    background_function = std::make_shared<RooExponential>("lambda_background","background",*mass,background_parameters[0]);

    fit_function = std::make_shared<RooAddPdf>("lambda_model","combined signal and background",RooArgList(*signal_function,*background_function),RooArgList(*n_signal,*n_background));
  }

  std::shared_ptr<TF1> fThresholdBackground;
};

#endif
