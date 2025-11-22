#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
//sphenix style
#include "sPhenixStyle.C"
#include "sPhenixStyle.h"

//root includes
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <TDirectory.h>
#include <TFitResult.h>
#include <TList.h>
#include <TLine.h>
#include <TBox.h>

//c++ includes
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <map>
#include <sstream>
std::map<std::string, int> color_map {{"truth", TColor::GetColor(252,106,223)}, {"reco", TColor::GetColor(91,206,250)}, {"emcal", TColor::GetColor(156,89,209)}, {"ihcal", TColor::GetColor(252,244,52)}, {"ohcal", TColor::GetColor(61,165,66)}};
void PlotJetCompDist(std::string input_type)
{
	SetsPhenixStyle();
	TDirectory* event_categorization=(TDirectory*)gFile->Get("event_categorization");
	event_categorization->cd();
	TH1F* h_jet_truth_R=(TH1F*)event_categorization->Get("h_jet_truth_R");	
	TH1F* h_total_E=(TH1F*)event_categorization->Get("h_total_E");	
	TCanvas* jet_truth=new TCanvas("jet_truth", "Jet energy component distribution");
	jet_truth->cd();
	h_jet_truth_R->Draw();
	jet_truth->Update();
	float y_up = jet_truth->GetUymax();
	TLine* l1=new TLine(h_jet_truth_R->GetBinCenter(11), 0, h_jet_truth_R->GetBinCenter(11), y_up);
	TLine* l2=new TLine(h_jet_truth_R->GetBinCenter(50), 0, h_jet_truth_R->GetBinCenter(50), y_up);
	TLine* l3=new TLine(h_jet_truth_R->GetBinCenter(60), 0, h_jet_truth_R->GetBinCenter(60), y_up);
	TLine* l4=new TLine(h_jet_truth_R->GetBinCenter(100), 0, h_jet_truth_R->GetBinCenter(100), y_up);
	TBox* b1=new TBox(l1->GetX1(), l1->GetY1(), l2->GetX2(), l2->GetY2());
	TBox* b2=new TBox(l3->GetX1(), l3->GetY1(), l4->GetX2(), l4->GetY2());
	l1->SetLineColor(color_map["truth"]);
	l2->SetLineColor(color_map["truth"]);
	l1->SetLineStyle(kDashDotted);
	l2->SetLineStyle(kDashDotted);
	b1->SetFillColor(color_map["truth"]);
	b1->SetFillStyle(3002);
	l3->SetLineColor(color_map["reco"]);
	l4->SetLineColor(color_map["reco"]);
	l3->SetLineStyle(kDashDotted);
	l4->SetLineStyle(kDashDotted);
	b2->SetFillColor(color_map["reco"]);
	b2->SetFillStyle(3002);
	l1->Draw("same");
	l2->Draw("same");
	b1->Draw("same");
	l3->Draw("same");
	l4->Draw("same");
	b2->Draw("same");
	TLegend* l=new TLegend(0, y_up*0.45, 1, y_up);
	l->SetFillStyle(0);
	l->SetFillColor(0);
	l->SetBorderSize(0);
	l->SetTextSize(0.03f);
	l->AddEntry("", "#bf{#it{sPHENIX}} Internal", "");
	l->AddEntry("", "p+p Pythia8 p_{T}^{lead} #geq 30 GeV", "");
	l->AddEntry("", "Pythia8 + GEANT + Noise", "");
	l->AddEntry("", input_type.c_str(), "");
	l->AddEntry("", Form("%d Dijet Events", (int)h_total_E->GetEntries()), "");
	l->AddEntry("", "p_{T}^{lead} #geq 20.9 GeV", "");
	l->AddEntry("", "p_{T}^{sublead} #geq 9.4 GeV", "");
	l->AddEntry("", "|#Delta #varphi|_{jet} #geq #frac{3 #pi}{4}", "");
	l->AddEntry("", "|#eta|_{jet} #leq 0.7", "");
	l->Draw();
	l->SetY1(y_up*0.45);
	l->SetY2(y_up);
	l->SetX1(-0.3);
	l->SetX2(1.5);
	l->Draw();
	TLegend* leg=new TLegend(h_jet_truth_R->GetBinCenter(55), y_up*0.45, h_jet_truth_R->GetBinCenter(75), y_up);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->AddEntry(b1, "Lead Jet Side", "f");
	float f1=h_jet_truth_R->Integral(11,50);
	int f1p=log10(f1);
	float f1r=f1*pow(10, -f1p);
	leg->AddEntry("", Form("#int E = %.3f #times 10^{%d} GeV",f1r, f1p), "");
	leg->AddEntry("", "", "");
	leg->AddEntry(b2, "Subleading Jet Side", "f");
	float f2=h_jet_truth_R->Integral(60,100);
	int f2p=log10(f2);
	float f2r=f2*pow(10, -f2p);
	leg->AddEntry("", Form("#int E = %.3f #times 10^{%d} GeV",f2r, f2p), "");
	leg->AddEntry("", "", "");
	leg->AddEntry("", Form("Sublead E/ Lead E = %.3f", h_jet_truth_R->Integral(60, 100)/(float) h_jet_truth_R->Integral(11, 50)), "");
	leg->Draw("same");
	leg->SetTextSize(0.03f);
	jet_truth->Update();
	l->SetY1(y_up*0.45);
	l->SetY2(y_up);
	l->SetX1(-0.3);
	l->SetX2(1.5);
	l->Draw();
	leg->SetX1(1.8);
	leg->SetX2(2.6);
	leg->SetY1(y_up*0.45);
	leg->SetY2(y_up*0.9);
	leg->SetY1(y_up*0.45);
	leg->SetY2(y_up*0.9);
	leg->Draw("same");
	std::cout<<leg->GetX1()<<", "<<leg->GetX2()<<", "<<leg->GetY1()<<", "<<leg->GetY2()<<std::endl;
	return;
}
#endif
