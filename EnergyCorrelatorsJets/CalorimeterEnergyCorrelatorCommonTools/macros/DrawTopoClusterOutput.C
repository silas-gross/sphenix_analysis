#pragma once 

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 00, 0)

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>

#include <string>
#include <vector>
#include <format>
#include <ctime>

#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

#include "Special_colors.h"

std::array<TColor*, 8> sp_cls {};

void ConfigureTLegend(TLegend* l, bool header=false)
{
	l->SetFillColor(0);
	l->SetFillStyle(0);
	l->SetBorderSize(0);
	if(header)
	{
		l->AddEntry("", "#bf{#it{sPHENIX}} Internal", "");
		l->AddEntry("", "Pythia8 p+p #sqrt{s} = 200 GeV","");

		time_t timestamp;
		time(&timestamp);
		l->AddEntry("", ctime(&timestamp).c_str(), "");
	}
	return;
}

void MakeColorSet()
{
	sp_cls[0]=new TColor(
			st->lesbian_Flag->at(0).at(0), 
			st->lesbian_Flag->at(0).at(1),
			st->lesbian_Flag->at(0).at(2)
		);
	sp_cls[1]=new TColor(
			st->lesbian_Flag->at(2).at(0), 
			st->lesbian_Flag->at(2).at(1),
			st->lesbian_Flag->at(2).at(2)
		);
	sp_cls[2]=new TColor(
			st->lesbian_Flag->at(3).at(0), 
			st->lesbian_Flag->at(3).at(1),
			st->lesbian_Flag->at(3).at(2)
		);
	sp_cls[3]=new TColor(
			st->lesbian_Flag->at(5).at(0), 
			st->lesbian_Flag->at(5).at(1),
			st->lesbian_Flag->at(5).at(2)
		);
	sp_cls[4]=new TColor(
			st->enby_Flag->at(0).at(0), 
			st->enby_Flag->at(0).at(1),
			st->enby_Flag->at(0).at(2)
		);
	sp_cls[5]=new TColor(
			st->enby_Flag->at(1).at(0), 
			st->enby_Flag->at(1).at(1),
			st->enby_Flag->at(1).at(2)
		);
	sp_cls[6]=new TColor(
			st->trans_Flag->at(0).at(0), 
			st->trans_Flag->at(0).at(1),
			st->trans_Flag->at(0).at(2)
		);
	sp_cls[7]=new TColor(
			st->trans_Flag->at(1).at(0), 
			st->trans_Flag->at(1).at(1),
			st->trans_Flag->at(1).at(2)
		);
	return;
}

void MatchToAll(TH1F* Real, TH1F* All, std::string label)
{
	TCanvas* c1 	= new TCanvas();
	c1->cd();
	TPad* p1	= new TPad("p1", "p1", 0, 0, 1, 0.33);
	TPad* p2	= new TPad("p2", "p2", 0, 0.35, 1,1);
	p2->cd();

	Real->SetLineColor(sp_cls[0]);
	All->SetLineColor(sp_cls[6]);
	Real->SetMarkerColor(sp_cls[0]);
	All>SetMarkerColor(sp_cls[6]);
	
	AllClusters->Draw();
	RealClusters->Draw("same");

	p2->SetLogx();
	p2->SetLogy();

	p1->cd();
	
	TH1F* ratio 	= (TH1F*) RealClusters->Clone();
	ratio->Divide(AllClusters);
	ratio->SetLineColor(sp_cls[2]);
	ratio->SetMarkerColor(sp_cls[2]);
	ratio->Draw();
	ratio->SetYTitle("Real / All ");

	c1->cd();
	p1->Draw();
	p2->Draw();

	p2->cd();
	
	TLegend* l 	= new TLegend(0.8, 0.8, 1, 1);
	TLegend* l1	= new TLegend(0.8, 0.4, 1, 0.6);
	ConfigureTLegend(l, true);
	ConfigureTLegend(l1);
	std::string lb {label};	
	std::string second_label="Truth particle";
	
	if(label.find("Truth") != std::string::npos)
	{
		second_label="Cluster";
		lb += " Particles";
	}
	
	if(label.find("_pair") != std::string::npos)
	{
		if( label.find("luster") != std::string::npos) lb="Cluster Pairs";
		else lb = "Truth Pairs";
	}
	
	l->AddEntry("", Form("%s E_{T}", lb.c_str()), "");
	l->Draw();
	l1->AddEntry(RealClusters, Form("'Real' %s (Matched to at least one %s)", lb.c_str(), second_label.c_str()));
	l1->AddEntry(AllClusters, Form("All %s (passing Kinematic cuts", lb.c_str()));
	l1->AddEntry(ratio, Form("Real / All %s", lb.c_str()));
	l1->Draw();
	
	c1->cd();
	c1->Print(Form("~/%s_ET_funct.pdf", label.c_str()));

	return;

}
void CompareSingles(TFile* raw_results)
{

	TH1F* RealClusters 	= (TH1F*) raw_results->Get("h_RealCluster");
	TH1F* AllClusters 	= (TH1F*) raw_results->Get("h_AllCluster");
	TH1F* MatchTruth 	= (TH1F*) raw_results->Get("h_MatchedTruth");
	TH1F* AllTruth 		= (TH1F*) raw_results->Get("h_AllTruth");
	MatchToAllSingle(RealClusters, AllClusters, "Clusters");
	MatchToAllSingle(MatchTruth, AllTruth, "Truth");
	return;
}

void ComparePairs(TFile* raw_results)
{
	TH1F* RealClusterPairs	= (TH1F*) raw_results->Get("h_ClusterPairET");
	TH1F* AllClusterPairs	= (TH1F*) raw_results->Get("h_ClusterPairETAll");
	TH1F* TruthPairs	= (TH1F*) raw_results->Get("h_TruthPairET");
	TH1F* AllTruthPairs	= (TH1F*) raw_results->Get("h_TruthPairETAll");
	MatchToAllSingle(RealClustersPairs, AllClustersPairs, "Cluster_pair");
	MatchToAllSingle(MatchTruthPairs, AllTruthPairs, "Truth_pair");
	return;	
}
void DrawResults(TFile* raw_results)
{
	SetsPhenixStyle();
	Skaydis_colors* st = new Skaydis_colors();
	gStyle->SetPalette(100, st->Lesbian_log_gradient_PT);
	MakeColorSet();
	CompareSingles(raw_results);
	ComparePairs(raw_results);
	return;
}
#endif;
