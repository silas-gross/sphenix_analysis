#pragma once 

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 00, 0)

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TKey.h>
#include <TSystem.h>

#include <string>
#include <vector>
#include <format>
#include <ctime>
#include <map>
#include <algorithm>

#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

#include "Special_colors.h"

std::array<Int_t, 8> sp_cls {};
Int_t x {0}, y {0};
UInt_t w = 0;
UInt_t h = 0;

void ConfigureTLegend(TLegend* l, bool header=false)
{
	l->SetFillColor(0);
	l->SetFillStyle(0);
	l->SetBorderSize(0);
	l->SetTextSize(0.035f);
	if(header)
	{
		l->AddEntry("", "#bf{#it{sPHENIX}} Internal", "");
		l->AddEntry("", "Pythia8 p+p #sqrt{s} = 200 GeV","");

		time_t timestamp;
		time(&timestamp);
		l->AddEntry("", ctime(&timestamp), "");
	}
	return;
}

void MakeColorSet(Skaydis_colors* st)
{
	sp_cls[0]=TColor::GetColor(
			st->lesbian_Flag->at(0).at(0), 
			st->lesbian_Flag->at(0).at(1),
			st->lesbian_Flag->at(0).at(2)
		);
	sp_cls[1]=TColor::GetColor(
			st->lesbian_Flag->at(2).at(0), 
			st->lesbian_Flag->at(2).at(1),
			st->lesbian_Flag->at(2).at(2)
		);
	sp_cls[2]=TColor::GetColor(
			st->lesbian_Flag->at(3).at(0), 
			st->lesbian_Flag->at(3).at(1),
			st->lesbian_Flag->at(3).at(2)
		);
	sp_cls[3]=TColor::GetColor(
			st->lesbian_Flag->at(5).at(0), 
			st->lesbian_Flag->at(5).at(1),
			st->lesbian_Flag->at(5).at(2)
		);
	sp_cls[4]=TColor::GetFreeColorIndex();	
	new TColor(
			sp_cls[4],
			st->enby_Flag->at(0).at(0), 
			st->enby_Flag->at(0).at(1),
			st->enby_Flag->at(0).at(2)
		);
	sp_cls[5]=TColor::GetFreeColorIndex();
	new TColor(
			sp_cls[5],
			st->enby_Flag->at(1).at(0), 
			st->enby_Flag->at(1).at(1),
			st->enby_Flag->at(1).at(2)
		);
	sp_cls[6]=TColor::GetFreeColorIndex();
	new TColor(
			sp_cls[6],
			st->trans_Flag->at(0).at(0), 
			st->trans_Flag->at(0).at(1),
			st->trans_Flag->at(0).at(2)
		);
	sp_cls[7]=TColor::GetFreeColorIndex();
	new TColor(
			sp_cls[7],
			st->trans_Flag->at(1).at(0), 
			st->trans_Flag->at(1).at(1),
			st->trans_Flag->at(1).at(2)
		);
	return;
}

void MatchToAll(TH1F* Real, TH1F* All, std::string label)
{
	if(!Real || !All) return;
	TCanvas* c1 	= new TCanvas(label.c_str(), label.c_str(), x, y, 0.9*w, 0.9*h);
	c1->cd();
	TPad* p1	= new TPad("p1", "p1", 0, 0, 1, 0.33);
	TPad* p2	= new TPad("p2", "p2", 0, 0.35, 1,1);
	p2->cd();

	Real->SetLineColor(sp_cls[5]);
	All->SetLineColor(sp_cls[6]);
	Real->SetMarkerColor(sp_cls[5]);
	All->SetMarkerColor(sp_cls[6]);
	
	All->Draw("e1");
	All->GetYaxis()->SetRangeUser(Real->GetMean()*0.5, All->GetMaximum()*1.25);
	Real->Draw("e1 same");

	p2->SetLogx();
	p2->SetLogy();

	p1->cd();
	
	TH1F* ratio 	= (TH1F*) Real->Clone();
	ratio->Divide(All);
	ratio->SetLineColor(sp_cls[2]);
	ratio->SetMarkerColor(sp_cls[2]);
	ratio->Draw("e1");
	ratio->SetYTitle("Real / All ");
	p1->SetLogx();
	p1->SetLogy();

	c1->cd();
	p1->Draw();
	p2->Draw();

	p2->cd();
	
	TLegend* l 	= new TLegend(0.2, 0.7, 0.35, 0.9);
	TLegend* l1	= new TLegend(0.2, 0.5, 0.35, 0.7);
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
	l1->AddEntry(Real, Form("'Real' %s (Matched to at least one %s)", lb.c_str(), second_label.c_str()));
	l1->AddEntry(All, Form("All %s (passing Kinematic cuts", lb.c_str()));
	l1->AddEntry(ratio, Form("Real / All %s", lb.c_str()));
	l1->AddEntry("", "", "");
	float integral_rat = (float) Real->Integral() / (float) All->Integral();
	l1->AddEntry("", std::format("#int Real / #int All = {:.3f}", integral_rat).c_str(), "");
	l1->Draw();
	
	c1->cd();
	c1->Print(Form("~/%s_ET_funct.pdf", label.c_str()));

	return;

}
void CompareSingles(TFile* raw_results)
{
	if(!raw_results) return;

	TH1F* RealClusters 	= (TH1F*) raw_results->Get("h_RealCluster");
	TH1F* AllClusters 	= (TH1F*) raw_results->Get("h_AllCluster");
	TH1F* MatchTruth 	= (TH1F*) raw_results->Get("h_MatchedTruth");
	TH1F* AllTruth 		= (TH1F*) raw_results->Get("h_AllTruth");
	if(RealClusters && AllClusters) MatchToAll(RealClusters, AllClusters, "Clusters");
	if(MatchTruth && AllTruth) MatchToAll(MatchTruth, AllTruth, "Truth");
	return;
}

void ComparePairs(TFile* raw_results)
{
	if(!raw_results) return;
	TH1F* RealClusterPairs	= (TH1F*) raw_results->Get("h_ClusterPairET");
	TH1F* AllClusterPairs	= (TH1F*) raw_results->Get("h_ClusterPairETAll");
	TH1F* MatchTruthPairs	= (TH1F*) raw_results->Get("h_TruthPairET");
	TH1F* AllTruthPairs	= (TH1F*) raw_results->Get("h_TruthPairETAll");
	if(RealClusterPairs && AllClusterPairs) MatchToAll(RealClusterPairs, AllClusterPairs, "Cluster_pair");
	if(MatchTruthPairs && AllTruthPairs) MatchToAll(MatchTruthPairs, AllTruthPairs, "Truth_pair");
	return;	
}

static std::string GetSubBinLabel(const std::string& name)
{
	std::string label = name;
	std::replace(label.begin(), label.end(), '_', ' ');
	return label;
}

static void AddLegendToPad(TPad* pad, const std::vector<std::pair<std::string, int>>& entries, const std::string& title)
{
	if(!pad) return;
	pad->cd();
	TLegend* leg = new TLegend(0.15, 0.7, 0.45, 0.9);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.035f);
	leg->SetHeader(title.c_str());
	for(const auto& [label, color] : entries)
	{
		TLine line;
		line.SetLineColor(color);
		leg->AddEntry((TObject*)nullptr, label.c_str(), "l");
	}
	leg->Draw();
}
void findPairedHistos(std::vector<TH1F*> Real, std::vector<TH1F*> All, std::array<std::pair<TH1F*, TH1F*>, 5>* hists)
{
	for(int i= 0; i<(int) Real.size(); i++)
	{
		//grab the name for the histogram and find the matching tags
		std::string real_name = Real.at(i)->GetName();
		std::string tag_second {};
		tag_second = real_name.substr(real_name.r_find("_"));
		for(int j=0; j<(int) All.size(); j++)
		{
			std::string all_name = All.at(j)->GetName();
			std::string all_second = all_name.substr(all_name.r_find("_"));
			if(all_second == tag_second)
			{
				hists.at(i)=std::make_pair(Real.at(i), All.at(i));
				break;
			}
		}
	}
	return;
}
// Draw the histograms from one jet_dep/lead_bin_N directory on a single page.
// Each sub-bin is shown in a distinct color; matched histograms use solid
// lines, all-object histograms use dashed lines, and the lower row shows
// matched/all ratios.
void CompareLeadBin(TFile* raw_results)
{
	std::array<TDirectory*, 10> dirs {};
	if(!raw_results) return;
	for(int i=0; i<(int)dirs.size(); i++)
	{
		TDirectory* dir = raw_results->GetDirectory(Form("jet_dep/lead_bin_%d", i));
		if(!dir){
			dirs.at(i)=nullptr;
		       	continue;
		}
		dirs.at(i)=dir;
	}
	std::array<std::array<std::pair<TH1F*, TH1F*>, 5>, 10> cluster_histograms {}; 
	std::array<std::array<std::pair<TH1F*, TH1F*>, 5>, 10> truth_histograms {}; 
	std::array<std::array<std::pair<TH1F*, TH1F*>, 5>, 10> cluster_pair_histograms {}; 
	std::array<std::array<std::pair<TH1F*, TH1F*>, 5>, 10> truth_pair_histograms {}; 
	std::array<std::vector<std::string>, 10> keys {};
	for(int i=0; i<(int)dirs.size(); i++)
	{
		std::vector<std::string> dir_key {};
		if(!dirs.at(i)){
			keys.at(i)=dir_key;
			continue;
		}
		TIter next_key(dirs.at(i)->GetListOfKeys());
		TKey* key {nullptr};
		std::vector<TH1F*> SingleClusterReal 	{};
		std::vector<TH1F*> SingleClusterAll  	{};
		std::vector<TH1F*> SingleTruthMatch  	{};
		std::vector<TH1F*> SingleTruthAll	{};
		std::vector<TH1F*> PairClusterReal 	{};
		std::vector<TH1F*> PairClusterAll  	{};
		std::vector<TH1F*> PairTruthMatch  	{};
		std::vector<TH1F*> PairTruthAll		{};

		while( (key=(TKey*)(next_key())))
		{
			std::string keyN = key->GetName();
			TH1F* h = dirs->at(i)->Get(keyN.c_str());
			if(!h) continue;
			if(keyN.find("Pair") == std::string::npos)
			{
				//this collects all the singles
				if(keyN.find("ClET") != std::string::npos)
				{
					if(keyN.find("All") != std::string::npos) SingleClusterAll->push_back(h);
					else SingleClusterReal->push_back(h);
				}
				else
				{
					if(keyN.find("All") != std::string::npos) SingleTruthAll->push_back(h);
					else SingleTruthMatch->push_back(h);
				}
			}
			else
			{
				if(keyN.find("ClET") != std::string::npos)
				{
					if(keyN.find("All") != std::string::npos) PairClusterAll->push_back(h);
					else PairClusterReal->push_back(h);
				}
				else
				{
					if(keyN.find("All") != std::string::npos) PairTruthAll->push_back(h);
					else PairTruthMatch->push_back(h);
				}
			}
		}
		std::array<std::pair<TH1F*, TH1F*>,5> cls 	{};
		std::array<std::pair<TH1F*, TH1F*>,5> ts 	{};
		std::array<std::pair<TH1F*, TH1F*>,5> clp 	{};
		std::array<std::pair<TH1F*, TH1F*>,5> tp	{};
		findPairedHistos(SingleClusterReal, SingleClusterAll, &cls);
		findPairedHistos(SingleTruthReal, SingleTruthAll, &ts);
		findPairedHistos(PairClusterReal, PairClusterAll, &clp);
		findPairedHistos(PairTruthReal, PairTruthAll, &tp);
		cluster_histograms.at(i)=cls;
		truth_histograms.at(i)=ts;
		cluster_pair_histograms.at(i)=clp;
		truth_pair_histograms.at(i)=tp;

	}
	//I need to rethink the way this is structured and just do it as 4 pages of a grid
	/*
	struct HistogramPair
	{
	TIter next_key(dir->GetListOfKeys());
	TKey* key = nullptr;
	while((key = static_cast<TKey*>(next_key())))
		keys.push_back(key->GetName());
	std::sort(keys.begin(), keys.end());

	for(const auto& key_name : keys)
	{
		for(std::size_t prefix_index = 0; prefix_index < prefixes.size(); ++prefix_index)
		{
			const std::string& prefix = prefixes[prefix_index];
			if(key_name.compare(0, prefix.size(), prefix) != 0) continue;
			TH1F* hist = dynamic_cast<TH1F*>(dir->Get(key_name.c_str()));
			if(!hist) break;
			const std::string sub_bin = key_name.substr(prefix.size());
			HistogramPair& pair = histograms[categories[prefix_index]][sub_bin];
			if(is_all[prefix_index]) pair.all = hist;
			else pair.matched = hist;
			break;
		}
	}

	const std::array<const char*, 4> category_names {
		"Cluster E_{T}", "Cluster-pair E_{T}", "Truth E_{T}", "Truth-pair E_{T}"
	};*/
	TCanvas* truth_canvas = new TCanvas("Truth_Singles", "Truth Singles", 0, 0, 0.9*w, 0.9*h);
	TCanvas* cluster_canvas = new TCanvas("Clusters_Singles", "Cluster Singles", 0, 0, 0.9*w, 0.9*h);
	TCanvas* truth_pair_canvas = new TCanvas("Truth_Pair", "Truth Pair", 0, 0, 0.9*w, 0.9*h);
	TCanvas* cluster_pair_canvas = new TCanvas("Clusters_Pair", "Cluster Pair", 0, 0, 0.9*w, 0.9*h);

	truth_canvas->Divide(5, 10);
	cluster_canvas->Divide(5, 10);
	truth_pair_canvas->Divide(5, 10);
	cluster_pair_canvas->Divide(5, 10);
	

	return;
}

void DrawResults(TFile* raw_results)
{
	SetsPhenixStyle();
	Skaydis_colors* st = new Skaydis_colors();
	gStyle->SetPalette(100, st->Lesbian_log_gradient_PT);
	gVirtualX->GetWindowSize(gClient->GetRoot()->GetId(), x, y, w, h);
	MakeColorSet(st);
	CompareSingles(raw_results);
	ComparePairs(raw_results);
	CompareBins(raw_resutls);
	return;
}
#endif
