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
std::string global_tag {"Er_09-dR_01"};
std::string global_name {};

void ProcessGlobalTag()
{
	std::size_t Er_pos 	= global_tag.find( "_" );
	std::size_t sep		= global_tag.find( "-" );
	std::size_t dR_pos 	= global_tag.rfind( "_" );
	std::string Er_val 	= global_tag.substr( Er_pos + 2, Er_pos - sep - 2 );
	std::string dR_val	= global_tag.substr( dR_pos + 2 );
	global_name = "E_{Cluster} / E_{Truth} = 0."+Er_val + " , #Delta R_{Truth, Cluster} = 0."+dR_val;
	return;	
}	
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
		l->AddEntry("", global_name.c_str(), "");
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
	TCanvas* c1 	= new TCanvas(label.c_str(), label.c_str(), x, y, w, h);
	c1->cd();
	TPad* p1	= new TPad("p1", "p1", 0, 0, 1, 0.33);
	TPad* p2	= new TPad("p2", "p2", 0, 0.35, 1,1);
	p2->cd();

	int offset = 0; 
	if(label.find("Tr") != std::string::npos) offset=6;
	Real->SetLineColor(sp_cls[0+offset]);
	All->SetLineColor(sp_cls[(5+offset)%8]);
	Real->SetMarkerColor(sp_cls[0+offset]);
	All->SetMarkerColor(sp_cls[(5+offset)%8]);
	
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
	l1->AddEntry("", std::format("#int Real / #int All = {:.3e}", integral_rat).c_str(), "");
	l1->Draw();
	
	c1->cd();
	c1->Print(Form("~/%s_%s_ET_funct.pdf", global_tag.c_str(), label.c_str()));

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

std::string DropUnderscore(std::string name)
{
	std::string label = name;
	std::replace(label.begin(), label.end(), '_', ' ');
	return label;
}

void findPairedHistos(std::vector<TH1F*> Real, std::vector<TH1F*> All, std::array<std::pair<TH1F*, TH1F*>, 5>* hists)
{
	for(int i= 0; i<(int) Real.size(); i++)
	{
		//grab the name for the histogram and find the matching tags
		std::string real_name = Real.at(i)->GetName();
		std::string tag_second {};
		tag_second = real_name.substr(real_name.rfind("_"));
		for(int j=0; j<(int) All.size(); j++)
		{
			std::string all_name = All.at(j)->GetName();
			std::string all_second = all_name.substr(all_name.rfind("_")+1);
			if(all_second == tag_second)
			{
				hists->at(i)=std::make_pair(Real.at(i), All.at(i));
				break;
			}
		}
	}
	return;
}
void DrawBinComps(TCanvas* c, int canvas_numb, std::pair<TH1F*, TH1F*> hists)
{
	c->cd();
	c->cd(canvas_numb);
	TH1F* Real = hists.first;
	TH1F* All = hists.second;
	int offset=0;
	if(((std::string)Real->GetName()).find("Tr")!= std::string::npos) offset=6;
	if( !Real || !All) return;
	Real->SetLineColor(sp_cls[0+offset]);
	All->SetLineColor(sp_cls[(5+offset)%8]);
	Real->SetMarkerColor(sp_cls[0+offset]);
	All->SetMarkerColor(sp_cls[(5+offset)%8]);
	
	All->Draw("e1");
	All->GetYaxis()->SetRangeUser(Real->GetMean()*0.5, All->GetMaximum()*1.25);
	Real->Draw("e1 same");
	float all_int = All->Integral();
	float real_int = Real->Integral();
	real_int = real_int / all_int;
	TLegend* l1 = new TLegend(0.15, 0.7, 0.3, 0.9);
	ConfigureTLegend(l1);
	l1->AddEntry("", std::format("#int Matched / #int All = {:.3e}", real_int).c_str(), "");
	l1->Draw();
	c->cd();
	return;
}
void DrawBinned(TCanvas* c, std::array<std::array<std::pair<TH1F*, TH1F*>, 5>, 10> histos, std::string var_tag)
{
	std::array<std::string, 10> bin_label {};
	std::array<std::string, 5> subbin_label {}; 
	for(int i = 0; i<(int)histos.size(); i++)
	{
		for(int j=0; j<(int)histos.at(i).size(); j++)
		{
			if(!histos[i][j].first) continue;
			else{
				std::string hTitle = histos[i][j].first->GetTitle();
				std::size_t first_comma = hTitle.find(",");
				std::size_t second_comma = hTitle.find(",", first_comma);
				std::string tag_second = hTitle.substr(second_comma+1);
				if(subbin_label[j].find("p") == std::string::npos)
					subbin_label[j]=tag_second;
				std::string tag_first = hTitle.substr(first_comma, second_comma - first_comma - 1);
				if(bin_label[i].find("p") == std::string::npos) 
					bin_label[i]=tag_first;
			}
		}
	}
	c->cd();
	c->cd(1);
	TLegend* l1= new TLegend(0,0,1,1);
	ConfigureTLegend(l1, true);
	l1->SetTextSize(0.04f);
	l1->AddEntry("", DropUnderscore(var_tag).c_str(), "");
	c->cd();
	for(int i = 2; i<(int)bin_label.size() + 2; i++)
	{
		c->cd(i);
		TLegend* l = new TLegend(0,0,1,1);
		ConfigureTLegend(l);
		l->SetTextSize(0.05f);
		l->AddEntry("", bin_label[i-2].c_str(), "");
		l->Draw();
	}
	c->cd();
	for(int i = 0; i<(int)subbin_label.size(); i++)
	{
		int j = 1+(i+1)*12;
		c->cd(j);
		TLegend* l = new TLegend(0,0,1,1);
		ConfigureTLegend(l);
		l->SetTextSize(0.05f);
		l->AddEntry("", subbin_label[i].c_str(), "");
		l->Draw();
	}
	c->cd();
	for(int i=0; i<(int)histos.size(); i++)
	{
		for(int j=0; j<(int)histos[i].size(); j++)
		{
			int canvas_column = i+2;
			int canvas_row = j+1;
			int canvas_number = canvas_row*11 + canvas_column;
			DrawBinComps(c, canvas_number, histos[i][j]);
		}
	}
	c->Print(Form("~/%s_%s_binned.pdf", global_tag.c_str(), var_tag.c_str()));
	return;
}
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
			TH1F* h = (TH1F*)(dirs.at(i)->Get(keyN.c_str()));
			if(!h) continue;
			if(keyN.find("Pair") == std::string::npos)
			{
				//this collects all the singles
				if(keyN.find("ClET") != std::string::npos)
				{
					if(keyN.find("All") != std::string::npos) SingleClusterAll.push_back(h);
					else SingleClusterReal.push_back(h);
				}
				else
				{
					if(keyN.find("All") != std::string::npos) SingleTruthAll.push_back(h);
					else SingleTruthMatch.push_back(h);
				}
			}
			else
			{
				if(keyN.find("ClET") != std::string::npos)
				{
					if(keyN.find("All") != std::string::npos) PairClusterAll.push_back(h);
					else PairClusterReal.push_back(h);
				}
				else
				{
					if(keyN.find("All") != std::string::npos) PairTruthAll.push_back(h);
					else PairTruthMatch.push_back(h);
				}
			}
		}
		std::array<std::pair<TH1F*, TH1F*>,5> cls 	{};
		std::array<std::pair<TH1F*, TH1F*>,5> ts 	{};
		std::array<std::pair<TH1F*, TH1F*>,5> clp 	{};
		std::array<std::pair<TH1F*, TH1F*>,5> tp	{};
		findPairedHistos(SingleClusterReal, SingleClusterAll, &cls);
		findPairedHistos(SingleTruthMatch, SingleTruthAll, &ts);
		findPairedHistos(PairClusterReal, PairClusterAll, &clp);
		findPairedHistos(PairTruthMatch, PairTruthAll, &tp);
		cluster_histograms.at(i)=cls;
		truth_histograms.at(i)=ts;
		cluster_pair_histograms.at(i)=clp;
		truth_pair_histograms.at(i)=tp;

	}
	TCanvas* truth_canvas = new TCanvas("Truth_Singles", "Truth Singles", 0, 0, w, h);
	TCanvas* cluster_canvas = new TCanvas("Clusters_Singles", "Cluster Singles", 0, 0, w, h);
	TCanvas* truth_pair_canvas = new TCanvas("Truth_Pair", "Truth Pair", 0, 0, w, h);
	TCanvas* cluster_pair_canvas = new TCanvas("Clusters_Pair", "Cluster Pair", 0, 0, w, h);

	truth_canvas->Divide(6, 11);
	cluster_canvas->Divide(6, 11);
	truth_pair_canvas->Divide(6, 11);
	cluster_pair_canvas->Divide(6, 11);
	DrawBinned(truth_canvas, truth_histograms, "Single_Truth");
	DrawBinned(cluster_canvas, cluster_histograms, "Single_Cluster");
	DrawBinned(truth_pair_canvas, truth_pair_histograms, "Pair_Truth");
	DrawBinned(cluster_pair_canvas, cluster_pair_histograms, "Pair_Cluster");

	return;
}

void DrawResults(TFile* raw_results, std::string tag = "Er_09-dR_01")
{
	SetsPhenixStyle();
	global_tag = tag;
	ProcessGlobalTag();
	Skaydis_colors* st = new Skaydis_colors();
	gStyle->SetPalette(100, st->Lesbian_log_gradient_PT);
	gVirtualX->GetWindowSize(gClient->GetRoot()->GetId(), x, y, w, h);
	MakeColorSet(st);
	CompareSingles(raw_results);
	ComparePairs(raw_results);
	CompareLeadBin(raw_results);
	return;
}
#endif
