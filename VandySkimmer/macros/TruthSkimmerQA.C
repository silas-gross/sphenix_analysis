#pragma once 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <TTree.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TExec.h>

#include "Special_colors.h"
#include "PlottingConfigs.h"
#include "sPhenixStyle.h"
#include "sPhenixStyle.C"
#include "DrawingMacro.C"

#include <vandyclasses/Event.h>
#include <vandyclasses/JetInfo.h>
#include <vandyclasse/Tower.h>

#include <vector>
#include <array>
#include <math.h>
#include <string>
#include <numbers>

R__LOAD_LIBRARY(libVandyClasses.so);

PlottingConfigs* conf {nullptr};
Skaydis_colors* stylepoints {nullptr};

TExec* les {nullptr};
TExec* bi {nullptr};
TExec* trans {nullptr};
TExec* enby {nullptr};
std::string tag{""};
bool isData {false}; 

TH1F* wEECtruthtowers 	{nullptr};
TH1F* wEECtruthparticles{nullptr};
TH1F* wEECrecotowers	{nullptr};


int TruthSkimmerQA(std::string input_file_name, std::string outputQA_name, bool isD=false, std::string gen="Pythia")
{
	tag 		= gen;
	stylepoints 	= new Skaydis_colors();
	isData 		= isD;
	//setting up the pallets
	trans=new TExec("trans_pallet", "gStyle->SetPalette(100, style_points->Trans_gradient_PT);");
	bi=new TExec("bi_pallet", "gStyle->SetPalette(100, style_points->Trans_gradient_PT);");
	les=new TExec("les_pallet", "gStyle->SetPalette(100, style_points->Trans_gradient_PT);");
	enby=new TExec("enby_pallet", "gStyle->SetPalette(100, style_points->Trans_gradient_PT);");
	
	TFile* input_file = new TFile(input_file_name.c_str(), "READ");
	input_file->cd();
	//make the branches
	
	//truth jets
	JetInfo truthjet02 {nullptr};
	JetInfo truthjet03 {nullptr};
	JetInfo truthjet04 {nullptr};
	JetInfo truthjet05 {nullptr};
	//accompaning histograms 
	TH1F* truthjet02pt {nullptr};	
	TH1F* truthjet03pt {nullptr};	
	TH1F* truthjet04pt {nullptr};	
	TH1F* truthjet05pt {nullptr};	
	
	//tower jets
	JetInfo towerjet02 {nullptr};
	JetInfo towerjet03 {nullptr};
	JetInfo towerjet04 {nullptr};
	JetInfo towerjet05 {nullptr};
	//accompaning histograms 
	TH1F* towerjet02pt {nullptr};	
	TH1F* towerjet03pt {nullptr};	
	TH1F* towerjet04pt {nullptr};	
	TH1F* towerjet05pt {nullptr};	
	
	//truthtower jets
	JetInfo truthtowerjet02 {nullptr};
	JetInfo truthtowerjet03 {nullptr};
	JetInfo truthtowerjet04 {nullptr};
	JetInfo truthtowerjet05 {nullptr};
	//accompaning histograms 
	TH1F* truthtowerjet02pt {nullptr};	
	TH1F* truthtowerjet03pt {nullptr};	
	TH1F* truthtowerjet04pt {nullptr};	
	TH1F* truthtowerjet05pt {nullptr};	

	//to have and to hold
	std::array< std::pair<JetInfo, TH1F*>*, 4>* truthjets = new std::array<std::pair<JetInfo, TH1F*>, 4> {};
	std::array< std::pair<JetInfo, TH1F*>*, 4>* towerjets = new std::array<std::pair<JetInfo, TH1F*>, 4> {};
	std::array< std::pair<JetInfo, TH1F*>*, 4>* truthtowerjets = new std::array<std::pair<JetInfo, TH1F*>, 4> {};

	InitHistos(truthjets, truthjet02, truthjet02pt, 0, "truth");
	InitHistos(truthjets, truthjet03, truthjet03pt, 1, "truth");
	InitHistos(truthjets, truthjet04, truthjet04pt, 2, "truth");
	InitHistos(truthjets, truthjet05, truthjet05pt, 3, "truth");

	InitHistos(towerjets, towerjet02, towerjet02pt, 0, "reco_tower");
	InitHistos(towerjets, towerjet03, towerjet03pt, 1, "reco_tower");
	InitHistos(towerjets, towerjet04, towerjet04pt, 2, "reco_tower");
	InitHistos(towerjets, towerjet05, towerjet05pt, 3, "reco_tower");

	InitHistos(truthtowerjets, truthtowerjet02, truthtowerjet02pt, 0, "truth_tower");
	InitHistos(truthtowerjets, truthtowerjet03, truthtowerjet03pt, 1, "truth_tower");
	InitHistos(truthtowerjets, truthtowerjet04, truthtowerjet04pt, 2, "truth_tower");
	InitHistos(truthtowerjets, truthtowerjet05, truthtowerjet05pt, 3, "truth_tower");

	TTree* t = (TTree*) input_file->Get("T");
	SetJetBranches(t, truthjets, towerjets, truthtowerjets);
	SetTowerBranches(t, truthparticles, truthtowers, recotowers);
	EventInfo* event = new EventInfo();
	t->SetBranchAddress(&event, "EventInfo");
	int n = t->GetEntries();
	for(int i=0; i<n; i++) 
		eventloop(i, t, 
			truthjets, towerjets, truthtowerjets, 
			truthparticles, truthtowers, recotowers, 
			event);	
	DrawPlots(
			truthjets, towerjets, truthtowerjets, 
			truthparticles, truthtowers, recotowers
		 );
	return 0;
};
void InitJetHistos(
		std::array<std::pair<JetInfo, TH1F*>, 4>* holding_array,
		JetInfo* jt,
		TH1F* pt, 
		int type,
		std::string nametag,
		std::string tag
		)
{
	ha=new TH1F(std::format("h_{}_jet_r0{}", nametag, type+2).c_str(), 
			std::format("{} Jet R=0.{}; p_{{T}} [GeV]; #sigma", tag, type+2).c_str(),
		       25, 
	       		-0.5,
	 		99.5 );
	std::pair<JetInfo, TH1F*> ha = std::make_pair(*jt, pt);
	holding_array->at(type)=ha;
	return;	
}
void InitTowerHistos(

		)
{
}
void InitwEECHistos()
{
	wEECtruthtowers 	= new TH1F("h_truth_wEEC", "Truth tower wEEC; #Delta #varphi; wEEC";
	wEECtruthparticles	= new TH1F("h_part_wEEC", "Truth Particles wEEC; #Delta #varphi; wEEC";

	wEECrecotowers		= new TH1F("h_reco_wEEC", "Reco tower wEEC; #Delta #varphi; wEEC";

}
void SetJetBranches(
		TTree* t, 
		std::array<std::pair<JetInfo, TH1F*>*, 4>* truthjets, 
		std::array<std::pair<JetInfo, TH1F*>*, 4>* towerjets, 
		std::array<std::pair<JetInfo, TH1F*>*, 4>* truthtowerjets
		)
{
	if(!isData)
	{
		for(int i=0; i<(int)truthjets->size(); i++)
		{
			if(!truthjets->at(i)) continue;
			t->SetBranchAddress(&(truthjets->at(i)->first), std::format("TruthJetInfo_r0{}", i+2).c_str());
		}
		for(int i=0; i<(int)truthtowerjets->size(); i++)
		{
			if(!t->FindBranch(std::format("TruthTowerJetInfo_r0{}", i+2).c_str())) break;
			if(!truthtowerjets->at(i) ) continue;
			t->SetBranchAddress(&(truthtowerjets->at(i)->first), std::format("TruthTowerJetInfo_r0{}", i+2).c_str());
		}
	}
	for(int i=0; i<(int)towerjets->size(); i++)
	{
		if(!towerjets->at(i)) continue;
		t->SetBranchAddress(&(towerjets->at(i)->first), std::format("JetInfo_r0{}", i+2).c_str());
	}
	return;	
}
void SetTowerBranches(
		TTree* t,
		std::vector<Tower>* particle,
		std::vector<Tower>* truth_tower,
		std::vector<Tower>* reco_tower
		)
{
	if(isData)
	{
		if(particle)	t->SetBranchAddress("TruthParticles", &particle);
		if(truth_tower) t->SetBranchAddress("TruthTowers", &truth_tower);
	}
	if(reco_tower)	t->SetBranchAddress("TowerInfo", &reco_tower);
	return;
}
void eventLoop(
		int event_n,	
		TTree* t, 
		std::array<std::pair<JetInfo, TH1F*>*, 4>* 	truthjets, 
		std::array<std::pair<JetInfo, TH1F*>*, 4>* 	towerjets, 
		std::array<std::pair<JetInfo, TH1F*>*, 4>*	truthtowerjets,
		std::pair<std::vector<Tower>*, TH1F*>*	truthparticles,
		std::pair<std::vector<Tower>*, TH1F*>*	truthtowers,
		std::pair<std::vector<Tower>*, TH1F*>*	recotowers,
		EventInfo*					event
	      )
{
	t->GetEntry(event_n);
	bool isDijet = true;
	for(int i=0; i<(int)event->dijet_event->size(); i++)
	{
		if(!event->is_dijet_event(i)) isDijet=false;
	}
	if(!isDijet) return;

	if(!isData)
	{
		bool isTruthDijet = true;
		for(int i=0; i<(int)event->dijetTruth_event->size(); i++)
		{
			if(!event->is_dijetTruth_event(i)) isTruthDijet=false;
		}
		if(!isTruthDijet) return;

		for(int i=0; i<(int)truthjets->size(); i++)
		{
			if(!truthjets->at(i) || 
					!truthjets->at(i)->second
					) continue;
			truthjets->at(i)->second->Fill(truthjets->at(i)->first->pt());
		}
		getwEEC(truthparticles->first, wEECtruthparticles);

		for(int i=0; i<(int)truthtowerjets->size(); i++)
		{
			if(!t->FindBranch(std::format("TruthTowerJetInfo_r0{}", i+2).c_str())) break;
			if(!truthtowerjets->at(i) || 
					!truthtowerjets->at(i)->second
					) continue;
			truthtowerjets->at(i)->second->Fill(truthtowerjets->at(i)->first->pt());
		}

		for(int i=0; i<(int)truthtowers->first->size(); i++)
		{
			if(!truthtowers || 
					!truthtowers->first ||
					!truthtowers->second
			  ) continue;
			truthtowers->second->Fill(truthtowers->first->at(i)->e());
		}

		getwEEC(truthtowers->first, wEECtruthtowers);

	}

	for(int i=0; i<(int)towerjets->size(); i++)
	{
		if(!towerjets->at(i) || 
				!towerjets->at(i)->second
				) continue;
		towerjets->at(i)->second->Fill(towerjets->at(i)->first->pt());
	}
	for(int i=0; i<(int)recotowers->first->size(); i++)
	{
		if(!recotowers || 
				!recotowers->first ||
				!recotowers->second
		  ) continue;
		recotowers->second->Fill(recotowers->first->at(i)->e());
		getwEEC(recotowers->first, wEECrecotowers);
	}
	return;	
}

void getwEEC(std::vector<Tower>* towers, TH1F* output_hist, float Q2)
{
	for(int i=0; i<(int)towers->size()-1; i++)
	{
		float pt1 = std::pow(towers->at(i)->px(), 2) + std::pow(towers->at(i)->py(), 2);
		float phi1 = std::atan2(towers->at(i)->py(), towers->at(i)->px());
		pt1 = std::sqrt(pt1);
		for(int j=i+1; j<(int)towers->size(); j++)
		{
			float pt2 = std::pow(towers->at(j)->px(), 2) + std::pow(towers->at(j)->py(), 2);
			pt2 = std::sqrt(pt2);
			float phi2 = std::atan2(towers->at(j)->py(), towers->at(j)->px());
			float pw = pt1 * pt2/Q2;
			float dphi = phi1 - phi2;
			if(dphi < -pi) dphi+=2*pi;
			if(dphi > pi ) dphi-=2*pi;
			output_hist->Fill(dphi, pw);
		}
	}
	return;
}

#endif
