#pragma once 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include "sPhenixStyle.h"
#include "sPhenixStyle.C"
#include "Special_colors.h"
#include <TLegend.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TMultiGraph.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <ctime>
#include <iomanip>

struct LeakageDataPoint
{
	std::string timestamp;
	float imeas_emcal;
	float imeas_ihcal;
	float imeas_ohcal;
};

std::vector<LeakageDataPoint> ReadLeakageData(const std::string& filename)
{
	std::vector<LeakageDataPoint> data;
	std::ifstream infile(filename);
	
	if (!infile.is_open())
	{
		std::cerr << "Error: Could not open file " << filename << std::endl;
		return data;
	}
	
	std::string line;
	// Skip header line
	std::getline(infile, line);
	int nskip = 24;
	int n=0;
	while (std::getline(infile, line))
	{
		if (line.empty()) continue;
		n++;
		if(n % nskip != 0) continue;
		
		std::stringstream ss(line);
		std::string time_str="", temp;
		float emcal=0., ihcal=0., ohcal=0.;
		LeakageDataPoint point;

		// Parse CSV: time,emcal imeas,ihcal imeas,ohcal imeas
		while (std::getline(ss, temp, ','))
		{
			if(n==24) std::cout<<temp<<std::endl;
			if(time_str.find("-")==std::string::npos) time_str = temp;
			else if (emcal == 0) emcal = std::stof(temp);
			else if (ihcal == 0) ihcal = std::stof(temp);
			else if (ohcal == 0) ohcal = std::stof(temp);
		}
			point.timestamp	  = time_str;
			point.imeas_emcal = emcal;
			point.imeas_ihcal = ihcal;
			point.imeas_ohcal = ohcal;
			
			data.push_back(point);
		
	}
	
	infile.close();
	return data;
}

int PlotLeakageCurrents(const std::string& csv_filename = "Leakage_Currents_per_12.csv")
{
	SetsPhenixStyle();
	Skaydis_colors* sc = new Skaydis_colors();
	gStyle->SetPalette(100, sc->Bi_gradient_PT);
	Int_t emcolor = sc->Bi_gradient_PT[0];
	Int_t ihcolor = sc->Bi_gradient_PT[49];
	Int_t ohcolor  = sc->Bi_gradient_PT[99];
	// Read data from CSV
	std::vector<LeakageDataPoint> data = ReadLeakageData(csv_filename);
	
	if (data.empty())
	{
		std::cerr << "Error: No data read from CSV file" << std::endl;
		return 1;
	}
	
	// Convert timestamps to hours from start
//	time_t start_time = data.back().timestamp; // Earliest time
	std::vector<std::string> time_hours;
	std::vector<double> emcal_imeas;
	std::vector<double> ihcal_imeas;
	std::vector<double> ohcal_imeas;
	
	for (const auto& point : data)
	{
		time_hours.push_back(point.timestamp);
		emcal_imeas.push_back(point.imeas_emcal);
		ihcal_imeas.push_back(point.imeas_ihcal);
		ohcal_imeas.push_back(point.imeas_ohcal);
	}
	int div_it= 1;
	int nbins = time_hours.size() / div_it;
	// Create graphs
	TH1F* g_emcal = new TH1F("hem", "EMCAL Leakage Current; Date; Leakage Current [mA]", nbins, 0, time_hours.size());
	TH1F* g_ihcal = new TH1F("hih", "IHCAL Leakage Current; Date; Leakage Current [mA]", nbins, 0, time_hours.size());
	TH1F* g_ohcal = new TH1F("hoh", "OHCAL Leakage Current; Date; Leakage Current [mA]", nbins, 0, time_hours.size());
	for(int i = 0; i<nbins; i++)
	{
		g_emcal->Fill(i+1, emcal_imeas.at(div_it*i));	
		if(i%50== 0) g_emcal->GetXaxis()->SetBinLabel(i+1, time_hours.at(div_it*i).c_str());
		g_ihcal->Fill(i+1, ihcal_imeas.at(div_it*i));	
		if(i%50== 0)g_ihcal->GetXaxis()->SetBinLabel(i+1, time_hours.at(div_it*i).c_str());
		g_ohcal->Fill(i+1, ohcal_imeas.at(div_it*i));	
		if(i%50== 0)g_ohcal->GetXaxis()->SetBinLabel(i+1, time_hours.at(div_it*i).c_str());
	}
	g_emcal->LabelsDeflate("X");
	g_ihcal->LabelsDeflate("X");
	g_ohcal->LabelsDeflate("X");
	// Style the graphs with sPHENIX colors
	// EMCAL - Enby Purple
	g_emcal->SetLineWidth(2);
	g_emcal->SetMarkerStyle(20);
	g_emcal->SetMarkerSize(2);
	g_emcal->SetLineColor(emcolor);
	g_emcal->SetMarkerColor(emcolor);
	
	// IHCAL - Black
	g_ihcal->SetLineWidth(2);
	g_ihcal->SetMarkerStyle(21);
	g_ihcal->SetMarkerSize(2);
	g_ihcal->SetLineColor(ihcolor);
	g_ihcal->SetMarkerColor(ihcolor);
	
	
	// OHCAL - Enby Yellow
	g_ohcal->SetLineWidth(2);
	g_ohcal->SetMarkerStyle(33);
	g_ohcal->SetMarkerSize(3);
	g_ohcal->SetLineColor(ohcolor);
	g_ohcal->SetMarkerColor(ohcolor);
	// Create multi-graph
//	TMultiGraph* mg = new TMultiGraph();
//	mg->Add(g_emcal);
//	mg->Add(g_ihcal);
//	mg->Add(g_ohcal);
	
	// Create canvas
	TCanvas* c = new TCanvas("c_leakage", "Calorimeter Leakage Currents", 1000, 600);
	c->SetLeftMargin(0.12);
	c->SetRightMargin(0.05);
	c->SetTopMargin(0.08);
	c->SetBottomMargin(0.12);
	
	// Draw multi-graph
//	mg->Draw("alp pmc plc plf");
	THStack* s=new THStack("s", "s");
	float hc_max = 4 * g_ohcal->GetMaximum();
	float scale = g_emcal->GetMaximum() / hc_max;
	std::cout<<scale<<std::endl;
	g_ihcal->Scale(scale);
	g_ohcal->Scale(scale);
	s->Add(g_emcal);
	s->Add(g_ihcal);
	s->Add(g_ohcal);
	//s->Draw("HIST P PMC PLC PLF");
	g_emcal->Draw("HIST P");
	g_ohcal->Draw("HIST P same");
	g_ihcal->SetLineColor(ihcolor);
	g_ihcal->SetMarkerColor(ihcolor);
	g_ihcal->Draw("same HIST P");
	float xmax = g_emcal->GetNbinsX()+1;
	float ymax = g_emcal->GetMaximum();
	TGaxis* axis = new TGaxis(xmax, 0, xmax, ymax, 0, hc_max, 510, "+L");
	axis->SetLineColor(ihcolor);
	axis->SetTitle("Leakage Currents [mA]--HCAL Scale");
	axis->Draw("same");
	// Style axes
	TAxis* xaxis = g_emcal->GetXaxis();
	TAxis* yaxis = g_emcal->GetYaxis();
	
	xaxis->SetTitle("Date");
	yaxis->SetTitle("Leakage Current [mA]--EMCAL");
//	yaxis->SetLineColor(g_emcal->GetLineColor());
	
	xaxis->SetTitleSize(0.045);
	yaxis->SetTitleSize(0.045);
	xaxis->SetLabelSize(0.04);
	yaxis->SetLabelSize(0.04);
	xaxis->SetTitleOffset(1.2);
	yaxis->SetTitleOffset(1.3);
	
	// Create legend
	TLegend* leg = new TLegend(0.15, 0.75, 0.32, 0.92);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->AddEntry(g_emcal, "EMCAL", "lp");
	leg->AddEntry(g_ihcal, "IHCAL", "lp");
	leg->AddEntry(g_ohcal, "OHCAL", "lp");
	leg->Draw("same");
	
	// Add sPHENIX label
	TPaveText* pt = new TPaveText(0.05, 0.92, 0.35, 0.98, "NB NDC");
	pt->SetFillColorAlpha(kWhite, 0.0);
	pt->SetBorderSize(0);
	pt->SetTextSize(0.04);
	pt->AddText("#it{#bf{sPHENIX}} Internal");
	pt->Draw("same");
	
	c->RedrawAxis();
	c->SaveAs("LeakageCurrents.pdf");
	c->SaveAs("LeakageCurrents.png");
	
	std::cout << "Plot saved as LeakageCurrents.pdf and LeakageCurrents.png" << std::endl;
	
	return 0;
}
#endif
