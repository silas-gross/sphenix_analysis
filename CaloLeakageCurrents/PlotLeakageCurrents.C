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
	time_t timestamp;
	double imeas_emcal;
	double imeas_ihcal;
	double imeas_ohcal;
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
	
	while (std::getline(infile, line))
	{
		if (line.empty()) continue;
		
		std::stringstream ss(line);
		std::string time_str;
		double emcal, ihcal, ohcal;
		char comma;
		
		// Parse CSV: time,emcal imeas,ihcal imeas,ohcal imeas
		if (std::getline(ss, time_str, ',') && 
		    ss >> emcal >> comma && 
		    ss >> ihcal >> comma && 
		    ss >> ohcal)
		{
			LeakageDataPoint point;
			
			// Parse ISO 8601 timestamp
			struct tm tm = {};
			std::stringstream ts(time_str);
			ts >> std::get_time(&tm, "%Y-%m-%d %H:%M:%S");
			point.timestamp = mktime(&tm);
			
			point.imeas_emcal = emcal;
			point.imeas_ihcal = ihcal;
			point.imeas_ohcal = ohcal;
			
			data.push_back(point);
		}
	}
	
	infile.close();
	return data;
}

int PlotLeakageCurrents(const std::string& csv_filename = "leakage_currents.csv")
{
	SetsPhenixStyle();
	Skaydis_colors* sc = new Skaydis_colors();
	gStyle->SetPalette(100, sc->Enby_gradient_PT);
	// Read data from CSV
	std::vector<LeakageDataPoint> data = ReadLeakageData(csv_filename);
	
	if (data.empty())
	{
		std::cerr << "Error: No data read from CSV file" << std::endl;
		return 1;
	}
	
	// Convert timestamps to hours from start
	time_t start_time = data.back().timestamp; // Earliest time
	std::vector<double> time_hours;
	std::vector<double> emcal_imeas;
	std::vector<double> ihcal_imeas;
	std::vector<double> ohcal_imeas;
	
	for (const auto& point : data)
	{
		double hours = static_cast<double>(point.timestamp - start_time) / 3600.0;
		time_hours.push_back(hours);
		emcal_imeas.push_back(point.imeas_emcal);
		ihcal_imeas.push_back(point.imeas_ihcal);
		ohcal_imeas.push_back(point.imeas_ohcal);
	}
	
	// Create graphs
	TGraph* g_emcal = new TGraph(time_hours.size(), time_hours.data(), emcal_imeas.data());
	TGraph* g_ihcal = new TGraph(time_hours.size(), time_hours.data(), ihcal_imeas.data());
	TGraph* g_ohcal = new TGraph(time_hours.size(), time_hours.data(), ohcal_imeas.data());
	// Style the graphs with sPHENIX colors
	// EMCAL - Enby Purple
	g_emcal->SetLineWidth(2);
	g_emcal->SetMarkerStyle(20);
	g_emcal->SetMarkerSize(0.8);
	
	// IHCAL - Black
	g_ihcal->SetLineWidth(2);
	g_ihcal->SetMarkerStyle(21);
	g_ihcal->SetMarkerSize(0.8);
	
	
	// OHCAL - Enby Yellow
	g_ohcal->SetLineWidth(2);
	g_ohcal->SetMarkerStyle(22);
	g_ohcal->SetMarkerSize(0.8);
	
	// Create multi-graph
	TMultiGraph* mg = new TMultiGraph();
	mg->Add(g_emcal);
	mg->Add(g_ihcal);
	mg->Add(g_ohcal);
	
	// Create canvas
	TCanvas* c = new TCanvas("c_leakage", "Calorimeter Leakage Currents", 1000, 600);
	c->SetLeftMargin(0.12);
	c->SetRightMargin(0.05);
	c->SetTopMargin(0.08);
	c->SetBottomMargin(0.12);
	
	// Draw multi-graph
	mg->Draw("ALP  PMC PLC PLF");
	
	// Style axes
	TAxis* xaxis = mg->GetXaxis();
	TAxis* yaxis = mg->GetYaxis();
	
	xaxis->SetTitle("Time (hours from start)");
	yaxis->SetTitle("Leakage Current [mA]");
	
	xaxis->SetTitleSize(0.045);
	yaxis->SetTitleSize(0.045);
	xaxis->SetLabelSize(0.04);
	yaxis->SetLabelSize(0.04);
	xaxis->SetTitleOffset(1.2);
	yaxis->SetTitleOffset(1.3);
	
	// Create legend
	TLegend* leg = new TLegend(0.65, 0.75, 0.92, 0.92);
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
