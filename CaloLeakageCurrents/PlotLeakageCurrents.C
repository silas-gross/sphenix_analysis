#pragma once 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include "sPHENIXStyle.C"
#include "sPHENIXStyle.h"
#include "Special_colors.h"
#include <TLegend.h>
#include <TCanvas.h>
#include <TH1.h>

#include <vector>
#include <string>
#include <format>

void SetHeaderLegend(TLegend* l1)
{

}
int PlotLeakageCurrents()
{
	SetsPhenixStyle();
	Skaydis_colors* sc = new Skaydis_colors();
	
	gStyle->SetPalette(100, sc->Enby_gradient_PT);
	return 0;
}
#endif
