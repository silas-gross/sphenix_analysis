// Tell emacs that this is a C++ source
//  -*- C++ -*-.

#ifndef __SOFTDROPCOMP_H__
#define __SOFTDROPCOMP_H__

#include <fastjet/PseduoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/contrib/SoftDrop.hh>

#include <vector>
#include <map> 
#include <string>
#include <math.h>
#include <array>

class SoftDropComp
{
	public:
		SoftDropComp(){};
		~SoftDropComp(){};
	private:
		int nJets {0};
		
};
#endif


