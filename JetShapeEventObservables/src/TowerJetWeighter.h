#ifndef __TOWERJETWEIGHTER_H__
#define __TOWERJETWEIGHTER_H__

#include <TH1.h>
#include <TTree.h>

#include <array>
#include <vector>
#include <map>
#include <string>
#include <thread>

struct WeightedTower
{
	float r {-1};
	float phi {-100}; 
	float eta {-100};
	float weight {0};
	float ET_orig {-1};
	std::map<float, float> ET_cut {};
	float pt_min {-1};

};
class ShapeTrim
{
	public:
		ShapeTrim(){};
		~ShapeTrim(){};

	private:
		float cutET(float pt_min=1.);//this function runs the skiming to send ET to ET_iR 
};	
#endif
