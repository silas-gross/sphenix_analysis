#ifndef __DISCRETIAZATIONCORRECTION_H__
#define __DISCRETIAZATIONCORRECTION_H__

//////////////////////////////////////////////////////////////////////////
//		Discretization Correction for two point correlators 	//
//									//
//		Skaydi 							//
//		27 May 2025 						//
//									//
//	discretization correction and jet size corrections		//
//	see attached documentation for derivation and discussion 	//
//	Useful to compare calo only to truth/track			//
//	Works for 2pt and integrated 3pt, No idea on full 3		//
//////////////////////////////////////////////////////////////////////////

#include <math.h>

#define Pi 3.1415926535

class DiscretizationCorrection
{	
	public:
		DiscretizationCorrection();
		~DiscretizationCorrection();
		
		float getLatticeSize();
		float getJetRadius();
		bool getIsSizeCorrection();
		bool getIsDiscretizationCorrection();

		void setLatticeSize(float binsize);
		void setJetRadius(float jetradius);
		void setIsSizeCorrection();
		void setIsDiscretizationCorrection();
	
		float CalcluateCorrectionFactor(float RL);
		
	private:
		
		float m_latticesize		=0.;
		float m_jetradius		=0.;
		bool isSizeCorrection		=false;
		bool isDiscretizationCorrection	=false;

		int correction_order = 3; //how many powers of R_L to expand the elliptic integrals to
		float first_order_cont(float RL)
		{
			float x = 4
}
#endif
