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
#include <TF1.h>

#define Pi 3.1415926535

class DiscretizationCorrection
{	
	public:
		DiscretizationCorrection()
		{
			fillr2s();
		}			
		
		~DiscretizationCorrection(){};
		
		float getLatticeSize();
		float getJetRadius();
		bool getIsSizeCorrection() { return this->isSizeCorrection;}
		bool getIsDiscretizationCorrection(){ return this->isDiscretizationCorrection;}

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
		
		float first_order_cont_jet(float RL)
		{
			float ep = 2*PI*std::pow(m_jetradius, 2)*std::sqrt(RL)*( 
		}
		float all_order_continous_jet(float RL)
		{
			//This is the solution of the integral for geometric pairs at distance RL where both points lie within the jet
			float I_inc=R_L*std::pow(2*PI*std::pow(m_jetradius, 2), -1)*std::sqrt(4*std::pow(m_jetradius, 2)-std::pow(RL,2)) * all_order_continous_full(RL);
			float elipint=std::ellint_1f(std::arccos(RL*std::pow(2*m_jetradius, -1)), std::pow(RL*std::pow(m_jetradius, -1), 2)) *all_order_continuous_full(RL);
			I_inc+=elipint;
			return I_inc;
		}
		float all_order_continous_full(float RL)
		{
			//this is the total pairs in the phase space defined by the points lying within the jet and their pairs at distance R_l irrespective of the pair being in the jet radius
			float I_Total=std::pow(2*PI, 2)*RL*m_jetradius;
			return I_Total;
		}
		float all_order_discrete_jet(float RL) 
		{
			//pairs with both in jet
			float c=0;
			int i=floor(std::pow((RL+m_jetradius)/m_latticesize, 2));
			for(int j=0; j<=i; j++) 
			{
				if(j < std::pow(m_jetradius/m_latticesize, 2)) continue;
				c+=this->r2s.at(j);
			}
			float ct=all_order_discrete_full(RL) - c;
			return ct;
		}
		float all_order_discrete_full(float RL)
		{
			//pairs with one in jet
			float c=0;
			int i=floor(std::pow(RL/m_latticesize, 2));
			int counter=0;
			for(auto n=this->r2s.begin(); n!=this->r2s.end(); ++n)
			{
				c+=*n;
				counter++;
				if(counter >=i) break;
			}
			c=c*r2s.at(i);
			return c;
		}
		int sum_of_squares(int r) 
		{
			int r2=0.;
			std::vector<int> divs=getDivisors(r);
			std::pair<int, int> d1_3=divisors_mod_1_3(divs);
			r2=4*(d1_3.first-d1_3.second);
			return r2;
		}

		std::pair<int,int> divisors_mod_1_3(std::vector<int> divisors)
		{
			std::pair<int, int> d1_3 {0,0};
			for(auto i:divisors)
			{
				if(i%4==1) d1_3.first++;
				else if (i%4==3) d1_3.second++;
				else continue;
			}
			return d1_3;
		}

		std::vector<int> getDivisors(int r)
		{
			std::vector<int> divisors; 
			for(int i=1; i<=r; i++)
			{
				if(r % i ==0 ) divisors.push_back(i);
				else continue;
			}
			return divisors;
		}
		void fillr2s()
		{
			float max=0.;
			if(isSizeCorrection) max=m_jetradius/m_latticesize;
			else max=std::sqrt(std::pow(2.2, 2) +std::pow(2*PI, 2))/m_latticesize; //set to sPHENIX phase space
			for(int i=0; i<std::pow(max, 2); i++)
			{
				this->r2s.push_back(sum_of_squares(r));
			}
			return;
		}
		std::vector<int> r2s;
		TF1 *Continous_jet_f1, *Continuous_full_f1; 
		TF1 *Discrete_jet_f1,  *Discete_full_f1;
		TF1 *Corretion_funct;
		
}
#endif
