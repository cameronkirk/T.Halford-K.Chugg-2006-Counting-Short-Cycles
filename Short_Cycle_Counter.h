/* Short_Cycle_Counter.h

   Copyright (c) 2005 Thomas R. Halford 
   All rights reserved.
 
   Developed by: Thomas R. Halford
                 Communication Sciences Institute
                 University of Southern California
                 http://csi.usc.edu
 
   Permission is hereby granted, free of charge, to any person obtaining a copy 
   of this software and associated documentation files (the "Software"), to deal 
   with the Software without restriction, including without limitation the rights 
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
   copies of the Software, and to permit persons to whom the Software is furnished 
   to do so, subject to the following conditions:
  
    * Redistributions of source code must retain the above copyright notice, this 
      list of conditions and the following disclaimers.
    * Redistributions in binary form must reproduce the above copyright notice, 
      this list of conditions and the following disclaimers in the documentation 
      and/or other materials provided with the distribution.
    * Neither the names of Thomas R. Halford, the Communication Sciences Institute, 
      the University of Southern California nor the names of its contributors may 
      be used to endorse or promote products derived from this Software without 
      specific prior written permission. 

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
   INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
   PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR 
   COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
   AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION 
   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
*/

// Short_Cycle_Counter.h implements the short cycle counting
// algorithm described in Halford and Chugg's paper
// "An Algorithm for Counting Short Cycles in Bipartite Graphs".

#ifndef SHORT_CYCLE_COUNTER
#define SHORT_CYCLE_COUNTER

#include "Short_Cycle_Matrix.h"

typedef Short_Cycle_Matrix SCM;

class Short_Cycle_Counter
{
  public:
	// Constructors.
	Short_Cycle_Counter( void );
	Short_Cycle_Counter( const Short_Cycle_Matrix& E );
	~Short_Cycle_Counter( void );
	void initialize( const Short_Cycle_Matrix& E );

	// Accessors for girth and number cycles.
	int girth( void ) { return g_;   };
	int Ng( void )    { return Ng_;  };
	int Ng2( void )   { return Ng2_; };
	int Ng4( void )   { return Ng4_; };
	
	// Compute the mean and standard deviation of the
	// cycle 
	void cycle_dist( double* mu_g,  double* sdev_g, 
					 double* mu_g2, double* sdev_g2,
					 double* mu_g4, double* sdev_g4 );
	
	// Count cycles of length g, g+2 and g+4.  This method
	// determines g.
	void count( void );
	
  private:
	// Cycle counting helpers.  
	void process_P_U_2( void );
	void process_P_W_2( void );

	int  count_four_cycles( void );
	void count_six_eight_cycles( void );
	void count_ten_cycles_g_6( void );
	void count_twelve_cycles_g_8( void );
	void count_longer_cycles( void );
	
	int U_;				// |\mathcal{U}|
	int W_;				// |\mathcal{W}|
	int g_;				// girth
	int Ng_;			// N_g
	int Ng2_;			// N_{g+2}
	int Ng4_;			// N_{g+4}
	
	// Vectors that store the number of cycles of length
	// g, g+2 and g+4 incident on each vertex in U.
	double* Ng_per_u_;
	double* Ng2_per_u_;
	double* Ng4_per_u_;
	
	// Matrices required by the cycle counter.  Note that
	// the memory required by each matrix is allocated only when
	// the matrix is first required and freed as soon as possible
	// in order to minimize the memory footprint.
	
	SCM L_U_temp_, L_W_temp_;
	SCM E_, ET_;
	
	// Matrices used to count short cycles (4,6,8).
	SCM P_U_2_, P_W_2_;				 // P_2^\mathcal{U,W}
	SCM P_U_2_c2_, P_W_2_c2_;		 // \binom{P_2^\mathcal{U,W}}{2}
	SCM P_U_3_, P_W_3_;				 // P_3^\mathcal{U,W}
	SCM P_U_4_, P_W_4_;				 // P_4^\mathcal{U,W}
	SCM P_U_5_, P_W_5_;				 // P_5^\mathcal{U,W}
	SCM P_U_6_, P_W_6_;				 // P_6^\mathcal{U,W}
	SCM P_U_7_, P_W_7_;				 // P_7^\mathcal{U,W}
	SCM L_U_0_2_m1_, L_W_0_2_m1_;    // \max{L_{(0,2)}^\mathcal{U,W}-1,0}
	SCM	L_U_0_2_m2_, L_W_0_2_m2_;    // \max{L_{(0,2)}^\mathcal{U,W}-2,0}
	SCM L_U_1_2_, L_W_1_2_;			 // L_{(1,2)}^\mathcal{U,W}
	SCM L_U_0_4_, L_W_0_4_;			 // L_{(0,4)}^\mathcal{U,W}
	SCM L_U_2_2_, L_W_2_2_;			 // L_{(2,2)}^\mathcal{U,W}
	SCM L_U_1_4_, L_W_1_4_;			 // L_{(1,4)}^\mathcal{U,W}
	SCM L_U_3_2_, L_W_3_2_;			 // L_{(3,2)}^\mathcal{U,W}
	SCM L_U_0_6_, L_W_0_6_;			 // L_{(0,6)}^\mathcal{U,W}
	SCM L_U_2_4_, L_W_2_4_;			 // L_{(2,4)}^\mathcal{U,W}
	SCM L_U_4_2_, L_W_4_2_;			 // L_{(4,2)}^\mathcal{U,W}
	SCM L_U_1_6_, L_W_1_6_;			 // L_{(1,6)}^\mathcal{U,W}
	SCM L_U_3_4_;					 // L_{(3,4)}^\mathcal{U}
	SCM L_U_5_2_, L_W_5_2_;			 // L_{(5,2)}^\mathcal{U,W}
	SCM L_U_0_8_, L_W_0_8_;			 // L_{(0,8)}^\mathcal{U,W}
	
	SCM Big_Term_A_, Big_Term_B_, Big_Term_C_, Big_Term_D_;
	
	// Long cycle matrices.
	SCM P_U_gm3_, P_W_gm3_;		 // P_{g-3}^\mathcal{U,W}
	SCM P_U_gm2_, P_W_gm2_;		 // P_{g-2}^\mathcal{U,W}
	SCM P_U_gm1_, P_W_gm1_;      // P_{g-1}^\mathcal{U,W}
	SCM P_U_g_, P_W_g_;			 // P_g^\mathcal{U,W}
	SCM P_U_g1_, P_W_g1_;		 // P_{g+1}^\mathcal{U,W}
	SCM P_U_g2_, P_W_g2_;		 // P_{g+2}^\mathcal{U,W}
	SCM P_U_g3_, P_W_g3_;		 // P_{g+3}^\mathcal{U,W}
	SCM L_U_0_g_, L_W_0_g_;		 // L_{(0,g)}^\mathcal{U,W}
	SCM L_U_0_g2_, L_W_0_g2_;    // L_{(0,g+2)}^\mathcal{U,W}
	SCM L_U_0_g4_, L_W_0_g4_;    // L_{(0,g+4)}^\mathcal{U,W}
	SCM L_U_1_g_, L_W_1_g_;      // L_{(1,g)}^\mathcal{U,W}
	SCM L_U_1_g2_, L_W_1_g2_;    // L_{(1,g+2)}^\mathcal{U,W}
	SCM L_U_2_g_, L_W_2_g_;      // L_{(2,g)}^\mathcal{U,W}
	SCM L_U_3_g_, L_W_3_g_;		 // L_{(3,g)}^\mathcal{U,W}
	SCM L_U_gm6_2_, L_W_gm6_2_;  // L_{(g-6,2)}^\mathcal{U,W}
	SCM L_U_gm5_2_, L_W_gm5_2_;  // L_{(g-5,2)}^\mathcal{U,W}
	SCM L_U_gm4_2_, L_W_gm4_2_;  // L_{(g-4,2)}^\mathcal{U,W}
	SCM L_U_gm3_2_, L_W_gm3_2_;  // L_{(g-3,2)}^\mathcal{U,W}
	SCM L_U_gm2_2_, L_W_gm2_2_;  // L_{(g-2,2)}^\mathcal{U,W}
	SCM L_U_gm1_2_, L_W_gm1_2_;  // L_{(g-1,2)}^\mathcal{U,W}
	SCM L_U_g_2_, L_W_g_2_;      // L_{(g,2)}^\mathcal{U,W}
	SCM L_U_g1_2_, L_W_g1_2_;    // L_{(g+1,2)}^\mathcal{U,W}
};

#endif

