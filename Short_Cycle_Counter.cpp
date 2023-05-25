/* Short_Cycle_Counter.cpp

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

#include <math.h>
// NOTE: Location of cblas.h is machine dependent.
//#include "openblas/cblas.h"
#include "vecLib/cblas.h"
#include "Short_Cycle_Counter.h"

// Constructors.
Short_Cycle_Counter::Short_Cycle_Counter( void ) 
 : U_(0), W_(0), g_(4), Ng_(0), Ng2_(0), Ng4_(0),  
   Ng_per_u_(NULL), Ng2_per_u_(NULL), Ng4_per_u_(NULL)
{
	return;
} 

Short_Cycle_Counter::Short_Cycle_Counter( const Short_Cycle_Matrix& E )
{
	initialize(E);
}	
	
void Short_Cycle_Counter::initialize( const Short_Cycle_Matrix& E )
{
	U_   = E.e_nr();
	W_   = E.e_nc();
	g_   = 1000000;
	Ng_  = 0;
	Ng2_ = 0;
	Ng4_ = 0;
	Ng_per_u_  = new double[U_];
	memset(Ng_per_u_,0,U_*sizeof(double));
	Ng2_per_u_ = new double[U_];
	memset(Ng_per_u_,0,U_*sizeof(double));
	Ng4_per_u_ = new double[U_];
	memset(Ng_per_u_,0,U_*sizeof(double));
	E_   = E;
	ET_.transpose(E_);
}

Short_Cycle_Counter::~Short_Cycle_Counter( void )
{
	if( Ng_per_u_ )  delete [] Ng_per_u_;
	if( Ng2_per_u_ ) delete [] Ng2_per_u_;
	if( Ng4_per_u_ ) delete [] Ng4_per_u_;
}

// The main counting method.  Different counting helper
// methods are called depending on the girth.
void Short_Cycle_Counter::count( void )
{
	g_   = 1000000;
	Ng_  = 0;
	Ng2_ = 0;
	Ng4_ = 0;
	// Count 4 cycles first to determine girth.
	if( count_four_cycles() )
	{
		// girth == 4
		count_six_eight_cycles();
	}
	
	else
	{
		// Count 6 and 8 cycles to determine the girth.
		count_six_eight_cycles();
		if( g_ == 6 )
		{
			Ng_  = L_U_0_6_.int_trace()/6;
			L_U_0_6_.diagonal(Ng_per_u_);
			
			Ng2_ = L_U_0_8_.int_trace()/8;
			L_U_0_8_.diagonal(Ng2_per_u_);
		
			// Count the 10 cycles with girth = 6.
			count_ten_cycles_g_6();			
			Ng4_ = L_U_0_g4_.int_trace()/10;
			L_U_0_g4_.diagonal(Ng4_per_u_);
		}
		
		else if( g_ == 8 )
		{
			Ng_ = L_U_0_8_.int_trace()/8;
			L_U_0_8_.diagonal(Ng_per_u_);
			// Count 10 and 12 cycles when the girth = 8.
			count_ten_cycles_g_6(); 
			count_twelve_cycles_g_8();
			Ng2_ = L_U_0_g2_.int_trace()/10;
			L_U_0_g2_.diagonal(Ng2_per_u_);
			Ng4_ = L_U_0_g4_.int_trace()/12;
			L_U_0_g4_.diagonal(Ng4_per_u_);
		}
		
		else count_longer_cycles();
	}

}

// Returns 1 if the girth is 4, 0 otherwise.
int Short_Cycle_Counter::count_four_cycles( void )
{
	// Compute P_U_2, P_W_2, L_U_0_2_m1, L_W_0_2_m1, 
	// L_U_0_2_m2, L_W_0_2_m2, P_U_2_c2, and P_W_2_c2
	// simulataneously for speed.
	process_P_U_2();
	process_P_W_2();
	
	// Compute L_U_1_2 and L_W_1_2.
	L_U_1_2_.matrix_mult(E_,L_W_0_2_m1_);
	L_W_1_2_.matrix_mult(ET_,L_U_0_2_m1_);

	// Compute P_U_3, P_W_3.
	P_U_3_.matrix_mult(P_U_2_,E_);
	P_U_3_ -= L_U_1_2_;
	P_W_3_.transpose(P_U_3_);
															
	// Compute L_U_0_4, L_W_0_4.
	L_U_0_4_.mx_mult_diag(P_U_3_,ET_);
	L_W_0_4_.mx_mult_diag(P_W_3_,E_);

	int temp = L_U_0_4_.int_trace();
	if( temp == 0 ) return 0;
	else
	{
		g_  = 4;
		Ng_ = temp/4;
		L_U_0_4_.diagonal(Ng_per_u_);
		return 1;
	}
}

void Short_Cycle_Counter::count_six_eight_cycles( void )
{
	// Compute L_U_2_2, L_W_2_2.
	L_U_2_2_.mx_mult_zero(E_,L_W_1_2_);
	L_W_2_2_.mx_mult_zero(ET_,L_U_1_2_);
	
	// Comput P_U_4, P_W_4.
	P_U_4_.matrix_mult(P_U_3_,ET_); 
	P_U_4_ -= L_U_0_4_; 
	P_U_4_ -= L_U_2_2_;
	
	P_W_4_.matrix_mult(P_W_3_,E_);  
	P_W_4_ -= L_W_0_4_; 
	P_W_4_ -= L_W_2_2_;
	
	// Compute L_U_1_4, L_W_1_4.
	L_U_1_4_.matrix_mult(E_,L_W_0_4_);  
	L_U_temp_ = P_U_3_*E_;
	L_U_1_4_ -= L_U_temp_;
	L_U_1_4_ -= L_U_temp_;
	
	L_W_1_4_.matrix_mult(ET_,L_U_0_4_); 
	L_W_temp_ = P_W_3_*ET_;
	L_W_1_4_ -= L_W_temp_;
	L_W_1_4_ -= L_W_temp_;
	
	// Compute L_U_3_2, L_W_3_2.
	L_U_3_2_.matrix_mult(P_U_3_,L_W_0_2_m1_); 
	L_U_3_2_ -= L_U_temp_;
	
	L_W_3_2_.matrix_mult(P_W_3_,L_U_0_2_m1_); 
	L_W_3_2_ -= L_W_temp_;
	
	// Compute P_U_5, P_W_5.  
	P_U_5_.matrix_mult(P_U_4_,E_);  
	P_U_5_ -= L_U_1_4_; 
	P_U_5_ -= L_U_3_2_; 
	
	P_W_5_.transpose(P_U_5_); 
	
	// Compute L_U_0_6, L_W_0_6.
	L_U_0_6_.mx_mult_diag(P_U_5_,ET_);
	L_W_0_6_.mx_mult_diag(P_W_5_,E_);
	
	int temp = L_U_0_6_.int_trace();
	if( g_ == 4 ) 
	{
		Ng2_ = temp/6; 
		L_U_0_6_.diagonal(Ng2_per_u_);
	}
	
	else if( temp > 0 ) g_ = 6;
	
	// Compute L_U_2_4, L_W_2_4.
	L_U_2_4_.mx_mult_zero(E_,L_W_1_4_);  
	if( g_ == 4 ) 
	{
		L_U_temp_ = P_U_2_.mx_choose_3(6.0);
		L_U_2_4_ -= L_U_temp_;
	}
	
	L_W_1_4_.delete_data();
	
	L_W_2_4_.mx_mult_zero(ET_,L_U_1_4_); 
	if( g_ == 4 ) 
	{
		L_W_temp_ = P_W_2_.mx_choose_3(6.0);
		L_W_2_4_ -= L_W_temp_;
	}

	// Compute L_U_4_2, L_W_4_2.
	L_U_4_2_.mx_mult_zero(E_,L_W_3_2_);  
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_2_2_); 
	L_U_4_2_ -= L_U_temp_; 
	if( g_ == 4 )
	{
		L_U_4_2_ += P_U_2_c2_;
		L_U_4_2_ += P_U_2_c2_;
	}
	
	L_U_2_2_.delete_data();
	
	L_W_4_2_.mx_mult_zero(ET_,L_U_3_2_); 
	L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_2_2_); 
	L_W_4_2_ -= L_W_temp_; 
	if( g_ == 4 )
	{
		L_W_4_2_ += P_W_2_c2_;
		L_W_4_2_ += P_W_2_c2_;
	}
	
	L_W_2_2_.delete_data();
	
	// Compute P_U_6, P_W_6.
	P_U_6_.matrix_mult(P_U_5_,ET_); 
	P_U_6_ -= L_U_0_6_; 
	P_U_6_ -= L_U_2_4_; 
	P_U_6_ -= L_U_4_2_;
	L_U_2_4_.delete_data();
	
	P_W_6_.matrix_mult(P_W_5_,E_); 
	P_W_6_ -= L_W_0_6_; 
	P_W_6_ -= L_W_2_4_; 
	P_W_6_ -= L_W_4_2_;
		
	// Compute L_U_1_6, L_W_1_6.
	L_U_1_6_.matrix_mult(E_,L_W_0_6_); 
	L_U_temp_  = P_U_5_;
	L_U_temp_ *= E_;
	L_U_temp_ *= 2.0;
	L_U_1_6_  -= L_U_temp_;
	if( g_ == 4 )
	{
		L_U_temp_  = P_U_3_.mx_choose_2(2.0);
		L_U_temp_ *= E_;
		L_U_1_6_  -= L_U_temp_;
		L_U_temp_.matrix_mult(P_U_2_c2_,E_);
		L_U_temp_ -= P_U_3_;
		L_U_temp_ *= E_;
		L_U_temp_ *= 2.0;
		L_U_1_6_  += L_U_temp_;
		L_U_temp_.matrix_mult(E_,P_W_2_c2_);
		L_U_temp_ -= P_U_3_;
		L_U_temp_ *= E_;
		L_U_temp_ *= 2.0;
		L_U_1_6_  += L_U_temp_;
	}

	L_W_1_6_.matrix_mult(ET_,L_U_0_6_); 
	L_W_temp_  = P_W_5_;
	L_W_temp_ *= ET_;
	L_W_temp_ *= 2.0;
	L_W_1_6_  -= L_W_temp_;
	if( g_ == 4 )
	{
		L_W_temp_  = P_W_3_.mx_choose_2(2.0);
		L_W_temp_ *= ET_;
		L_W_1_6_  -= L_W_temp_;
		L_W_temp_.matrix_mult(P_W_2_c2_,ET_);
		L_W_temp_ -= P_W_3_;
		L_W_temp_ *= ET_;
		L_W_temp_ *= 2.0;
		L_W_1_6_  += L_W_temp_;
		L_W_temp_.matrix_mult(ET_,P_U_2_c2_);
		L_W_temp_ -= P_W_3_;
		L_W_temp_ *= ET_;
		L_W_temp_ *= 2.0;
		L_W_1_6_  += L_W_temp_;
	}

	// Compute L_U_3_4.
	if( g_ == 4 )
	{
		L_U_3_4_.matrix_mult(E_,L_W_2_4_); 
		L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_1_4_); 
		L_U_3_4_  -= L_U_temp_;
		L_U_temp_  = P_U_3_.mx_choose_2(4.0);
		L_U_temp_ *= E_;
		L_U_3_4_  -= L_U_temp_;
		L_U_temp_.matrix_mult(P_U_2_c2_,E_);
		L_U_temp_ -= P_U_3_;
		L_U_temp_ *= E_;
		L_U_temp_ *= 4.0;
		L_U_3_4_  += L_U_temp_;
		L_U_temp_.matrix_mult(E_,P_W_2_c2_);
		L_U_temp_ -= P_U_3_;
		L_U_temp_ *= E_;
		L_U_temp_ *= 6.0;
		L_U_3_4_  += L_U_temp_;
	}
	
	else{ L_U_3_4_ = L_U_1_4_; } // Both 0 if g_ > 4.
	
	L_U_1_4_.delete_data();

	// Compute L_U_5_2 and L_W_5_2.
	L_U_5_2_.matrix_mult(E_,L_W_4_2_);
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_3_2_);
	L_U_5_2_ -= L_U_temp_;
	L_U_temp_ = P_U_5_*E_;
	L_U_5_2_ -= L_U_temp_;
	if( g_ == 4 )
	{
		L_U_temp_  = P_U_3_*E_;
		L_U_5_2_  += L_U_temp_;
		L_U_5_2_  += L_U_temp_; 
		L_U_temp_.matrix_mult(L_U_0_4_,L_U_1_2_);
		L_U_5_2_  -= L_U_temp_;
		L_U_temp_  = L_U_3_2_*E_;
		L_U_5_2_  += L_U_temp_;
		L_W_temp_.matrix_mult(P_U_3_,L_W_0_2_m2_);
		L_U_temp_  = L_W_temp_*E_;
		L_U_5_2_  += L_U_temp_;
		L_U_5_2_  += L_U_temp_;
		L_W_temp_.matrix_mult(P_U_2_c2_,E_);
		L_W_temp_ -= P_U_3_;
		L_U_temp_  = L_W_temp_*E_;
		L_U_5_2_  += L_U_temp_;
		L_U_5_2_  += L_U_temp_;
	}

	P_U_2_c2_.delete_data(); 

	L_W_5_2_.matrix_mult(ET_,L_U_4_2_);
	L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_3_2_);
	L_W_5_2_ -= L_W_temp_;
	L_W_temp_ = P_W_5_*ET_;
	L_W_5_2_ -= L_W_temp_;
	if( g_ == 4 )
	{
		L_W_temp_.matrix_mult(L_W_0_4_,L_W_1_2_);
		L_W_5_2_  -= L_W_temp_;
		L_W_temp_  = L_W_3_2_;
		L_U_temp_.matrix_mult(P_W_3_,L_U_0_2_m2_);
		L_W_temp_ += L_U_temp_;
		L_W_temp_ += L_U_temp_;
		L_U_temp_.matrix_mult(P_W_2_,ET_);
		L_W_temp_ += L_U_temp_;
		L_W_temp_ += L_U_temp_;
		L_W_5_2_  += L_W_temp_;	
	}	
	
	P_W_2_c2_.delete_data();
	
	// Compute P_U_7 and P_W_7.
	P_U_7_.matrix_mult(P_U_6_,E_); 
	P_U_7_ -= L_U_1_6_; 
	P_U_7_ -= L_U_3_4_; 
	P_U_7_ -= L_U_5_2_;	
	L_U_3_4_.delete_data();
	
	P_W_7_.transpose(P_U_7_);
	
	// Compute L_U_0_8 and L_W_0_8.
	L_U_0_8_.mx_mult_diag(P_U_7_,ET_);
	L_W_0_8_.mx_mult_diag(P_W_7_,E_);
	
	temp = L_U_0_8_.int_trace();
	if( g_ == 4 ) 
	{
		Ng4_ = temp/8;
		L_U_0_8_.diagonal(Ng4_per_u_);
	}
	
	else if( temp && g_ != 6 ) g_ = 8;	 
}

void Short_Cycle_Counter::count_ten_cycles_g_6( void )
{	
	L_U_0_4_.delete_data(); 
	L_W_0_4_.delete_data();

	// Compute L_U_2_6, L_W_2_6.
	L_U_2_g_.mx_mult_zero(E_,L_W_1_6_);
	L_W_2_g_.mx_mult_zero(ET_,L_U_1_6_);

	// Compute L_U_6_2, L_W_6_2.
	L_U_g_2_.mx_mult_zero(E_,L_W_5_2_);  
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_4_2_); 
	L_U_g_2_ -= L_U_temp_;
	L_U_temp_ = P_U_4_*P_U_2_; 
	L_U_g_2_ += L_U_temp_;
	L_U_4_2_.delete_data(); 
	P_U_4_.delete_data();
	
	L_W_g_2_.mx_mult_zero(ET_,L_U_5_2_);  
	L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_4_2_); 
	L_W_g_2_ -= L_W_temp_;
	L_W_temp_ = P_W_4_*P_W_2_; 
	L_W_g_2_ += L_W_temp_;
	L_W_4_2_.delete_data();  
	P_W_4_.delete_data(); 
	
	// Compute P_U_8.
	P_U_g2_.matrix_mult(P_U_7_,ET_); 
	P_U_g2_ -= L_U_0_8_; 
	P_U_g2_ -= L_U_2_g_; 
	P_U_g2_ -= L_U_g_2_;

	// Compute L_U_1_8.
	L_U_1_g2_.matrix_mult(E_,L_W_0_8_);  
	L_U_temp_ = 2.0*P_U_7_; L_U_temp_ *= E_;  
	L_U_1_g2_ -= L_U_temp_;
	
	// Compute L_U_3_6.
	L_U_3_g_.matrix_mult(E_,L_W_2_g_); 
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_1_6_); 
	L_U_3_g_ -= L_U_temp_; 
	L_U_temp_ = P_U_3_.mx_choose_3(6.0);
	L_U_3_g_ -= L_U_temp_;
	P_U_3_.delete_data();
	L_U_1_6_.delete_data(); 

	// Compute L_U_7_2.
	L_U_g1_2_.matrix_mult(E_,L_W_g_2_); 
	L_U_temp_ = P_U_7_*E_; 
	L_U_g1_2_ -= L_U_temp_; 
	L_U_temp_ = L_U_5_2_*E_; 
	L_U_g1_2_ += L_U_temp_;
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_5_2_); 
	L_U_g1_2_ -= L_U_temp_; 
	L_U_temp_.matrix_mult(L_U_0_6_,L_U_1_2_); 
	L_U_g1_2_ -= L_U_temp_;
	L_U_temp_ = 2.0*P_U_5_; 
	L_U_temp_ *= E_; 
	L_U_g1_2_ += L_U_temp_; 
	L_U_temp_.matrix_mult(P_U_5_,L_W_0_2_m2_);
	L_U_temp_ *= 2.0; 
	L_U_temp_ *= E_;
	L_U_g1_2_ += L_U_temp_;
	L_U_5_2_.delete_data();
	P_U_5_.delete_data();

	// Compute P_U_9.
	P_U_g3_.matrix_mult(P_U_g2_,E_); 
	P_U_g3_ -= L_U_1_g2_; 
	P_U_g3_ -= L_U_3_g_; 
	P_U_g3_ -= L_U_g1_2_;

	// Compute L_U_0_10.
	L_U_0_g4_.mx_mult_diag(P_U_g3_,ET_);
}

// Count 12 cycles when the girth is known to be 8.
void Short_Cycle_Counter::count_twelve_cycles_g_8( void )
{
	L_U_0_6_.delete_data();
	// Compute the matrices that weren't computed in count_ten_cycles_g_6().
	P_W_7_.transpose(P_U_7_);

	P_W_g2_.matrix_mult(P_W_7_,E_); 
	P_W_g2_ -= L_W_0_8_; 
	P_W_g2_ -= L_W_2_g_; 
	P_W_g2_ -= L_W_g_2_;
	L_W_1_g2_.matrix_mult(ET_,L_U_0_8_);
	//THIS LINE CAUSES ERROR, SPLITTING INTO TWO LINES FIXES IT, INVESTIGATE LATER?
	//L_W_temp_ = 2.0*P_W_7_; 
	L_W_temp_ = P_W_7_;
	L_W_temp_ *= 2.0; 
	L_W_temp_ *= ET_;
	L_W_1_g2_ -= L_W_temp_;
	L_W_g1_2_.matrix_mult(ET_,L_U_g_2_); 
	L_W_temp_ = P_W_7_*ET_; 
	L_W_g1_2_ -= L_W_temp_; 
	L_W_temp_ = L_W_5_2_*ET_; 
	L_W_g1_2_ += L_W_temp_;
	L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_5_2_); 
	L_W_g1_2_ -= L_W_temp_; 
	L_W_temp_.matrix_mult(L_W_0_6_,L_W_1_2_); 
	L_W_g1_2_ -= L_W_temp_;
	L_W_temp_ = 2.0*P_W_5_; 
	L_W_temp_ *= ET_; 
	L_W_g1_2_ += L_W_temp_; 
	L_W_temp_.matrix_mult(P_W_5_,L_U_0_2_m2_);
	L_W_temp_ *= 2.0; 
	L_W_temp_ *= ET_; 
	L_W_g1_2_ += L_W_temp_;
	L_W_0_6_.delete_data();
	P_W_7_.delete_data();
	L_W_1_2_.delete_data();
	P_W_5_.delete_data(); 
	L_U_0_2_m2_.delete_data();
	L_W_3_g_.matrix_mult(ET_,L_U_2_g_); 
	L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_1_6_);
	L_W_3_g_ -= L_W_temp_;
	L_W_temp_ = P_W_3_.mx_choose_3(6.0); 
	L_W_3_g_ -= L_W_temp_;
	P_W_g3_.matrix_mult(P_W_g2_,ET_); 
	P_W_g3_ -= L_W_1_g2_; 
	P_W_g3_ -= L_W_3_g_; 
	P_W_g3_ -= L_W_g1_2_;
	L_W_0_g4_.mx_mult_diag(P_W_g3_,E_);
	P_W_g2_.delete_data();
	L_W_1_6_.delete_data();
	// Copy the needed matrices from count_ten_cycles_g_6() results.
	P_U_g1_    = P_U_g3_;
	P_W_g1_    = P_W_g3_;
	L_U_0_g2_  = L_U_0_g4_;
	L_W_0_g2_  = L_W_0_g4_;
	L_U_gm2_2_ = L_U_g_2_;
	L_W_gm2_2_ = L_W_g_2_; 
	L_U_gm1_2_ = L_U_g1_2_;
	L_W_gm1_2_ = L_W_g1_2_;
	L_U_1_g_   = L_U_1_g2_;
	L_W_1_g_   = L_W_1_g2_;
	
	P_W_g3_.delete_data(); 
	L_W_g1_2_.delete_data(); 
	L_W_1_g2_.delete_data(); 
	L_W_0_g4_.delete_data();

	// Compute L_U_2_8, L_W_2_8.
	L_U_2_g_.mx_mult_zero(E_,L_W_1_g_);
	
	L_W_2_g_.mx_mult_zero(ET_,L_U_1_g_);
	L_W_1_g_.delete_data();

	// Compute L_U_8_2, L_W_8_2.
	L_U_g_2_.mx_mult_zero(E_,L_W_gm1_2_);  
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_gm2_2_); 
	L_U_g_2_ -= L_U_temp_;
	L_U_temp_ = P_U_6_*P_U_2_; 
	L_U_g_2_ += L_U_temp_;
	L_W_gm1_2_.delete_data();
	L_U_gm2_2_.delete_data();
	P_U_2_.delete_data();
	P_U_6_.delete_data();
	
	L_W_g_2_.mx_mult_zero(ET_,L_U_gm1_2_);  
	L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_gm2_2_); 
	L_W_g_2_ -= L_W_temp_;
	L_W_temp_ = P_W_6_*P_W_2_; 
	L_W_g_2_ += L_W_temp_;
	L_W_0_2_m1_.delete_data();  
	L_W_gm2_2_.delete_data();
	P_W_2_.delete_data();  
	P_W_6_.delete_data(); 
	
	// Compute P_U_10.
	P_U_g2_.matrix_mult(P_U_g1_,ET_); 
	P_U_g2_ -= L_U_0_g2_; 
	P_U_g2_ -= L_U_2_g_; 
	P_U_g2_ -= L_U_g_2_;
	L_U_2_g_.delete_data(); 
	L_U_g_2_.delete_data(); 
	
	// Compute L_U_1_10.
	L_U_1_g2_.matrix_mult(E_,L_W_0_g2_);  
	L_U_temp_ = 2.0*P_U_g1_; 
	L_U_temp_ *= E_;  
	L_U_1_g2_ -= L_U_temp_;
	L_W_0_g2_.delete_data(); 

	// Compute L_U_3_8.
	L_U_3_g_.matrix_mult(E_,L_W_2_g_);  
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_1_g_); 
	L_U_3_g_ -= L_U_temp_; 
	L_W_2_g_.delete_data(); 
	L_U_1_g_.delete_data();

	// Compute L_U_9_2.
	L_U_g1_2_.matrix_mult(E_,L_W_g_2_); 
	L_U_temp_ = P_U_g1_*E_; 
	L_U_g1_2_ -= L_U_temp_; 
	L_U_temp_ = L_U_gm1_2_*E_; 
	L_U_g1_2_ += L_U_temp_;
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_gm1_2_); 
	L_U_g1_2_ -= L_U_temp_; 
	L_U_temp_.matrix_mult(L_U_0_8_,L_U_1_2_); 
	L_U_g1_2_ -= L_U_temp_;
	L_U_temp_ = 2.0*P_U_7_; 
	L_U_temp_ *= E_; 
	L_U_g1_2_ += L_U_temp_; 
	L_U_temp_.matrix_mult(P_U_7_,L_W_0_2_m2_);
	L_U_temp_ *= 2.0; 
	L_U_temp_ *= E_; 
	L_U_g1_2_ += L_U_temp_;
	L_W_g_2_.delete_data(); 
	P_U_g1_.delete_data(); 
	L_U_gm1_2_.delete_data(); 
	L_U_0_2_m1_.delete_data();
	L_U_1_2_.delete_data(); 
	P_U_7_.delete_data(); 
	L_W_0_2_m2_.delete_data(); 

	// Compute P_U_11.
	P_U_g3_.matrix_mult(P_U_g2_,E_); 
	P_U_g3_ -= L_U_1_g2_; 
	P_U_g3_ -= L_U_3_g_; 
	P_U_g3_ -= L_U_g1_2_;
	P_U_g2_.delete_data(); 
	L_U_1_g2_.delete_data(); 
	L_U_3_g_.delete_data(); 
	L_U_g1_2_.delete_data(); 

	// Compute L_U_0_12.
	L_U_0_g4_.mx_mult_diag(P_U_g3_,ET_);
	P_U_g3_.delete_data();
}

// Assumes that count_six_eight_cycles has been called.
void Short_Cycle_Counter::count_longer_cycles( void )
{
	int temp;
	
	// Set the maximum girth.
	int max_girth = 2*(U_ > W_ ? U_ : W_);
	g_ = max_girth+2;
	
	// Free as much memory as possible.  Memory is freed for all matrices
	// not longer needed.
	L_U_0_2_m2_.delete_data(); 
	L_W_1_2_.delete_data(); 
	L_U_0_4_.delete_data(); 
	L_W_0_4_.delete_data(); 
	P_U_4_.delete_data(); 
	P_W_4_.delete_data();
	P_U_5_.delete_data();
	P_W_5_.delete_data(); 
	L_U_0_6_.delete_data(); 
	L_W_0_6_.delete_data();
	P_U_6_.delete_data(); 
	P_W_6_.delete_data(); 
	L_U_1_6_.delete_data(); 
	L_W_1_6_.delete_data();
	L_U_0_8_.delete_data(); 
	L_W_0_8_.delete_data();
	
	// Starting with the assumption that g_ = 10, initialize the
	// search for the girth.
	
	// Compute L_U_6_2, L_W_6_2.  L_U_4_2, L_W_4_2 not needed after this.
	L_U_gm4_2_.matrix_mult(E_,L_W_5_2_);  
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_4_2_); 
	L_U_gm4_2_ -= L_U_temp_;
	L_U_4_2_.delete_data(); 
	
	L_W_gm4_2_.matrix_mult(ET_,L_U_5_2_); 
	L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_4_2_); 
	L_W_gm4_2_ -= L_W_temp_;
	L_W_4_2_.delete_data();

	// Compute P_U_8, P_W_8.  P_U_7, P_W_7 not needed after this.
	P_U_gm2_.matrix_mult(P_U_7_,ET_); 
	P_U_gm2_ -= L_U_gm4_2_;
	P_U_7_.delete_data();
	
	P_W_gm2_.matrix_mult(P_W_7_,E_);  
	P_W_gm2_ -= L_W_gm4_2_;
	P_W_7_.delete_data();

	// Compute L_U_7_2, L_W_7_2.  L_U_5_2, L_W_5_2 not needed after this.
	L_U_gm3_2_.matrix_mult(E_,L_W_gm4_2_);  
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_5_2_); 
	L_U_gm3_2_ -= L_U_temp_;
	L_U_5_2_.delete_data();
	
	L_W_gm3_2_.matrix_mult(ET_,L_U_gm4_2_); 
	L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_5_2_); 
	L_W_gm3_2_ -= L_W_temp_;
	L_W_5_2_.delete_data();
	
	// Compute P_U_9, P_W_9.
	P_U_gm1_.matrix_mult(P_U_gm2_,E_);  
	P_U_gm1_ -= L_U_gm3_2_;
	
	P_W_gm1_.matrix_mult(P_W_gm2_,ET_); 
	P_W_gm1_ -= L_W_gm3_2_;
	
	// Search for the girth.
	for( int gtry = 10; gtry <= max_girth; gtry += 2 )
	{
		// Assume g_ = grty and calculate L_U_0_g_, L_W_0_g_.
		L_U_0_g_.mx_mult_diag(P_U_gm1_,ET_);
		L_W_0_g_.mx_mult_diag(P_W_gm1_,E_);
		
		if( temp = L_U_0_g_.int_trace() )
		{
			// Cycles of length gtry exist.
			g_  = gtry;
			Ng_ = temp/g_;
			L_U_0_g_.diagonal(Ng_per_u_);

			// The following no longer needed.
			P_U_gm3_.delete_data(); 
			P_W_gm3_.delete_data(); 
			L_U_gm6_2_.delete_data();
			L_W_gm6_2_.delete_data(); 
			L_U_gm5_2_.delete_data(); 
			L_W_gm5_2_.delete_data();
			
			break; 
		}
		
		else
		{
			// Prepare for next recursion by updating P_U_gm1_, P_W_gm1_, P_U_gm2_, P_W_gm2_.
			L_U_gm6_2_ = L_U_gm4_2_;
			L_W_gm6_2_ = L_W_gm4_2_;
			L_U_gm5_2_ = L_U_gm3_2_;
			L_W_gm5_2_ = L_W_gm3_2_;
			P_U_gm3_   = P_U_gm1_;
			P_W_gm3_   = P_W_gm1_;
		
			// Compute L_U_gm4_2, L_W_gm4_2.
			L_U_gm4_2_.matrix_mult(E_,L_W_gm5_2_);  
			L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_gm6_2_); 
			L_U_gm4_2_ -= L_U_temp_;
		
			L_W_gm4_2_.matrix_mult(ET_,L_U_gm5_2_); 
			L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_gm6_2_); 
			L_W_gm4_2_ -= L_W_temp_;
		
			// Compute P_U_gm2, P_W_gm2.
			P_U_gm2_.matrix_mult(P_U_gm3_,ET_); 
			P_U_gm2_ -= L_U_gm4_2_;
		
			P_W_gm2_.matrix_mult(P_W_gm3_,E_);  
			P_W_gm2_ -= L_W_gm4_2_;

			// Compute L_U_gm3_2, L_W_gm3_2.
			L_U_gm3_2_.matrix_mult(E_,L_W_gm4_2_);  
			L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_gm5_2_); 
			L_U_gm3_2_ -= L_U_temp_;
			
			L_W_gm3_2_.matrix_mult(ET_,L_U_gm4_2_); 
			L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_gm5_2_); 
			L_W_gm3_2_ -= L_W_temp_;
	
			// Compute P_U_gm1, P_W_gm1.
			P_U_gm1_.matrix_mult(P_U_gm2_,E_);  
			P_U_gm1_ -= L_U_gm3_2_;
			
			P_W_gm1_.matrix_mult(P_W_gm2_,ET_); 
			P_W_gm1_ -= L_W_gm3_2_;		
		}
	}
	
	// Check if g > max_girth.  If so, this must be
	// a tree so set g_ to 1000000 and exit.
	if( g_ == max_girth+2 )
	{
		g_ = 1000000;
		Ng_ = Ng2_ = Ng4_ = 0;
		return;
	}

	// If g = max_girth then there can be no cycles
	// of length g+2 or g+4.
	else if( g_ == max_girth )
	{
		Ng2_ = Ng4_ = 0;
		return;
	}
	
	// Compute Ng2 via L_U_0_g2.
	
	// Compute L_U_gm2_2, L_W_gm2_2.  L_U_gm4_2, L_W_gm4_2 no longer needed.
	L_U_gm2_2_.matrix_mult(E_,L_W_gm3_2_);  
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_gm4_2_); 
	L_U_gm2_2_ -= L_U_temp_; 
	L_U_gm4_2_.delete_data();
	
	L_W_gm2_2_.matrix_mult(ET_,L_U_gm3_2_); 
	L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_gm4_2_); 
	L_W_gm2_2_ -= L_W_temp_; 
	L_W_gm4_2_.delete_data();
	
	// Compute P_U_g, P_W_g.
	P_U_g_.matrix_mult(P_U_gm1_,ET_); 
	P_U_g_ -= L_U_0_g_; 
	P_U_g_ -= L_U_gm2_2_;
	
	P_W_g_.matrix_mult(P_W_gm1_,E_);  
	P_W_g_ -= L_W_0_g_; 
	P_W_g_ -= L_W_gm2_2_;
	
	// Compute L_U_1_g, L_W_1_g.  L_W_0_g_ is no longer needed.
	L_U_1_g_.matrix_mult(E_,L_W_0_g_);  
	L_U_temp_ = 2.0*P_U_gm1_; 
	L_U_temp_ *= E_; 
	L_U_1_g_ -= L_U_temp_;
	
	L_W_1_g_.matrix_mult(ET_,L_U_0_g_); 
	L_W_temp_ = 2.0*P_W_gm1_; 
	L_W_temp_ *= ET_; 
	L_W_1_g_ -= L_W_temp_;
	L_W_0_g_.delete_data();
	
	// Compute L_U_gm1_2, L_W_gm1_2.  L_U_gm3_2, L_W_gm3_2, P_W_gm1
	// are no longer needed.
	L_U_gm1_2_.matrix_mult(E_,L_W_gm2_2_); 
	L_U_temp_ = P_U_gm1_*E_; 
	L_U_gm1_2_ -= L_U_temp_;
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_gm3_2_);
	L_U_gm1_2_ -= L_U_temp_;
	L_U_gm3_2_.delete_data();
	
	L_W_gm1_2_.matrix_mult(ET_,L_U_gm2_2_); 
	L_W_temp_ = P_W_gm1_*ET_; 
	L_W_gm1_2_ -= L_W_temp_;
	L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_gm3_2_); 
	L_W_gm1_2_ -= L_W_temp_;
	P_W_gm1_.delete_data(); 
	L_W_gm3_2_.delete_data();

	// Compute P_U_g1, P_W_g1.  P_U_g, P_W_g no longer needed.
	P_U_g1_.matrix_mult(P_U_g_,E_);  
	P_U_g1_ -= L_U_1_g_; 
	P_U_g1_ -= L_U_gm1_2_;
	P_U_g_.delete_data(); 
	
	P_W_g1_.matrix_mult(P_W_g_,ET_); 
	P_W_g1_ -= L_W_1_g_; 
	P_W_g1_ -= L_W_gm1_2_;
	P_W_g_.delete_data();
	
	// Compute L_U_0_g2, L_W_0_g2.  P_W_g1 no longer needed.
	L_U_0_g2_.mx_mult_diag(P_U_g1_,ET_);
	L_W_0_g2_.mx_mult_diag(P_W_g1_,E_);
	P_W_g1_.delete_data();
	
	Ng2_ = L_U_0_g2_.int_trace()/(g_+2);
	L_U_0_g2_.diagonal(Ng2_per_u_);
	
	// If g = max_girth-2, there can exist no
	// cycles of length g+4 so exit.
	if( g_ == max_girth-2 )
	{
		Ng4_ = 0;
		return;
	}
	
	// Compute Ng4 via L_U_0_g4.  Free memory ASAP.
	
	// Compute L_U_2_g, L_W_2_g.
	L_U_2_g_.mx_mult_zero(E_,L_W_1_g_);
	L_W_2_g_.mx_mult_zero(ET_,L_U_1_g_);	
	L_W_1_g_.delete_data();
	
	// Compute L_U_g_2, L_W_g_2.
	L_U_g_2_.mx_mult_zero(E_,L_W_gm1_2_); 
	L_U_temp_ = P_U_gm2_*P_U_2_; 
	L_U_g_2_ += L_U_temp_;
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_gm2_2_); 
	L_U_g_2_ -= L_U_temp_;
	P_U_2_.delete_data(); 
	P_U_gm2_.delete_data();
	L_W_gm1_2_.delete_data(); 
	L_U_gm2_2_.delete_data();
	
	L_W_g_2_.mx_mult_zero(ET_,L_U_gm1_2_); 
	L_W_temp_ = P_W_gm2_*P_W_2_; 
	L_W_g_2_ += L_W_temp_;
	L_W_temp_.matrix_mult(L_W_0_2_m1_,L_W_gm2_2_); 
	L_W_g_2_ -= L_W_temp_;
	P_W_2_.delete_data(); 
	P_W_gm2_.delete_data(); 
	L_W_0_2_m1_.delete_data(); 
	L_W_gm2_2_.delete_data();

	// Compute P_U_g2.
	P_U_g2_.matrix_mult(P_U_g1_,ET_); 
	P_U_g2_ -= L_U_0_g2_; 
	P_U_g2_ -= L_U_2_g_; 
	P_U_g2_ -= L_U_g_2_;
	L_U_2_g_.delete_data();
	
	// Compute L_U_1_g2.
	L_U_1_g2_.matrix_mult(E_,L_W_0_g2_); 
	L_U_temp_ = 2.0*P_U_g1_; 
	L_U_temp_ *= E_; 
	L_U_1_g2_ -= L_U_temp_;
	L_W_0_g2_.delete_data(); 
	
	// Compute L_U_3_g.
	L_U_3_g_.matrix_mult(E_,L_W_2_g_); 
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_1_g_); 
	L_U_3_g_ -= L_U_temp_;
	L_W_2_g_.delete_data(); 
	L_U_1_g_.delete_data();
	
	// Compute L_U_g1_2.
	L_U_g1_2_.matrix_mult(E_,L_W_g_2_); 
	L_U_temp_ = P_U_g1_*E_; 
	L_U_g1_2_ -= L_U_temp_; 
	L_U_temp_ = L_U_gm1_2_*E_; 
	L_U_g1_2_ += L_U_temp_; 
	L_U_temp_.matrix_mult(L_U_0_2_m1_,L_U_gm1_2_); 
	L_U_g1_2_ -= L_U_temp_;
	L_U_temp_.matrix_mult(L_U_0_g_,L_U_1_2_); 
	L_U_g1_2_ -= L_U_temp_; 
	L_U_temp_ = 2.0*P_U_gm1_; 
	L_U_temp_ *= E_;
	L_U_g1_2_ += L_U_temp_; 
	L_U_temp_.matrix_mult(P_U_gm1_,L_W_0_2_m2_); 
	L_U_temp_ *= 2.0; 
	L_U_temp_ *= E_;
	L_U_g1_2_ += L_U_temp_;
	L_U_0_2_m1_.delete_data();
	L_W_0_2_m2_.delete_data();
	L_U_1_2_.delete_data();
	L_W_g_2_.delete_data(); 
	P_U_g1_.delete_data();
	L_U_gm1_2_.delete_data(); 
	P_U_gm1_.delete_data(); 
	
	// Compute P_U_g3.
	P_U_g3_.matrix_mult(P_U_g2_,E_); 
	P_U_g3_ -= L_U_1_g2_; 
	P_U_g3_ -= L_U_3_g_; 
	P_U_g3_ -= L_U_g1_2_;
	P_U_g2_.delete_data(); 
	L_U_1_g2_.delete_data(); 
	L_U_3_g_.delete_data(); 
	L_U_g1_2_.delete_data();
	
	// Compute L_U_0_g4.
	L_U_0_g4_.mx_mult_diag(P_U_g3_,ET_);
	P_U_g3_.delete_data();
	
	Ng4_ = L_U_0_g4_.int_trace()/(g_+4);
	L_U_0_g4_.diagonal(Ng4_per_u_);
}

void Short_Cycle_Counter::process_P_U_2( void )
{
	P_U_2_.matrix_mult(E_,ET_);
	
	P_U_2_c2_.copy_size(P_U_2_);
	L_U_0_2_m1_.copy_size(P_U_2_);
	L_U_0_2_m2_.copy_size(P_U_2_);
	
	P_U_2_c2_.reset_data();
	L_U_0_2_m1_.reset_data();
	L_U_0_2_m2_.reset_data();	
				
	for( int rr = 0, oo = 0; rr < P_U_2_.i_nr(); rr++, oo += P_U_2_.i_nc() )
	{
		for( int cc = 0, pp = oo; cc < P_U_2_.i_nc(); cc++, pp++ )
		{
			double val = P_U_2_[pp];
			if( val > 1 ) P_U_2_c2_.set_el(pp,val*(val-1.0)/2.0);

			if( rr == cc )
			{
				P_U_2_.set_el(pp,0); P_U_2_c2_.set_el(pp,0);
				if( val > 1 )
				{
					L_U_0_2_m1_.set_el(pp,val-1.0);
					if( val > 2 ) L_U_0_2_m2_.set_el(pp,val-2.0);
				} 				
			}
		}
	}
}

void Short_Cycle_Counter::process_P_W_2( void )
{
	P_W_2_.matrix_mult(ET_,E_);
	
	P_W_2_c2_.copy_size(P_W_2_);
	L_W_0_2_m1_.copy_size(P_W_2_);
	L_W_0_2_m2_.copy_size(P_W_2_);
	
	P_W_2_c2_.reset_data();
	L_W_0_2_m1_.reset_data();
	L_W_0_2_m2_.reset_data();	
				
	for( int rr = 0, oo = 0; rr < P_W_2_.i_nr(); rr++, oo += P_W_2_.i_nc() )
	{
		for( int cc = 0, pp = oo; cc < P_W_2_.i_nc(); cc++, pp++ )
		{
			double val = P_W_2_[pp];
			if( val > 1 )
			{
				P_W_2_c2_.set_el(pp,val*(val-1.0)/2.0);
			}

			if( rr == cc )
			{
				P_W_2_.set_el(pp,0); P_W_2_c2_.set_el(pp,0);
				if( val > 1 )
				{
					L_W_0_2_m1_.set_el(pp,val-1.0);
					if( val > 2 ) L_W_0_2_m2_.set_el(pp,val-2.0);
				} 				
			}
		}
	}
}

void Short_Cycle_Counter::cycle_dist( double* mu_g,  double* sdev_g, 
									  double* mu_g2, double* sdev_g2,
									  double* mu_g4, double* sdev_g4 )
{
	// Compute the means.
	double t_g = 0.0, t_g2 = 0.0, t_g4 = 0.0;
	for( int ii = 0; ii < U_; ii++ )
	{
		t_g  += Ng_per_u_[ii];
		t_g2 += Ng2_per_u_[ii];
		t_g4 += Ng4_per_u_[ii];
	}
	
	*mu_g  = t_g/(U_+0.0);
	*mu_g2 = t_g2/(U_+0.0);
	*mu_g4 = t_g4/(U_+0.0);
	
	// Compute the standard deviations.
	t_g = t_g2 = t_g4 = *sdev_g = *sdev_g2 = *sdev_g4 = 0.0;
	for( int jj = 0; jj < U_; jj++ )
	{
		t_g  += Ng_per_u_[jj]-*mu_g;
		t_g2 += Ng2_per_u_[jj]-*mu_g2;
		t_g4 += Ng4_per_u_[jj]-*mu_g4;

		*sdev_g  += (Ng_per_u_[jj]-*mu_g)*(Ng_per_u_[jj]-*mu_g); 
		*sdev_g2 += (Ng2_per_u_[jj]-*mu_g2)*(Ng2_per_u_[jj]-*mu_g2);
		*sdev_g4 += (Ng4_per_u_[jj]-*mu_g4)*(Ng4_per_u_[jj]-*mu_g4);
	}

	*sdev_g  = sqrt((*sdev_g-t_g*t_g/(U_+0.0))/(U_-1.0));
	*sdev_g2 = sqrt((*sdev_g2-t_g2*t_g2/(U_+0.0))/(U_-1.0));
	*sdev_g4 = sqrt((*sdev_g4-t_g4*t_g4/(U_+0.0))/(U_-1.0));
}