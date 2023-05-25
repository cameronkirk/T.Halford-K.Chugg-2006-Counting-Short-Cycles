/* main.cpp

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

// main.cpp provides a sample usage of the Short_Cycle_Counter and 
// Short_Cycle_Matrix classes.  

#include <iostream>
#include <fstream>
#include "Short_Cycle_Matrix.h"
#include "Short_Cycle_Counter.h"
#include <time.h>
#include <windows.h>
using namespace std;

int main( int argc, const char* argv[] ) 
{
    if( (argc != 2) && (argc != 4) ) 
	{
		cout << "TWO USAGES: " << argv[0] << " nc nr inicidence_matrix_filename" << endl
		     << "            " << argv[0] << " alist_filename" << endl;
	}
																			
	else
	{
		Short_Cycle_Matrix E;
		if( argc == 4 ) E.read_incidence_matrix_file(atoi(argv[1]),atoi(argv[2]),argv[3]);
		else		    E.read_alist_file(argv[1]);
		Short_Cycle_Counter E_counter(E);

		LARGE_INTEGER freq, t1, t2;
		double elapsedTime;
		QueryPerformanceFrequency(&freq);
		QueryPerformanceCounter(&t1);

		E_counter.count();
		int g = E_counter.girth();
		std::cout << endl
		     << "Cycle Count:" << endl
			 << "girth = " << g << endl
		     << "N_" << g   << " = " << E_counter.Ng()  << endl
			 << "N_" << g+2 << " = " << E_counter.Ng2() << endl
			 << "N_" << g+4 << " = " << E_counter.Ng4() << endl;
		
		double mg = 0.0, sg = 0.0, mg2 = 0.0, sg2 = 0.0, mg4 = 0.0, sg4 = 0.0;
		E_counter.cycle_dist(&mg,&sg,&mg2,&sg2,&mg4,&sg4);
		cout << endl	 
			 << "Cycle Distribution:" << endl
			 << "mu_" << g   << " = " << mg  << ", sigma_" << g   << " = " << sg  << endl
			 << "mu_" << g+2 << " = " << mg2 << ", sigma_" << g+2 << " = " << sg2 << endl
			 << "mu_" << g+4 << " = " << mg4 << ", sigma_" << g+4 << " = " << sg4 << endl
			 << endl;

		QueryPerformanceCounter(&t2);
		elapsedTime = (t2.QuadPart - t1.QuadPart) * 1.0 / freq.QuadPart;
		std::cout << elapsedTime << "\n";
	}
}
