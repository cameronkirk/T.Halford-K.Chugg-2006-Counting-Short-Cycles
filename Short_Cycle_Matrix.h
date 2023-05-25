/* Short_Cycle_Matrix.h

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

// Short_Cycle_Matrix.h defines the matrix class used by the
// short cycle counter.  The distinction between "internal" and "external"
// exists in order to take advantage of Apple's AltiVec instruction set
// that requires matrix dimensions to be multiples of 4.
 
#ifndef SHORT_CYCLE_MATRIX
#define SHORT_CYCLE_MATRIX

#include <iostream>
 
class Short_Cycle_Matrix
{
  public:
	// Constructors.
	Short_Cycle_Matrix( void );
	Short_Cycle_Matrix( const Short_Cycle_Matrix& copy_mx );
	~Short_Cycle_Matrix( void );
	Short_Cycle_Matrix& operator=( const Short_Cycle_Matrix& copy_mx ); 
	
	// Methods for reading matrices from files.
	void read_incidence_matrix_file( int nc, int nr, const char* filename );
	void read_alist_file( const char* filename );
	
	// Accessors.
	int set( void )  const { return set_; };
	int i_nc( void ) const { return i_nc_; };
	int i_nr( void ) const { return i_nr_; };
	int e_nc( void ) const { return e_nc_; };
	int e_nr( void ) const { return e_nr_; };
	double get_el( int p ) { return data_[p]; };
	double get_el( int r, int c ) const { return data_[i_nc_*r+c]; };
	double operator[]( int p ) const { return data_[p]; };
	double operator()( int r, int c ) const { return data_[i_nc_*r+c]; };
	
	// Data accessor needed by BLAS multiplication methods.
	double* data( void ) { return data_; }; 
	
	// Dimension setting methods.
	void set_e_nc( int e_nc ) { e_nc_ = e_nc; };
	void set_e_nr( int e_nr ) { e_nr_ = e_nr; };
	void set_i_nc( void ); // Set i_nc_ using e_nc_.
	void set_i_nr( void ); // Set i_nr_ using e_nr_.
	void set_i_nc( int i_nc ) { i_nc_ = i_nc; };
	void set_i_nr( int i_nr ) { i_nr_ = i_nr; };
	void copy_size( const Short_Cycle_Matrix& copy_mx );
	void copy_transpose_size( const Short_Cycle_Matrix& copy_mx );
	
	// Elementing setting methods.
	void set_el( int r, int c, double v ) { data_[r*i_nc_+c] = v; };
	void dec_el( int p, double v ) { data_[p] -= v; };
	void set_el( int p, double v ) { data_[p] = v; };
	void reset_data( void );
	
	// Free the matrix memory.
	void delete_data( void );
	
	// Matrix trace.
	double trace( void ) const;
	int int_trace( void ) const { return (int)trace(); };
	
	// Matrix operations.
	void transpose( const Short_Cycle_Matrix& source );   // *this = source^T 
	void matrix_mult( Short_Cycle_Matrix& left,
					  Short_Cycle_Matrix& right );        // *this = left \times right
	void mx_mult_diag( const Short_Cycle_Matrix& left,
					   const Short_Cycle_Matrix& right ); // *this = (left \times right) o I
	void mx_mult_zero( Short_Cycle_Matrix& left,
					   Short_Cycle_Matrix& right );       // *this = Z(left \times right)
	void operator+=( Short_Cycle_Matrix& right );         // *this = *this + right
	void operator-=( Short_Cycle_Matrix& right );		  // *this = *this - right
	void operator*=( const Short_Cycle_Matrix& right );   // *this = *this o right
	void operator*=( double right );					  // *this = right*(*this)
	friend Short_Cycle_Matrix operator*( const Short_Cycle_Matrix& left,
										 const Short_Cycle_Matrix& right );
	friend Short_Cycle_Matrix operator*( double left, const Short_Cycle_Matrix& right );
	
	Short_Cycle_Matrix mx_choose_2( double mult_fac );
	Short_Cycle_Matrix mx_choose_3( double mult_fac );
	
	// Place the matrix diagonal in d.  The external matrix dimension is used. 
	void diagonal( double* d );
	
	// Output matrix for debugging.
	friend std::ostream& operator<<(std::ostream& os, Short_Cycle_Matrix& mx );
	
  private:
	double*     data_;
	int			i_nc_;		// "Internal" matrix dimensions.
	int			i_nr_;		
	int			e_nc_;      // "External" matrix dimensions. 
	int			e_nr_;
	int			set_;		// Flag indicating if memory allocated.
};
 
#endif



