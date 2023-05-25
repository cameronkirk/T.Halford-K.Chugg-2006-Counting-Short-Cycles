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

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
// NOTE: Location of cblas.h is machine dependent.
//#include "openblas/cblas.h"
#include "vecLib/cblas.h"
#include "Short_Cycle_Matrix.h"

using namespace std;
// Constructors.
Short_Cycle_Matrix::Short_Cycle_Matrix( void )
  : data_(NULL), i_nc_(0), i_nr_(0), e_nc_(0), e_nr_(0), set_(0)
{
	return;
}
	
Short_Cycle_Matrix::Short_Cycle_Matrix( const Short_Cycle_Matrix& copy_mx )
{
	*this = copy_mx;
}

Short_Cycle_Matrix::~Short_Cycle_Matrix( void )
{
	if( set_ ) 
	{ 
		delete [] data_; 
		set_ = 0; 
	}
}
	
Short_Cycle_Matrix& Short_Cycle_Matrix::operator=( const Short_Cycle_Matrix& copy_mx )
{
	if( set_ ) delete [] data_;
	
	// Copy the matrix dimensions.
	e_nc_ = copy_mx.e_nc_;
	e_nr_ = copy_mx.e_nr_;
	i_nc_ = copy_mx.i_nc_;
	i_nr_ = copy_mx.i_nr_;
	
	// Copy the matrix data.
	data_ = new double[i_nr_*i_nc_];
	memcpy(data_,copy_mx.data_,i_nr_*i_nc_*sizeof(double));
	set_ = 1;
	
	return* this;
}


// Read an incidence matrix from a file.
void Short_Cycle_Matrix::read_incidence_matrix_file( int nc, int nr, const char* filename )
{
	// Set the matrix dimensions.
	e_nc_ = nc;
	e_nr_ = nr;	
	set_i_nc();
	set_i_nr();

	// Allocate the data memory.
	reset_data();
	
	// Read the matrix data from the file.
	ifstream fin(filename);
	int tmp;
	for( int rr = 0, ww = 0; rr < e_nr_; rr++, ww += i_nc_ )
	{
		for( int cc = 0, pp = ww; cc < e_nc_; cc++, pp++ )
		{
			fin >> tmp; 
			data_[pp] = tmp+0.0;
		}
	}
	
	fin.close();
}

// Read an incidence matrix from an alist file.
void Short_Cycle_Matrix::read_alist_file( const char* filename )
{
	// Read the matrix dimensions from the file along
	// with the extraneous maximum vertex degree information. 
	ifstream fin(filename);
	int nc, nr, tmp, tmp2;
	fin >> nr >> nc >> tmp >> tmp2;
	cout << " " << nr << " " << nc << " " << tmp << " " << tmp2 << "\n";
	// Set the matrix dimensions.
	e_nc_ = nc;
	e_nr_ = nr;
	set_i_nc();
	set_i_nr();
	
	// Allocate the data memory.
	reset_data();
	

	// Read in the number of 1's per row.
	// Skip over the number of 1's per column.
	int ii, jj, rp;
	int* row_weights = new int[nr];
	for( ii = 0; ii < nr; ii++ ) { fin >> tmp; row_weights[ii] = tmp; }
	for( ii = 0; ii < nc; ii++ ) { fin >> tmp; }
	
	// Read in the 1's in each row and set the matrix data.
	for( ii = 0, rp = 0; ii < nr; ii++, rp+=i_nc_ )
	{
		for( jj = 0; jj < row_weights[ii]; jj++ )
		{
			fin >> tmp;
			//NEW LINES TO HANDLE D.MACKAY IRREGULAR GRAPHS (WHERE ALIST FILE ZERO PADS EACH LINE)
			while (tmp == 0)
				fin >> tmp;
			data_[rp+tmp-1] = 1.0; // The alist format indexes rows from 1.
			
		}
	}
	
	delete [] row_weights;
	fin.close();
}

// Given that the external nc has been set,
// set the appropriate internal nc.
void Short_Cycle_Matrix::set_i_nc( void )
{
	if( e_nc_%4 == 0 )      i_nc_ = e_nc_;
	else if( e_nc_%4 == 1 ) i_nc_ = e_nc_+3;
	else if( e_nc_%4 == 2 ) i_nc_ = e_nc_+2;
	else 					i_nc_ = e_nc_+1;
}

// Given that the external nr has been set,
// set the appropriate internal nr.
void Short_Cycle_Matrix::set_i_nr( void )
{
	if( e_nr_%4 == 0 )      i_nr_ = e_nr_;
	else if( e_nr_%4 == 1 ) i_nr_ = e_nr_+3;
	else if( e_nr_%4 == 2 ) i_nr_ = e_nr_+2;
	else 					i_nr_ = e_nr_+1;
}

// Reset the data to all zero elements.
void Short_Cycle_Matrix::reset_data( void )
{
	if( set_ ) delete [] data_;
	data_ = new double[i_nr_*i_nc_];
	memset(data_,0,i_nr_*i_nc_*sizeof(double));
	set_ = 1;
}

// Free the matrix data.
void Short_Cycle_Matrix::delete_data( void )
{
	if( set_ ) delete [] data_;
	set_ = 0;
}

// Copy the dimensions, but not data, of copy_mx into *this.
void Short_Cycle_Matrix::copy_size( const Short_Cycle_Matrix& copy_mx )
{
	e_nc_ = copy_mx.e_nc();
	e_nr_ = copy_mx.e_nr();
	i_nc_ = copy_mx.i_nc();
	i_nr_ = copy_mx.i_nr();
}

// Copy the dimensions, but not data, of copy_mx^T into *this. 
void Short_Cycle_Matrix::copy_transpose_size( const Short_Cycle_Matrix& copy_mx )
{
	e_nc_ = copy_mx.e_nr();
	e_nr_ = copy_mx.e_nc();
	i_nc_ = copy_mx.i_nr();
	i_nr_ = copy_mx.i_nc();
}

// Matrix trace.
double Short_Cycle_Matrix::trace( void ) const
{
	double sum = 0.0;
	for( int ii = 0, pp = 0; ii < e_nr_; ii++, pp+=(i_nc_+1) ) sum += data_[pp];
	return sum;
}

// Set *this to source^T.
void Short_Cycle_Matrix::transpose( const Short_Cycle_Matrix& source )
{
	if( set_ ) delete [] data_;
	copy_transpose_size(source);
	
	data_ = new double[i_nr_*i_nc_];
	for( int rr = 0, oo = 0; rr < i_nr_; rr++, oo += i_nc_ )
	{
		for( int cc = 0, pp = oo; cc < i_nc_; cc++, pp++ ) data_[pp] = source(cc,rr);
	}
	
	set_ = 1;
}

// *this = left \times right
void Short_Cycle_Matrix::matrix_mult( Short_Cycle_Matrix& left,
									  Short_Cycle_Matrix& right )
{
	// Set dimensions and allocate data memory.
	e_nc_ = right.e_nc();
	e_nr_ = left.e_nr();
	i_nc_ = right.i_nc();
	i_nr_ = left.i_nr();

	reset_data();
	
	// Use BLAS matrix multiplication.
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,i_nr_,i_nc_,left.i_nc(),1.0,
				left.data(),left.i_nc(),right.data(),i_nc_,0.0,data_,i_nc_);
}

// *this = (left \times right) o I
// Assumes a square matrix results.	
void Short_Cycle_Matrix::mx_mult_diag( const Short_Cycle_Matrix& left,
									   const Short_Cycle_Matrix& right )
{
	// Set dimensions and allocate memory.
	e_nc_ = right.e_nc();
	e_nr_ = left.e_nr();
	i_nc_ = right.i_nc();
	i_nr_ = left.i_nr();

	reset_data();
	
	// Perform multiplication only for diagonal elements.
	double sum;
	int inner = left.i_nc();
	for( int rr = 0; rr < i_nr_; rr++ )
	{
		sum = 0.0;
		for( int kk = 0; kk < inner; kk++ ) sum += left(rr,kk)*right(kk,rr);
		data_[rr*i_nc_+rr] = sum;
	}
}

// *this = Z[left \times right]
// Assumes a square matrix results.
void Short_Cycle_Matrix::mx_mult_zero( Short_Cycle_Matrix& left,
									   Short_Cycle_Matrix& right ) 
{
	matrix_mult(left,right);
	for( int r = 0; r < i_nr_; r++ ) data_[r*i_nc_+r] = 0.0;	
}

// Use BLAS matrix addition.
void Short_Cycle_Matrix::operator+=( Short_Cycle_Matrix& right )
{
	cblas_daxpy(i_nr_*i_nc_,1.0,right.data(),1,data_,1);
}

// Use BLAS matrix addition.
void Short_Cycle_Matrix::operator-=( Short_Cycle_Matrix& right )
{
	cblas_daxpy(i_nr_*i_nc_,-1.0,right.data(),1,data_,1);
}
	
// Direct matrix product.	
void Short_Cycle_Matrix::operator*=( const Short_Cycle_Matrix& right )
{		
	for( int rr = 0, oo = 0; rr < i_nr_; rr++, oo += i_nc_ )
	{
		for( int cc = 0, pp = oo; cc < i_nc_; cc++, pp++ )
		{
			data_[pp] *= right[pp];
		}
	}
}

// Multiplication by a constant.
void Short_Cycle_Matrix::operator*=( double right )
{
	for( int rr = 0, oo = 0; rr < i_nr_; rr++, oo += i_nc_ )
	{
		for( int cc = 0, pp = oo; cc < i_nc_; cc++, pp++ )
		{
			data_[pp] *= right;
		}
	}
}

// Direct matrix product.
Short_Cycle_Matrix operator*( const Short_Cycle_Matrix& left,
							  const Short_Cycle_Matrix& right )
{
	Short_Cycle_Matrix out;
	out.copy_size(left);
	out.reset_data();
	
	for( int rr = 0, oo = 0; rr < out.i_nr(); rr++, oo += out.i_nc() )
	{
		for( int cc = 0, pp = oo; cc < out.i_nc(); cc++, pp++ )
		{
			out.set_el(pp,left[pp]*right[pp]);
		}
	}
	
	return out;
}

// Multiplication by a constant.
Short_Cycle_Matrix operator*( double left, const Short_Cycle_Matrix& right )
{
	Short_Cycle_Matrix out = right;
	out *= left;
	return out;
}

Short_Cycle_Matrix Short_Cycle_Matrix::mx_choose_2( double mult_fac )
{
	Short_Cycle_Matrix out = *this;
	for( int ii = 0; ii < out.i_nc()*out.i_nr(); ii++ )
	{
		out.set_el(ii,mult_fac*out[ii]*(out[ii]-1.0)/2.0);
	}
	
	return out;
}

Short_Cycle_Matrix Short_Cycle_Matrix::mx_choose_3( double mult_fac )
{
	Short_Cycle_Matrix out = *this;
	for( int ii = 0; ii < out.i_nc()*out.i_nr(); ii++ )
	{
		out.set_el(ii,mult_fac*out[ii]*(out[ii]-1.0)*(out[ii]-2.0)/6.0);
	}
	
	return out;
}

void Short_Cycle_Matrix::diagonal( double* d )
{
	for( int ii = 0, pp = 0; ii < e_nr_; ii++, pp+=(i_nc_+1) ) d[ii] = data_[pp];
}

ostream& operator<<( ostream& os, Short_Cycle_Matrix& mx )
{
	for( int rr = 0, oo = 0; rr < mx.e_nr_; rr++, oo += mx.i_nc_ )
	{
		for( int cc = 0, pp = oo; cc < mx.e_nc_; cc++, pp++ )
		{
			os << mx.data_[pp] << " ";
		}
		
		os << endl;
	}

	return os;
}