//
//  Preconditioned conjugate gradient solver
//
//  Created by Robert Bridson, Ryoichi Ando and Nils Thuerey and Fluidchen Team
//

#ifndef RCMATRIX3_H
#define RCMATRIX3_H

#include <iterator>
#include <cassert>
#include <vector>
#include <fstream>
#include <cmath>
#include <functional>
#include "Communication.hpp"
#include "Utilities.hpp"

// index type 
#define int_index long long

// parallelization disabled

#define parallel_for(size)  { int thread_number = 0; int_index parallel_index=0; for( int_index parallel_index=0; parallel_index<(int_index)size; parallel_index++ ) {
#define parallel_end 	} thread_number=parallel_index=0; }

#define parallel_block 
#define do_parallel 
#define do_end
#define block_end
// note - "Int" instead of "N" here, the latter is size!
template<class Int, class T>
struct InstantBLAS {
	static inline Int offset(Int N, Int incX) { return ((incX) > 0 ?  0 : ((N) - 1) * (-(incX))); }
	static T cblas_ddot( const Int N, const T *X, const Int incX, const T *Y, const Int incY) {
		double r  = 0.0; // always use double precision internally here...
		Int    i;
		Int    ix = offset(N,incX);
		Int    iy = offset(N,incY);
		for (i = 0; i < N; i++) {
			r  += X[ix] * Y[iy];
			ix += incX;
			iy += incY;
		}
		return (T)r;
	}
	static void cblas_daxpy( const Int N, const T alpha, const T *X, const Int incX, T *Y, const Int incY) {
		Int i;
		if (N     <= 0  ) return;
		if (alpha == 0.0) return;
		if (incX == 1 && incY == 1) {
			const Int m = N % 4;
			for (i = 0; i < m; i++)
				Y[i] += alpha * X[i];
			for (i = m; i + 3 < N; i += 4) {
				Y[i    ] += alpha * X[i    ];
				Y[i + 1] += alpha * X[i + 1];
				Y[i + 2] += alpha * X[i + 2];
				Y[i + 3] += alpha * X[i + 3];
			}
		} else {
			Int ix = offset(N, incX);
			Int iy = offset(N, incY);
			for (i = 0; i < N; i++) {
				Y[iy] += alpha * X[ix];
				ix    += incX;
				iy    += incY;
			}
		}
	}
	// dot products ==============================================================
	static inline T dot(const std::vector<T> &x, const std::vector<T> &y) {
		return cblas_ddot((int)x.size(), &x[0], 1, &y[0], 1);
	}
	
	// inf-norm (maximum absolute value: index of max returned) ==================
	static inline Int index_abs_max(const std::vector<T> &x) {
		int maxind = 0;
		T maxvalue = 0;
		for(Int i = 0; i < (Int)x.size(); ++i) {
			if(std::abs(x[i]) > maxvalue) {
				maxvalue = fabs(x[i]);
				maxind = i;
			}
		}
		return maxind;
	}
	
	// inf-norm (maximum absolute value) =========================================
	// technically not part of BLAS, but useful
	static inline T abs_max(const std::vector<T> &x)
	{ return std::abs(x[index_abs_max(x)]); }
	
	// saxpy (y=alpha*x+y) =======================================================
	static inline void add_scaled(T alpha, const std::vector<T> &x, std::vector<T> &y) {
		cblas_daxpy((Int)x.size(), alpha, &x[0], 1, &y[0], 1);
	}
};




template<class T>
void zero(std::vector<T> &v)
{ for(int i=(int)v.size()-1; i>=0; --i) v[i]=0; }

template<class T>
void insert(std::vector<T> &a, unsigned int index, T e)
{
   a.push_back(a.back());
   for(unsigned int i=(unsigned int)a.size()-1; i>index; --i)
      a[i]=a[i-1];
   a[index]=e;
}

template<class T>
void erase(std::vector<T> &a, unsigned int index)
{
   for(unsigned int i=index; i<a.size()-1; ++i)
      a[i]=a[i+1];
   a.pop_back();
}

//============================================================================
// Dynamic compressed sparse row matrix.

template<class T>
struct SparseMatrix
{
   int n; // dimension
   std::vector<std::vector<int> > index; // for each row, a list of all column indices (sorted)
   std::vector<std::vector<T> > value; // values corresponding to index

   explicit SparseMatrix(int n_=0, int expected_nonzeros_per_row=7)
      : n(n_), index(n_), value(n_)
   {
      for(int i=0; i<n; ++i){
         index[i].reserve(expected_nonzeros_per_row);
         value[i].reserve(expected_nonzeros_per_row);
      }
   }

   void clear(void)
   {
      n=0;
      index.clear();
      value.clear();
   }

   void zero(void)
   {
      for(int i=0; i<n; ++i){
         index[i].resize(0);
         value[i].resize(0);
      }
   }

   void resize(int n_)
   {
      n=n_;
      index.resize(n);
      value.resize(n);
   }

   T operator()(int i, int j) const
   {
      for(int k=0; k<(int)index[i].size(); ++k){
         if(index[i][k]==j) return value[i][k];
         else if(index[i][k]>j) return 0;
      }
      return 0;
   }

   void set_element(int i, int j, T new_value)
   {
      int k=0;
      for(; k<(int)index[i].size(); ++k){
         if(index[i][k]==j){
            value[i][k]=new_value;
            return;
         }else if(index[i][k]>j){
            insert(index[i], k, j);
            insert(value[i], k, new_value);
            return;
         }
      }
      index[i].push_back(j);
      value[i].push_back(new_value);
   }

   void add_to_element(int i, int j, T increment_value)
   {
      int k=0;
      for(; k<(int)index[i].size(); ++k){
         if(index[i][k]==j){
            value[i][k]+=increment_value;
            return;
         }else if(index[i][k]>j){
            insert(index[i], k, j);
            insert(value[i], k, increment_value);
            return;
         }
      }
      index[i].push_back(j);
      value[i].push_back(increment_value);
   }

   // assumes indices is already sorted
   void add_sparse_row(int i, const std::vector<int> &indices, const std::vector<T> &values)
   {
      int j=0, k=0;
      while(j<indices.size() && k<(int)index[i].size()){
         if(index[i][k]<indices[j]){
            ++k;
         }else if(index[i][k]>indices[j]){
            insert(index[i], k, indices[j]);
            insert(value[i], k, values[j]);
            ++j;
         }else{
            value[i][k]+=values[j];
            ++j;
            ++k;
         }
      }
      for(;j<indices.size(); ++j){
         index[i].push_back(indices[j]);
         value[i].push_back(values[j]);
      }
   }

   // assumes matrix has symmetric structure - so the indices in row i tell us which columns to delete i from
   void symmetric_remove_row_and_column(int i)
   {
      for(int a=0; a<index[i].size(); ++a){
         int j=index[i][a]; // 
         for(int b=0; b<index[j].size(); ++b){
            if(index[j][b]==i){
               erase(index[j], b);
               erase(value[j], b);
               break;
            }
         }
      }
      index[i].resize(0);
      value[i].resize(0);
   }

   void write_matlab(std::ostream &output, const char *variable_name)
   {
      output<<variable_name<<"=sparse([";
      for(int i=0; i<n; ++i){
         for(int j=0; j<index[i].size(); ++j){
            output<<i+1<<" ";
         }
      }
      output<<"],...\n  [";
      for(int i=0; i<n; ++i){
         for(int j=0; j<index[i].size(); ++j){
            output<<index[i][j]+1<<" ";
         }
      }
      output<<"],...\n  [";
      for(int i=0; i<n; ++i){
         for(int j=0; j<value[i].size(); ++j){
            output<<value[i][j]<<" ";
         }
      }
      output<<"], "<<n<<", "<<n<<");"<<std::endl;
   }
};

typedef SparseMatrix<float> SparseMatrixf;
typedef SparseMatrix<double> SparseMatrixd;

// perform result=matrix*x
template<class T>
void multiply(const SparseMatrix<T> &matrix, const std::vector<T> &x, std::vector<T> &result)
{
   assert(matrix.n==x.size());
   result.resize(matrix.n);
   //for(int i=0; i<matrix.n; ++i)
	parallel_for(matrix.n) {
		unsigned i (parallel_index);
		T value=0;
		for(int j=0; j<(int)matrix.index[i].size(); ++j){
			value+=matrix.value[i][j]*x[matrix.index[i][j]];
		}
		result[i]=value;
	} parallel_end
}

template <class T>
void mat_mat_multiply(const SparseMatrix<T> &m1, const SparseMatrix<T> &m2, SparseMatrix<T> &result) {
    for (int i = 0; i < m1.n; ++i) {
        for (int j = 0; j < m2.n; ++j) {
            T value = 0;
            int c = 0;
            for (const int k : m1.index[i]) {
                //value += m1.value[i][c] * m2(k, j);
                value += m1(i,k) * m2(k, j);
                c++;
            }
            if (value != 0) {
                result.set_element(i, j, value);
            }
        }
    }
}

template <class T>
void mat_transpose(const SparseMatrix<T> &m1, SparseMatrix<T> &result) {
    for (int i = 0; i < m1.n; ++i) {
        for (const int k : m1.index[i]) {
            result.set_element(k, i, m1(i, k));
        }
    }
}

// perform result=result-matrix*x
template<class T>
void multiply_and_subtract(const SparseMatrix<T> &matrix, const std::vector<T> &x, std::vector<T> &result)
{
   assert(matrix.n==x.size());
   result.resize(matrix.n);
   for(int i=0; i<(int)matrix.n; ++i){
      for(int j=0; j<(int)matrix.index[i].size(); ++j){
         result[i]-=matrix.value[i][j]*x[matrix.index[i][j]];
      }
   }
}

//============================================================================
// Fixed version of SparseMatrix. This is not a good structure for dynamically
// modifying the matrix, but can be significantly faster for matrix-vector
// multiplies due to better data locality.

template<class T>
struct FixedSparseMatrix
{
   int n; // dimension
   std::vector<T> value; // nonzero values row by row
   std::vector<int> colindex; // corresponding column indices
   std::vector<int> rowstart; // where each row starts in value and colindex (and last entry is one past the end, the number of nonzeros)

   explicit FixedSparseMatrix(int n_=0)
      : n(n_), value(0), colindex(0), rowstart(n_+1)
   {}

   void clear(void)
   {
      n=0;
      value.clear();
      colindex.clear();
      rowstart.clear();
   }

   void resize(int n_)
   {
      n=n_;
      rowstart.resize(n+1);
   }

   void construct_from_matrix(const SparseMatrix<T> &matrix)
   {
      resize(matrix.n);
      rowstart[0]=0;
      for(int i=0; i<n; ++i){
         rowstart[i+1]=rowstart[i]+matrix.index[i].size();
      }
      value.resize(rowstart[n]);
      colindex.resize(rowstart[n]);
      int j=0;
      for(int i=0; i<n; ++i){
         for(int k=0; k<(int)matrix.index[i].size(); ++k){
            value[j]=matrix.value[i][k];
            colindex[j]=matrix.index[i][k];
            ++j;
         }
      }
   }

   void write_matlab(std::ostream &output, const char *variable_name)
   {
      output<<variable_name<<"=sparse([";
      for(int i=0; i<n; ++i){
         for(int j=rowstart[i]; j<rowstart[i+1]; ++j){
            output<<i+1<<" ";
         }
      }
      output<<"],...\n  [";
      for(int i=0; i<n; ++i){
         for(int j=rowstart[i]; j<rowstart[i+1]; ++j){
            output<<colindex[j]+1<<" ";
         }
      }
      output<<"],...\n  [";
      for(int i=0; i<n; ++i){
         for(int j=rowstart[i]; j<rowstart[i+1]; ++j){
            output<<value[j]<<" ";
         }
      }
      output<<"], "<<n<<", "<<n<<");"<<std::endl;
   }
};



// perform result=matrix*x
template<class T>
void multiply(const FixedSparseMatrix<T> &matrix, const std::vector<T> &x, std::vector<T> &result)
{
   assert(matrix.n==x.size());
   result.resize(matrix.n);
   //for(int i=0; i<matrix.n; ++i)
	parallel_for(matrix.n) {
		unsigned i (parallel_index);
		T value=0;
		for(int j=matrix.rowstart[i]; j<matrix.rowstart[i+1]; ++j){
			value+=matrix.value[j]*x[matrix.colindex[j]];
		}
		result[i]=value;
	} parallel_end
}

template <class T>
void mat_mat_multiply(const FixedSparseMatrix<T> &m1, const FixedSparseMatrix<T> &m2, FixedSparseMatrix<T> &result) {
    result.resize(m1.n);
    result.rowstart[0] = 0;
    auto get_col = [&](const FixedSparseMatrix<T> &m, int row, int col) -> int {
        for (int i = m.rowstart[row]; i < m.rowstart[row + 1]; i++) {
            if (m.colindex[i] == col) {
                return i;
            }
        }
        return -1;
    };
    for (int i = 0; i < m1.n; ++i) {
        int added = 0;
        for (int j = 0; j < m2.n; ++j) {
            T value = 0;
            for (int k = m1.rowstart[i]; k < m1.rowstart[i + 1]; ++k) {
                auto col_m1 = m1.colindex[k];
                auto col = get_col(m2, col_m1, j);
                if (col == -1) {
                    continue;
                }
                value += m1.value[k] * m2.value[col];
            }
            if (value != 0) {
                result.value.push_back(value);
                result.colindex.push_back(j);
                added++;
            } 
        }
        result.rowstart[i + 1] = result.rowstart[i] + added;
    }
}

// perform result=result-matrix*x
template<class T>
void multiply_and_subtract(const FixedSparseMatrix<T> &matrix, const std::vector<T> &x, std::vector<T> &result)
{
   assert(matrix.n==x.size());
   result.resize(matrix.n);
   for(int i=0; i<matrix.n; ++i){
      for(int j=matrix.rowstart[i]; j<matrix.rowstart[i+1]; ++j){
         result[i]-=matrix.value[j]*x[matrix.colindex[j]];
      }
   }
}

//============================================================================
// A simple compressed sparse column data structure (with separate diagonal)
// for lower triangular matrices

template<class T>
struct SparseColumnLowerFactor
{
   int n;
   std::vector<T> invdiag; // reciprocals of diagonal elements
   std::vector<T> value; // values below the diagonal, listed column by column
   std::vector<int> rowindex; // a list of all row indices, for each column in turn
   std::vector<int> colstart; // where each column begins in rowindex (plus an extra entry at the end, of #nonzeros)
   std::vector<T> adiag; // just used in factorization: minimum "safe" diagonal entry allowed

   explicit SparseColumnLowerFactor(int n_=0)
      : n(n_), invdiag(n_), colstart(n_+1), adiag(n_)
   {}

   void clear(void)
   {
      n=0;
      invdiag.clear();
      value.clear();
      rowindex.clear();
      colstart.clear();
      adiag.clear();
   }

   void resize(int n_)
   {
      n=n_;
      invdiag.resize(n);
      colstart.resize(n+1);
      adiag.resize(n);
   }

   void write_matlab(std::ostream &output, const char *variable_name)
   {
      output<<variable_name<<"=sparse([";
      for(int i=0; i<n; ++i){
         output<<" "<<i+1;
         for(int j=colstart[i]; j<colstart[i+1]; ++j){
            output<<" "<<rowindex[j]+1;
         }
      }
      output<<"],...\n  [";
      for(int i=0; i<n; ++i){
         output<<" "<<i+1;
         for(int j=colstart[i]; j<colstart[i+1]; ++j){
            output<<" "<<i+1;
         }
      }
      output<<"],...\n  [";
      for(int i=0; i<n; ++i){
         output<<" "<<(invdiag[i]!=0 ? 1/invdiag[i] : 0);
         for(int j=colstart[i]; j<colstart[i+1]; ++j){
            output<<" "<<value[j];
         }
      }
      output<<"], "<<n<<", "<<n<<");"<<std::endl;
   }
};

//============================================================================
// Incomplete Cholesky factorization, level zero, with option for modified version.
// Set modification_parameter between zero (regular incomplete Cholesky) and
// one (fully modified version), with values close to one usually giving the best
// results. The min_diagonal_ratio parameter is used to detect and correct
// problems in factorization: if a pivot is this much less than the diagonal
// entry from the original matrix, the original matrix entry is used instead.

template<class T>
void factor_modified_incomplete_cholesky0(const SparseMatrix<T> &matrix, SparseColumnLowerFactor<T> &alpha,
                                          T modification_parameter=0.97, T min_diagonal_ratio=0.25)
{
   // first copy lower triangle of matrix into factor (Note: assuming A is symmetric of course!)
   alpha.resize(matrix.n);
   zero(alpha.invdiag); // important: eliminate old values from previous solves!
   alpha.value.resize(0);
   alpha.rowindex.resize(0);
   zero(alpha.adiag);
   for(int i=0; i<matrix.n; ++i){
      alpha.colstart[i]=(int)alpha.rowindex.size();
      for(int j=0; j<(int)matrix.index[i].size(); ++j){
         if(matrix.index[i][j]>i){
            alpha.rowindex.push_back(matrix.index[i][j]);
            alpha.value.push_back(matrix.value[i][j]);
         }else if(matrix.index[i][j]==i){
            alpha.invdiag[i]=alpha.adiag[i]=matrix.value[i][j];
         }
      }
   }
   alpha.colstart[matrix.n]=(int)alpha.rowindex.size();
   // now do the incomplete factorization (figure out numerical values)

   // MATLAB code:
   // L=tril(A);
   // for k=1:size(L,2)
   //   L(k,k)=sqrt(L(k,k));
   //   L(k+1:end,k)=L(k+1:end,k)/L(k,k);
   //   for j=find(L(:,k))'
   //     if j>k
   //       fullupdate=L(:,k)*L(j,k);
   //       incompleteupdate=fullupdate.*(A(:,j)~=0);
   //       missing=sum(fullupdate-incompleteupdate);
   //       L(j:end,j)=L(j:end,j)-incompleteupdate(j:end);
   //       L(j,j)=L(j,j)-omega*missing;
   //     end
   //   end
   // end

   for(int k=0; k<matrix.n; ++k){
      if(alpha.adiag[k]==0) continue; // null row/column
      // figure out the final L(k,k) entry
      if(alpha.invdiag[k]<min_diagonal_ratio*alpha.adiag[k])
         alpha.invdiag[k]=1/sqrt(alpha.adiag[k]); // drop to Gauss-Seidel here if the pivot looks dangerously small
      else
         alpha.invdiag[k]=1/sqrt(alpha.invdiag[k]);
      // finalize the k'th column L(:,k)
      for(int p=alpha.colstart[k]; p<alpha.colstart[k+1]; ++p){
         alpha.value[p]*=alpha.invdiag[k];
      }
      // incompletely eliminate L(:,k) from future columns, modifying diagonals
      for(int p=alpha.colstart[k]; p<alpha.colstart[k+1]; ++p){
         int j=alpha.rowindex[p]; // work on column j
         T multiplier=alpha.value[p];
         T missing=0;
         int a=alpha.colstart[k];
         // first look for contributions to missing from dropped entries above the diagonal in column j
         int b=0;
         while(a<alpha.colstart[k+1] && alpha.rowindex[a]<j){
            // look for factor.rowindex[a] in matrix.index[j] starting at b
            while(b<(int)matrix.index[j].size()){
               if(matrix.index[j][b]<alpha.rowindex[a])
                  ++b;
               else if(matrix.index[j][b]==alpha.rowindex[a])
                  break;
               else{
                  missing+=alpha.value[a];
                  break;
               }
            }
            ++a;
         }
         // adjust the diagonal j,j entry
         if(a<alpha.colstart[k+1] && alpha.rowindex[a]==j){
            alpha.invdiag[j]-=multiplier*alpha.value[a];
         }
         ++a;
         // and now eliminate from the nonzero entries below the diagonal in column j (or add to missing if we can't)
         b=alpha.colstart[j];
         while(a<alpha.colstart[k+1] && b<alpha.colstart[j+1]){
            if(alpha.rowindex[b]<alpha.rowindex[a])
               ++b;
            else if(alpha.rowindex[b]==alpha.rowindex[a]){
               alpha.value[b]-=multiplier*alpha.value[a];
               ++a;
               ++b;
            }else{
               missing+=alpha.value[a];
               ++a;
            }
         }
         // and if there's anything left to do, add it to missing
         while(a<alpha.colstart[k+1]){
            missing+=alpha.value[a];
            ++a;
         }
         // and do the final diagonal adjustment from the missing entries
         alpha.invdiag[j]-=modification_parameter*multiplier*missing;
      }
   }
}

//============================================================================
// Solution routines with lower triangular matrix.

// solve L*result=rhs
template<class T>
void solve_lower(const SparseColumnLowerFactor<T> &alpha, const std::vector<T> &rhs, std::vector<T> &result)
{
   assert(alpha.n==rhs.size());
   assert(alpha.n==result.size());
   result=rhs;
   for(int i=0; i<alpha.n; ++i){
      result[i]*=alpha.invdiag[i];
      for(int j=alpha.colstart[i]; j<alpha.colstart[i+1]; ++j){
         result[alpha.rowindex[j]]-=alpha.value[j]*result[i];
      }
   }
}

// solve L^T*result=rhs
template<class T>
void solve_lower_transpose_in_place(const SparseColumnLowerFactor<T> &alpha, std::vector<T> &x)
{
   assert(alpha.n==(int)x.size());
   assert(alpha.n>0);
   int i=alpha.n;
   do{
      --i;
      for(int j=alpha.colstart[i]; j<alpha.colstart[i+1]; ++j){
         x[i]-=alpha.value[j]*x[alpha.rowindex[j]];
      }
      x[i]*=alpha.invdiag[i];
   }while(i!=0);
}




//============================================================================
// Encapsulates the Conjugate Gradient algorithm with incomplete Cholesky
// factorization preconditioner.

template <class T>
struct SparsePCGSolver
{
   SparsePCGSolver(void)
   {
      set_solver_parameters(1e-5, 100, 0.97, 0.25);
   }

   void set_solver_parameters(T tolerance_factor_, int max_iterations_, T modified_incomplete_cholesky_parameter_=0.97, T min_diagonal_ratio_=0.25)
   {
      tolerance_factor=tolerance_factor_;
      if(tolerance_factor<1e-30) tolerance_factor=1e-30;
      max_iterations=max_iterations_;
      modified_incomplete_cholesky_parameter=modified_incomplete_cholesky_parameter_;
      min_diagonal_ratio=min_diagonal_ratio_;
   }

   bool solve(const SparseMatrix<T> &matrix, const std::vector<T> &rhs, std::vector<T> &result, T &relative_residual_out, int &iterations_out, int precondition=2)
   {
      int n=matrix.n;
      if((int)m.size()!=n){ m.resize(n); s.resize(n); z.resize(n); r.resize(n); }
      zero(result);
      r=rhs;
      double residual_out=InstantBLAS<int,T>::abs_max(r);
      if(residual_out==0) {
         iterations_out=0;
         return true;
      }
      //double tol=tolerance_factor*residual_out; // relative residual
      double tol=tolerance_factor;
	  double residual_0 = residual_out;

      form_preconditioner(matrix, precondition);
      apply_preconditioner( r, z, precondition);
      double rho=InstantBLAS<int,T>::dot(z, r);
      if(rho==0 || rho!=rho) {
         iterations_out=0;
         return false;
      }

      s=z;
      fixed_matrix.construct_from_matrix(matrix);
      int iteration;
      for(iteration=0; iteration<max_iterations; ++iteration){
         multiply(fixed_matrix, s, z);
         double alpha=rho/InstantBLAS<int,T>::dot(s, z);
         InstantBLAS<int,T>::add_scaled(-alpha, s, result);
         InstantBLAS<int,T>::add_scaled(-alpha, z, r);
         residual_out=InstantBLAS<int,T>::abs_max(r);
		 relative_residual_out = residual_out / residual_0;
         if(residual_out<=tol) {
            iterations_out=iteration+1;
            return true; 
         }
         apply_preconditioner(r, z, precondition);
         double rho_new=InstantBLAS<int,T>::dot(z, r);
         double beta=rho_new/rho;
         InstantBLAS<int,T>::add_scaled(beta, s, z); s.swap(z); // s=beta*s+z
         rho=rho_new;
       
      }
      iterations_out=iteration;
	  relative_residual_out = residual_out / residual_0;
      return false;
   }
public:
    FixedSparseMatrix<T> fixed_matrix; // used within loop
   protected:
   // internal structures
   SparseColumnLowerFactor<T> ic_factor; // modified incomplete cholesky factor
   std::vector<T> m, z, s, r; // temporary vectors for PCG


   // parameters
   T tolerance_factor;
   int max_iterations;
   T modified_incomplete_cholesky_parameter;
   T min_diagonal_ratio;

   void form_preconditioner(const SparseMatrix<T>& matrix, int precondition=2)
   {
		if(precondition==2) {
			// incomplete cholesky
			factor_modified_incomplete_cholesky0(matrix, ic_factor, modified_incomplete_cholesky_parameter, min_diagonal_ratio);
		
		} else if(precondition==1) {
			// diagonal
		    ic_factor.resize(matrix.n);
		    zero(ic_factor.invdiag);
			for(int i=0; i<matrix.n; ++i) {
				  for(int j=0; j<(int)matrix.index[i].size(); ++j){
					 if(matrix.index[i][j]==i){
						ic_factor.invdiag[i] = 1./matrix.value[i][j];
					 }
				  }
			}
		}
   }

   void apply_preconditioner(const std::vector<T> &x, std::vector<T> &result, int precondition=2)
   {
		if (precondition==2) {
			// incomplete cholesky
			solve_lower(ic_factor, x, result);
			solve_lower_transpose_in_place(ic_factor,result);
		} else if(precondition==1) {
			// diagonal
			for(int_index i=0; i<(int_index)result.size(); ++i) {
				result[i] = x[i] * ic_factor.invdiag[i];
		    }
		} else {
			// off
			result = x;
		}
   }
};



#undef parallel_for
#undef parallel_end             
#undef int_index

#undef parallel_block 
#undef do_parallel 
#undef do_end
#undef block_end

#endif