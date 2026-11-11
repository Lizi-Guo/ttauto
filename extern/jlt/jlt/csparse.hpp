//
// Copyright (c) 2004-2026 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

//
// Routines to interface with the CSparse library.
//
// CSparse is Copyright (c) Timothy A. Davis, 2006-2012
//
// See http://www.cise.ufl.edu/research/sparse/CSparse
//
// Compile with -lcsparse flag.
//

#ifndef JLT_CSPARSE_HPP
#define JLT_CSPARSE_HPP

#include <iostream>
#include <cstdio>
#include <memory>
#include <jlt/mathmatrix.hpp>

namespace csparse
{
extern "C"
{
#include "cs.h"
}
} // namespace csparse

namespace jlt {

// unique_ptr wrappers for csparse pointers that will take care of
// freeing the memory when we're done with a matrix.

// CSparse allocates with malloc/calloc. Release the complete object through
// its C API instead of invoking C++ delete on the allocation.
struct cs_deleter {
  void operator()(csparse::cs* value) const noexcept { csparse::cs_spfree(value); }
};
struct csd_deleter {
  void operator()(csparse::csd* value) const noexcept { csparse::cs_dfree(value); }
};
class cs_unique_ptr : public std::unique_ptr<csparse::cs, cs_deleter> {
public:
  cs_unique_ptr(csparse::cs* value = nullptr)
    : std::unique_ptr<csparse::cs, cs_deleter>(value) {}
  operator csparse::cs*() const { return get(); }
  bool operator!() const { return get() == nullptr; }
};
class csd_unique_ptr : public std::unique_ptr<csparse::csd, csd_deleter> {
public:
  csd_unique_ptr(csparse::csd* value = nullptr)
    : std::unique_ptr<csparse::csd, csd_deleter>(value) {}
  operator csparse::csd*() const { return get(); }
  bool operator!() const { return get() == nullptr; }
};

// Note that we return a pointer, not an unique_ptr, since the return
// value is a temporary.
template<class T>
csparse::cs* mathmatrix_to_cs_sparse_matrix(const mathmatrix<T>& M)
{
  int m = M.rows(), n = M.columns();
  cs_unique_ptr csM(csparse::cs_spalloc(0,0,1,1,1));

  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      if (M(i,j) != 0)
	csparse::cs_entry(csM,i,j,M(i,j));

  // Column-compress the sparse matrix.
  return csparse::cs_compress(csM);
}

template<class T>
mathmatrix<T> cs_sparse_matrix_to_mathmatrix(const cs_unique_ptr& csM)
{
  if (!csM) return mathmatrix<T>();

  int m = csM->m, n = csM->n;
  auto* Ap = csM->p; auto* Ai = csM->i;
  int nz = csM->nz;
  double *Ax = csM->x;

  mathmatrix<T> M(m,n);

  if (nz < 0)
    {
      // Matrix is in column-compressed form.
      for (int j = 0; j < n; ++j)
	  for (int p = Ap[j]; p < Ap[j+1]; ++p)
	      M(Ai[p],j) = T(Ax[p]);
    }
  else
    {
      // Matrix is in triplet form.
      for (int p = 0; p < nz; ++p)
	M(Ai[p],Ap[p]) = T(Ax[p]);
    }

  return M;
}

} // namespace jlt

#endif // JLT_CSPARSE_HPP
