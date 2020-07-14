///
/// \file classical_sort.h
/// \author Charles Bouillaguet & Jérôme Bonacchi
/// \date 2020-07-09
/// 
/// @copyright Copyright (c) 2020
///

#ifndef INCLUDE_CLASSICAL_SORT_H_W
#define INCLUDE_CLASSICAL_SORT_H_W

#include "mini_spasm_wang.h"

void classical_compress(const spasm_triplet * T, spasm * A, int *W);
void classical_transpose(const spasm * A, spasm * R, int *W);

//#endif /* INCLUDE_CLASSICAL_SORT_H */
#include <stdio.h>                                                               
#include <stdlib.h>                                                              
#include <assert.h>                                                              
#include <stdbool.h>                                                             
#include <err.h>                                                                 
                                                                                 
                                                      
                                                          
                                                                                 
// assemble the CSR representation of T                                          
// W = scratch space, size == #rows + 1                                          
void classical_compress(const spasm_triplet * T, spasm * A, int *W)              
{                                                                                
  const int *Ti = T->i;                                                          
  const int *Tj = T->j;                                                          
  const double *Tx = T->x;                                                       
  int *Ap = A->p;                                                                
  int *Aj = A->j;                                                                
  double *Ax = A->x;                                                             
                                                                                 
  int n = T->n;                                                                  
  int nnz = T->nz;                                                               
                                                                                 
  // setup                                                                       
  for (int i = 0; i < n; i++)  /* gcc replaces this with AVX2-optimized memset */
    W[i] = 0;                                                                    
                                                                                 
  // counting entries on each row                                                
  for (int k = 0; k < nnz; k++) {                                                
    int i = Ti[k];                                                               
    int w = W[i];                                                                
    w += 1;                                                                      
    W[i] = w;                                                                    
  }                                                                              
                                                                                 
  // prefix-sum + setup Ap                                                       
  int s = 0;                                                                     
  for (int i = 0; i < n; i++) {                                                  
    int w = W[i];                                                                
    Ap[i] = s;                                                                   
    W[i] = s;                                                                    
    s += w;                                                                      
  }                                                                              
  Ap[n] = s;                                                                     
                                                                                 
  // dispatch                                                                    
  for (int k = 0; k < nnz; k++) {                                                
    int i = Ti[k];                                                               
    int j = Tj[k];                                                               
    double x = Tx[k];                                                            
    int l = W[i];                                                                
    Aj[l] = j;                                                                   
    Ax[l] = x;
    l += 1;                                                                      
    W[i] = l;                                                                    
  }                                                                              
}
#endif
