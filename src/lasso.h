/* Copyright (C) 1998
   Berwin A Turlach <bturlach@stats.adelaide.edu.au> */

/* This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public License
   as published by the Free Software Foundation; either version 2 of
   the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details. */

/* You should have received a copy of the GNU Library General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston,
   MA 02111-1307, USA. */
#define S_or_R

#ifndef BT_LASSO_H
#define BT_LASSO_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#if defined(S_or_R)
# include <S.h>
# if defined(SPLUS_VERSION) && SPLUS_VERSION >= 4000
#  include <newredef.h>
# endif
# ifdef USING_R
   typedef int Sint;
# else
   typedef long Sint;
# endif

#else /* !defined(S_or_R) */
# define Calloc(n,t)       malloc((n)*sizeof(t))
# define Realloc(p,n,t)    realloc(p,(n)*sizeof(t))
# define Free(p)           free(p), p=NULL
# include <strings.h>
# define Memcpy(p,q,n)     memcpy(p,q,(n)*sizeof(*(p)))
#endif

#ifndef USING_R
# define Sprintf Rprintf
# else
# define Sprintf printf
#endif

void lasso(double *x, Sint *pn, Sint *pm, double *pt,
           double *beta, double *y, double *yhat1, double *r,
           double *lagrangian, Sint *psuc,  Sint *pverb, Sint *pas_sub);
void mult_lasso(double *x, Sint *pn, Sint *pm, double *pt, Sint *pl,
                double *beta, double *y, double *yhat1, double *r,
                double *lagrangian, Sint *psuc, Sint *pverb);
void rs_lasso(double *x, double *y, Sint *pn, Sint *porder,
           double *pt, double *beta, double *yhat1,  double *yhat2,
           double *r, Sint *psuc, Sint *pverb, Sint *pas_sub);
void rs_auto(double *x, double *y, Sint *pn, Sint *porder,
           double *pt, double *beta, double *yhat1, double *yhat2,
           double *r,  Sint *pcrit1, double *pprec,
           Sint *psuc, Sint *pverb,
           double *tgr, double *tgr_val, Sint *ptgr_len);

#define TRUE 1
#define FALSE 0
#define RMAT(i,j) (*(rmat+(j)*((j)+1)/2+(i)))
#define QMAT(i,j) (*(qmat+(j)*q_nrow+(i)))
#define RLAST(i)  (*(rmat+(r_ncol-1)*r_ncol/2+(i)))
#define RLTEL     (*(rmat+r_ncol*(r_ncol+1)/2-1))
#define RCOL(j)   (rmat+(j)*((j)+1)/2)
#define QCOL(j)   (qmat+(j)*q_nrow)
#define RLTCOL    (rmat+(r_ncol-1)*r_ncol/2)
#define QR_CHUNK 10
#if defined (S_or_R)
# define ERRMSG(where,msg) errmsg(msg)
#else
# define ERRMSG(where,msg) errmsg(where,msg)
#endif
#endif
