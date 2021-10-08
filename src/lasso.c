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

#include "lasso.h"

static void lasso_alloc(Sint n, Sint m);
static void lasso_free(void);
static void qr_init(int n);
static void qr_incr(void);
static void qr_free(void);
static void qr_del(int l, int aug);
static void qr_add(double *x, int swap);
#if defined (S_or_R)
static void errmsg(char* string);
#else
static void errmsg(char* where, char* string);
#endif


/* Global Variables [yuck!] : */

  static double *xtr=NULL, *btmp=NULL, *qtr=NULL,
    *rinvt_theta=NULL, *step=NULL, ytyd2=0.0;
  static int *theta=NULL, *nz_x=NULL, num_nz_x=0;
  static double *qmat=NULL, *rmat=NULL;
  static int qr_max_size=0, r_ncol=0, q_nrow=0, q_use_row=0;
  static char *no_dyn_mem_message="Cannot allocate dynamic memory";


void lasso(double *x, Sint *pn, Sint *pm, double *pt,
           double *beta, double *y, double *yhat1, double *r,
           double *lagrangian, Sint *psuc, Sint *pverb, Sint *pas_sub)
{

    double t = *pt, prec;
    Sint n = *pn, m = *pm, verb = *pverb, as_sub = *pas_sub;
    int not_solved;
    double *x_elem=NULL, tmp, max_val, b_1norm;
    int i, j, max_ind;
    double p_obj, d_obj;
    int num_iter=0, max_num;
    double b_1norm_old, mu;
    double rho_up, rho_low, rho;
    int to_del;
    double *q_elem, wtc, wtw;
    double *r_elem;
    int k;
    int add;
    if( !as_sub)
	lasso_alloc(n,m);

    prec = sqrt(DBL_EPSILON);
    if( as_sub ) {

	b_1norm = 0.0;
	for(j=0;j < num_nz_x;j++)
	    b_1norm += fabs(beta[nz_x[j]]);
    } else {

	b_1norm = 0.0;
	num_nz_x = 0;
	for(j=0; j < m; j++) {
	    if(fabs(beta[j]) > prec) {
		b_1norm += fabs(beta[j]);
		nz_x[num_nz_x] = j;
		num_nz_x++;
	    } else
		beta[j] = 0.0;
	}
    }
    if( b_1norm > t) {
	if(verb) {
	    Sprintf("******************************\n");
	    Sprintf("Rescaling beta from L1-norm %f to %f\n",b_1norm,t);
	}
	for(j=0; j < num_nz_x; j++)
	    beta[nz_x[j]] = beta[nz_x[j]] * t/b_1norm;
	b_1norm = t;
    }


    for(i=0; i <  n; i++)
	yhat1[i] = 0.0;
    for(j=0; j < num_nz_x; j++) {
	/* x_elem points to the first element in the column of X to which
	   the j-th entry in nz_x points  */
	x_elem = x + nz_x[j]*n;
	tmp = beta[nz_x[j]];
	for(i=0; i < n; i++) {
	    yhat1[i] += *x_elem*tmp;
	    x_elem++;           /* now we point to the next element in X */
	}
    }

    /* calculate the residual vector */
    for(i=0; i <  n; i++)
	r[i] = y[i]-yhat1[i];


    /* multiply X^T with the residual vector */
    x_elem = x;
    for(j=0; j < m; j++) {
	tmp = 0.0;
	for(i=0; i < n; i++) {
	    tmp += *x_elem*r[i];
	    x_elem++;
	}
	xtr[j] = tmp;
    }
    max_val = fabs(xtr[0]);
    max_ind = 0;
    for(j=1; j < m; j++)
	if( fabs(xtr[j]) > max_val ) {
	    max_val = fabs(xtr[j]);
	    max_ind = j;
	}

    if( !as_sub ) {

	qr_add(y,TRUE);
	ytyd2 = *rmat * *rmat/2.0;
	for(j=0;j < num_nz_x;j++) {
	    qr_add(x+nz_x[j]*n, TRUE);
	    if(fabs(beta[nz_x[j]]) < prec)
		theta[j] = xtr[nz_x[j]] < 0 ? -1 : 1;
	    else
		theta[j] = beta[nz_x[j]] < 0 ? -1 : 1;
	}
    }
    if( num_nz_x==0 ) {
	nz_x[0] = max_ind;
	num_nz_x = 1;
	if(verb) {
	    Sprintf("******************************\n");
	    Sprintf("  -->\tAdding variable: %d\n",max_ind+1);
	}
	qr_add(x+max_ind*n, TRUE);
	theta[0] = xtr[max_ind] < 0 ? -1 : 1;
    }
    *psuc=0;
    if(verb) {


	/* Find out how many times [[max_val]] is attained */
	tmp = (1.0-prec)*max_val;
	if( tmp < prec ) tmp = 0.0;
	max_num = 1;
	for(j=0; j < m; j++) {
	    if( tmp <= fabs(xtr[j]) && j!=max_ind) {
		/*we found another element equal to the (current)
		  maximal absolute value */
		max_num++;
	    }
	}

	/* for the value of the dual objective function we need
	   to calculate the L2 norm of the vector of fitted values */
	p_obj = 0.0;
	d_obj = 0.0;
	for(i=0;i < n;i++) {
	    p_obj += r[i]*r[i];
	    d_obj += yhat1[i]*yhat1[i];
	}
	p_obj /= 2.0;
	d_obj = ytyd2 - d_obj/2.0 - t*max_val ;
	Sprintf("******************************\n");
	Sprintf("\nIteration number: %d\n", num_iter);
	Sprintf("Value of primal object function      : %f\n", p_obj);
	Sprintf("Value of dual   object function      : %f\n", d_obj);
	Sprintf("L1 norm of current beta              : %f", b_1norm);
	Sprintf(" <= %f\n", t);
	Sprintf("Maximal absolute value in t(X)%%*%%r   : %e", max_val);
	Sprintf(" attained %d time(s)\n", max_num);
	Sprintf("Number of parameters allowed to vary : %d\n", num_nz_x);
	num_iter++;
    }
    while(1) {
	do {


	    q_elem = qmat;
	    for(j=0;j < num_nz_x;j++) {
		tmp = 0.0;
		for(i=0;i < n;i++,q_elem++)
		    tmp += *q_elem*r[i];
		qtr[j] = tmp;
	    }
	    if( b_1norm < (1.0-prec)*t )
		mu = 0.0;
	    else {

		/* z=R^{-T}\theta can be calculated by solving R^Tz=\theta */
		r_elem = rmat;
		for(j=0;j < num_nz_x;j++,r_elem++) {
		    tmp = theta[j];
		    for(k=0;k < j;k++,r_elem++)
			tmp -= *r_elem * rinvt_theta[k];
		    rinvt_theta[j] = tmp /  *r_elem;
		}

		wtc = 0.0;
		wtw = 0.0;
		for(j=0;j < num_nz_x;j++) {
		    wtc += rinvt_theta[j]*qtr[j];
		    wtw += rinvt_theta[j]*rinvt_theta[j];
		}
		mu = (wtc-(t-b_1norm))/wtw;
	    }

	    if(mu<=0.0)
		for(j=0;j < num_nz_x;j++)
		    step[j] = qtr[j];
	    else
		for(j=0;j < num_nz_x;j++)
		    step[j] = qtr[j] - mu*rinvt_theta[j];

	    /* h=R^{-1}z can be calculated by solving Rh=z */
	    for(j=num_nz_x-1;j>=0;j--) {
		tmp = step[j];
		for(k=num_nz_x-1;k>j;k--)
		    tmp -= RMAT(j,k) * step[k];
		step[j] = tmp /  RMAT(j,j);
	    }
	    for(j=0;j < num_nz_x;j++) {
		btmp[j] = beta[nz_x[j]];
		beta[nz_x[j]] += step[j];
	    }
	    b_1norm_old=b_1norm;

	    b_1norm = 0.0;
	    for(j=0;j < num_nz_x;j++)
		b_1norm += fabs(beta[nz_x[j]]);

	    not_solved=FALSE;
	    if( b_1norm > (1+prec)*t) {
		not_solved=TRUE;
		if(b_1norm_old < (1.0-prec)*t) {

		    if(verb)
			Sprintf("  -->\tStepping onto the border of the L1 ball.\n");
		    rho_up  = rho = 1.0;
		    rho_low = 0.0;
		    while( fabs(t-b_1norm) > prec*t ) {
			if( b_1norm > t) {
			    rho_up = rho;
			    rho = (rho+rho_low)/2.0;
			}
			if( b_1norm < t) {
			    rho_low = rho;
			    rho = (rho+rho_up)/2.0;
			}
			if(rho < prec) break;
			for(j=0; j < num_nz_x; j++)
			    beta[nz_x[j]] = btmp[j]+rho*step[j];

			b_1norm = 0.0;
			for(j=0;j < num_nz_x;j++)
			    b_1norm += fabs(beta[nz_x[j]]);
		    }
		    for(j=0;j < num_nz_x;j++) {
			if(fabs(beta[nz_x[j]]) < prec)
			    theta[j] = btmp[j] > 0 ? -1 : 1;
			else
			    theta[j] = beta[nz_x[j]] < 0 ? -1 : 1;
		    }
		} else {


		    rho = 1.0;
		    to_del = -1;
		    for(j=0; j < num_nz_x; j++) {
			if(fabs(step[j]) > prec) {
			    tmp = -btmp[j]/(step[j]);
			    if( 0.0 < tmp && tmp < rho ) {
				rho = tmp;
				to_del = j;
			    }
			}
		    }
		    if(to_del < 0 ) {
			*psuc= -1;
			goto EXIT_HERE;
		    }
		    for(j=0; j<num_nz_x; j++)
			beta[nz_x[j]] = btmp[j]+rho*step[j];
		    if(verb)
			Sprintf("  -->\tRemoving variable: %d",nz_x[to_del]+1);
		    beta[nz_x[to_del]] = 0.0;
		    qr_del(to_del,TRUE);
		    for(j=to_del+1; j< num_nz_x; j++) {
			nz_x[j-1] = nz_x[j];
			theta[j-1] = theta[j];
		    }
		    num_nz_x--;

		    b_1norm = 0.0;
		    for(j=0;j<num_nz_x;j++)
			b_1norm += fabs(beta[nz_x[j]]);
		    if(verb) {
			if(b_1norm < (1-prec)*t)
			    Sprintf(", and stepping into the interior of the L1 ball\n");
			else
			    Sprintf("\n");
		    }
		}
	    }

	    for(i=0; i< n; i++)
		yhat1[i] = 0.0;
	    for(j=0; j < num_nz_x; j++) {
		/* x_elem points to the first element in the column of X to which
		   the j-th entry in nz_x points  */
		x_elem = x + nz_x[j]*n;
		tmp = beta[nz_x[j]];
		for(i=0; i < n; i++) {
		    yhat1[i] += *x_elem*tmp;
		    x_elem++;           /* now we point to the next element in X */
		}
	    }

	    /* calculate the residual vector */
	    for(i=0; i< n; i++)
		r[i] = y[i]-yhat1[i];
	}while(not_solved);


	for(i=0; i< n; i++)
	    yhat1[i] = 0.0;
	for(j=0; j < num_nz_x; j++) {
	    /* x_elem points to the first element in the column of X to which
	       the j-th entry in nz_x points  */
	    x_elem = x + nz_x[j]*n;
	    tmp = beta[nz_x[j]];
	    for(i=0; i < n; i++) {
		yhat1[i] += *x_elem*tmp;
		x_elem++;           /* now we point to the next element in X */
	    }
	}

	/* calculate the residual vector */
	for(i=0; i< n; i++)
	    r[i] = y[i]-yhat1[i];


	/* multiply X^T with the residual vector */
	x_elem = x;
	for(j=0; j < m; j++) {
	    tmp = 0.0;
	    for(i=0; i<n; i++) {
		tmp += *x_elem*r[i];
		x_elem++;
	    }
	    xtr[j] = tmp;
	}
	max_val = fabs(xtr[0]);
	max_ind = 0;
	for(j=1; j<m; j++)
	    if( fabs(xtr[j]) > max_val ) {
		max_val = fabs(xtr[j]);
		max_ind = j;
	    }

	if(verb) {


	    /* Find out how many times [[max_val]] is attained */
	    tmp = (1.0-prec)*max_val;
	    if( tmp < prec ) tmp = 0.0;
	    max_num = 1;
	    for(j=0; j<m; j++) {
		if( tmp <= fabs(xtr[j]) && j!=max_ind) {
		    /*we found another element equal to the (current)
		      maximal absolute value */
		    max_num++;
		}
	    }

	    /* for the value of the dual objective function we need
	       to calculate the L2 norm of the vector of fitted values */
	    p_obj = 0.0;
	    d_obj = 0.0;
	    for(i=0;i<n;i++) {
		p_obj += r[i]*r[i];
		d_obj += yhat1[i]*yhat1[i];
	    }
	    p_obj /= 2.0;
	    d_obj = ytyd2 - d_obj/2.0 - t*max_val ;
	    Sprintf("******************************\n");
	    Sprintf("\nIteration number: %d\n", num_iter);
	    Sprintf("Value of primal object function      : %f\n", p_obj);
	    Sprintf("Value of dual   object function      : %f\n", d_obj);
	    Sprintf("L1 norm of current beta              : %f", b_1norm);
	    Sprintf(" <= %f\n", t);
	    Sprintf("Maximal absolute value in t(X)%%*%%r   : %e", max_val);
	    Sprintf(" attained %d time(s)\n", max_num);
	    Sprintf("Number of parameters allowed to vary : %d\n", num_nz_x);
	    num_iter++;
	}

	add = TRUE;
	for(i=0; i<num_nz_x; i++)
	    if(nz_x[i]==max_ind) {
		add = FALSE;
		break;
	    }
	if(add) {
	    qr_add(x+max_ind*n,TRUE);
	    nz_x[num_nz_x] = max_ind;
	    theta[num_nz_x] = xtr[max_ind]<0 ? -1 : 1;
	    num_nz_x++;
	} else {
	    *lagrangian = max_val;
	    break;
	}

	b_1norm = 0.0;
	for(j=0; j < num_nz_x; j++)
	    b_1norm += fabs(beta[nz_x[j]]);
	if(verb)
	    Sprintf("  -->\tAdding variable: %d\n",max_ind+1);
    }
EXIT_HERE:
    if( !as_sub)
	lasso_free();
}
void mult_lasso(double *x, Sint *pn, Sint *pm, double *pt, Sint *pl,
                double *beta, double *y, double *yhat1, double *r,
                double *lagrangian, Sint *psuc, Sint *pverb)
{

    double prec;
    Sint n = *pn, m = *pm, l = *pl, verb = *pverb, as_sub = TRUE;
    int i, j;
    lasso_alloc(n,m);

    qr_add(y,TRUE);
    ytyd2 = *rmat * *rmat/2.0;
    prec = sqrt(DBL_EPSILON);
    num_nz_x = 0;
    for(j=0; j<m; j++) {
	if(fabs(beta[j]) > prec) {
	    qr_add(x+j*n, TRUE);
	    nz_x[num_nz_x] = j;
	    num_nz_x++;
	} else
	    beta[j] = 0.0;
    }
    *psuc = 0;
    for(i=0; i<l; i++) {

	if(verb) {
	    Sprintf("\n\n++++++++++++++++++++++++++++++\n");
	    Sprintf("Solving problem number %d with bound %f\n", i+1, pt[i]);
	    Sprintf("++++++++++++++++++++++++++++++\n");
	}
	if(i>0) Memcpy(beta,beta-m,m);
	lasso(x, pn, pm, pt+i, beta, y, yhat1, r, lagrangian, psuc, pverb, &as_sub);
	if( *psuc < 0 ) {
	    goto EXIT_HERE;
	}

	beta += m;
	yhat1 += n;
	r += n;
	lagrangian++;
    }
EXIT_HERE:
    lasso_free();
}
static void lasso_alloc(Sint n, Sint m)
{

#if defined(S_or_R)
    if( nz_x != NULL || theta != NULL || xtr != NULL || btmp != NULL ||
	qtr != NULL || rinvt_theta != NULL || step != NULL ||
	num_nz_x != 0 || ytyd2 != 0.0) {
	LASSO2_MESSAGE "Possible memory corruption or memory leak.\n  We"
	    "advise to restart your S+ session"  LASSO2_WARNING(LASSO2_NULL_ENTRY);
	lasso_free();
    }
#endif
    nz_x = Calloc(m,int);
    if( nz_x == NULL )
	ERRMSG("lasso_alloc", no_dyn_mem_message);
    theta = Calloc(m,int);
    if( theta==NULL )
	ERRMSG("lasso_alloc", no_dyn_mem_message);
    xtr = Calloc(m,double);
    if( xtr==NULL )
	ERRMSG("lasso_alloc", no_dyn_mem_message);
    btmp = Calloc(m,double);
    if( btmp==NULL )
	ERRMSG("lasso_alloc", no_dyn_mem_message);
    qtr = Calloc(m,double);
    if( qtr==NULL )
	ERRMSG("lasso_alloc", no_dyn_mem_message);
    rinvt_theta = Calloc(m,double);
    if( rinvt_theta==NULL )
	ERRMSG("lasso_alloc", no_dyn_mem_message);
    step = Calloc(m,double);
    if( step==NULL )
	ERRMSG("lasso_alloc", no_dyn_mem_message);
    qr_init(n);
}
static void lasso_free(void) {
  num_nz_x=0;
  ytyd2 = 0.0;
  Free(nz_x);
  Free(theta);
  Free(xtr);
  Free(btmp);
  Free(qtr);
  Free(rinvt_theta);
  Free(step);
  qr_free();
}
static void qr_init(int n) {
#if defined (S_or_R)
    if(qr_max_size!=0 || r_ncol!=0 || q_nrow!=0 || q_use_row!=0 ||
       qmat!=NULL || rmat!=NULL) {
	LASSO2_MESSAGE "Possible memory corruption or memory leak.\n  We"
	    "advise to restart your S+ session"  LASSO2_WARNING(LASSO2_NULL_ENTRY);
	qr_free();
    }
#endif
    qr_max_size = 2*n; //QR_CHUNK;
    r_ncol = 0;
    q_nrow = n;
    qmat = Calloc(n*qr_max_size,double);
    if(qmat==NULL)
	ERRMSG("qr_init", no_dyn_mem_message);
    rmat = Calloc(qr_max_size*(qr_max_size+1)/2,double);
    if(rmat==NULL)
	ERRMSG("qr_init", no_dyn_mem_message);
}
static void qr_incr(void) {
  qr_max_size += QR_CHUNK;

  /* reallocate R always */
  rmat = Realloc(rmat,qr_max_size*(qr_max_size+1)/2,double);
  if(rmat==NULL)
    ERRMSG("qr_incr", no_dyn_mem_message);

  /* reallocate Q only if necessary and only to maximal necessary size */
    qmat = Realloc(qmat,q_nrow*qr_max_size,double);
    if(qmat==NULL)
      ERRMSG("qr_incr", no_dyn_mem_message);
}
static void qr_free(void)
{
  qr_max_size = 0;
  r_ncol = 0;
  q_nrow = 0;
  q_use_row = 0;
  Free(qmat);
  Free(rmat);
}
static void qr_del(int l, int aug) {

    double c, s, tau, nu, *col_k, *col_kp1, *a, *b, tmp;
    int i, j, k, l0;

    if( l<0 || l>=r_ncol )
	ERRMSG("qr_del", "Invalid column number");
    r_ncol--;
    if(l == r_ncol) {
	if( aug )
	    ERRMSG("qr_del", "Trying to delete last column of augmented matrix");
	else
	    return;
    }
    /* this is TRUE if $m\le n$ ($m-1\le n$ if the matrix is augmented) */
    if ( r_ncol < q_nrow ) {
	for(k=l; k<r_ncol; k++) /* Update the factorisation and be done */
	{

	    col_k = RCOL(k);  /* first element in column $k$ in R */
	    Memcpy(col_k,col_k+k+1,k+1);
	    a = col_k+k;
	    b = a+k+2;
	    tau = fabs(*a)+fabs(*b);
	    if( tau == 0.0 ) continue; /* both elements are zero
					  nothing to update */
	    nu = tau*sqrt((*a/tau)*(*a/tau)+(*b/tau)*(*b/tau));
	    c = *a/nu;
	    s = *b/nu;
	    *a = nu;
	    b += k+2;
	    a = b-1;
	    for(j=k+2;j<=r_ncol;j++, a+=j, b+=j) {
		tmp = c * *a + s * *b;
		*b  = c * *b - s * *a;
		*a  = tmp;
	    }
	    col_k = QCOL(k);
	    col_kp1 = col_k+q_nrow;
	    for(j=0;j<q_nrow;j++) {
		tmp = c*col_k[j]+s*col_kp1[j];
		col_kp1[j] = c*col_kp1[j]-s*col_k[j];
		col_k[j] = tmp;
	    }
	}
    } else
	if( l < q_nrow ) {
	    /* Update columns upto $m$ and than shift remaining columns */
	    for(k=l; k<q_nrow-1; k++) {

		col_k = RCOL(k);  /* first element in column $k$ in R */
		Memcpy(col_k,col_k+k+1,k+1);
		a = col_k+k;
		b = a+k+2;
		tau = fabs(*a)+fabs(*b);
		if( tau == 0.0 ) continue; /* both elements are zero
					      nothing to update */
		nu = tau*sqrt((*a/tau)*(*a/tau)+(*b/tau)*(*b/tau));
		c = *a/nu;
		s = *b/nu;
		*a = nu;
		b += k+2;
		a = b-1;
		for(j=k+2;j<=r_ncol;j++, a+=j, b+=j) {
		    tmp = c * *a + s * *b;
		    *b  = c * *b - s * *a;
		    *a  = tmp;
		}
		col_k = QCOL(k);
		col_kp1 = col_k+q_nrow;
		for(j=0;j<q_nrow;j++) {
		    tmp = c*col_k[j]+s*col_kp1[j];
		    col_kp1[j] = c*col_kp1[j]-s*col_k[j];
		    col_k[j] = tmp;
		}
	    }
	    l0 = q_nrow-1;

	    col_k   = RCOL(l0); /* first element in column $l_0$ in R */
	    for(i=l0; i<r_ncol; i++, col_k +=i)
		Memcpy(col_k,col_k+i+1,q_nrow);

	    /* [[rmat + (q_nrow-1)*q_nrow/2+q_nrow-1]] is last element in column
	       [[q_nrow]] in R */
	    a = rmat + (q_nrow-1)*(q_nrow+2)/2;
	    if( *a < 0 ) {
		/* [[j<r_ncol]] sufficient since we shifted the columns already */
		for(j=l0;j<r_ncol;j++,a+=j)
		    *a = - *a;
		col_k = qmat+(q_nrow-1)*q_nrow;
		for(i=0;i<q_nrow;i++,col_k++)
		    *col_k = -*col_k;
	    }
	} else { /* just shift last columns to the left */
	    l0 = l;

	    col_k   = RCOL(l0); /* first element in column $l_0$ in R */
	    for(i=l0; i<r_ncol; i++, col_k +=i)
		Memcpy(col_k,col_k+i+1,q_nrow);

	    /* [[rmat + (q_nrow-1)*q_nrow/2+q_nrow-1]] is last element in column
	       [[q_nrow]] in R */
	    a = rmat + (q_nrow-1)*(q_nrow+2)/2;
	    if( *a < 0 ) {
		/* [[j<r_ncol]] sufficient since we shifted the columns already */
		for(j=l0;j<r_ncol;j++,a+=j)
		    *a = - *a;
		col_k = qmat+(q_nrow-1)*q_nrow;
		for(i=0;i<q_nrow;i++,col_k++)
		    *col_k = -*col_k;
	    }
	}
}
static void qr_add(double *x, int swap) {

    double tmp, norm_orig, norm_init, norm_last, *q_new_col;
    int i;
    double *col_l, *col_lm1, *q_elem, *r_new_col;
    double norm_new;
    int j;
    double c, s, tau, nu, *a, *b;

    if( r_ncol == qr_max_size )
	qr_incr();

    norm_orig = 0.0;
    tmp = 0.0;
    for(i=0; i<q_nrow; i++) tmp += x[i]*x[i];
    norm_orig = sqrt(tmp);
    if(r_ncol<q_nrow) {
	q_new_col = QCOL(r_ncol); /* points to the first element of new column */
	tmp = 0.0;
	for(i=0; i<q_nrow; i++) {
	    q_new_col[i] = x[i]/norm_orig;
	    tmp += q_new_col[i]*q_new_col[i];
	}
	if(r_ncol==0) {
	    *rmat = norm_orig;
	    r_ncol++;
	    return;
	}
	norm_init = norm_last = sqrt(tmp);
    }
    r_new_col = RCOL(r_ncol);
    for(i=0; i<=r_ncol; i++) r_new_col[i] = 0.0;
    if( r_ncol >= q_nrow ) {
	/* Multiply new column with Q-transpose in [[qr_add]] */
	q_elem = qmat;
	for(i=0;i<q_nrow;i++) {
	    tmp = 0.0;
	    for(j=0;j<q_nrow;j++,q_elem++)
		tmp += *q_elem * x[j];
	    r_new_col[i] = tmp;
	}
	if( swap ) { /* swap last two columns */
	    col_l   = r_new_col;
	    col_lm1 = r_new_col-r_ncol;
	    for(i=0; i<q_nrow; i++,col_l++,col_lm1++) {
		tmp = *col_l;
		*col_l = *col_lm1;
		*col_lm1 = tmp;
	    }
	    col_lm1--;
	    col_l--;
	    if(q_nrow==r_ncol && *col_lm1<0) {
		/* if new R[$n$,$n$] is negative, than change signs */
		*col_lm1 = - *col_lm1;
		*col_l   = - *col_l;
		q_new_col = qmat+(q_nrow-1)*q_nrow;
		for(j=0;j<q_nrow;j++)
		    q_new_col[j] = -q_new_col[j];
	    }
	}
	r_ncol++;
	return;
    }
    while(1) {

	q_elem = qmat;
	for(j=0; j<r_ncol;j++) {
	    tmp = 0.0;
	    for(i=0;i<q_nrow;i++,q_elem++)
		tmp += *q_elem * q_new_col[i];
	    r_new_col[j] += tmp;
	    q_elem -= q_nrow;
	    for(i=0;i<q_nrow;i++,q_elem++)
		q_new_col[i] -= *q_elem*tmp;
	}
	tmp = 0.0;
	for(i=0;i<q_nrow;i++) tmp += q_new_col[i]*q_new_col[i];
	norm_new = sqrt(tmp);

	if( norm_new >= norm_last/2.0) break;
	if( norm_new > 0.1*norm_init*DBL_EPSILON )
	    norm_last = norm_new;
	else{
	    norm_init = norm_last = 0.1*norm_last*DBL_EPSILON;
	    for(i=0;i<q_nrow;i++)
		q_new_col[i] = 0.0;
	    if(q_use_row == q_nrow)
		ERRMSG("qr_add","Cannot orthogonalise new column");
	    q_new_col[q_use_row] = norm_new = norm_last;
	    q_use_row++;
	}
    }

    for(i=0;i<q_nrow;i++)
	q_new_col[i] /= norm_new;
    for(i=0;i<r_ncol;i++)
	r_new_col[i] *= norm_orig;
    r_new_col[r_ncol] = norm_new*norm_orig;
    if(swap) {

	col_l = r_new_col; /* first element in last column in R */
	col_lm1 = r_new_col-r_ncol; /* first element in column before last in R */
	for(j=0;j<r_ncol;j++,col_l++,col_lm1++) {
	    tmp = *col_l;
	    *col_l = *col_lm1;
	    *col_lm1 = tmp;
	}
	a = col_lm1-1;
	b = col_l;
	tau = fabs(*a)+fabs(*b);
	if( tau > 0.0 ) {
	    /*Calculate Givens rotation*/
	    nu = tau*sqrt((*a/tau)*(*a/tau)+(*b/tau)*(*b/tau));
	    c = *a/nu;
	    s = *b/nu;

	    *a = nu;   /* Fix column before last in R */
	    *b = -s* *(b-1);  /* Fix last element of last column in R */
	    *(b-1) = c * *(b-1); /* Fix second last element of last column in R */
	    if( *b < 0) {
		*b = -*b;

		col_lm1 = QCOL(r_ncol-1);
		col_l   = col_lm1+q_nrow;
		for(j=0;j<q_nrow;j++) {
		    tmp = c*col_lm1[j]+s*col_l[j];
		    col_l[j] = -c*col_l[j]+s*col_lm1[j];
		    col_lm1[j] = tmp;
		}

	    } else {

		col_lm1 = QCOL(r_ncol-1);
		col_l   = col_lm1+q_nrow;
		for(j=0;j<q_nrow;j++) {
		    tmp = c*col_lm1[j]+s*col_l[j];
		    col_l[j] = c*col_l[j]-s*col_lm1[j];
		    col_lm1[j] = tmp;
		}
	    }
	}

    }
    r_ncol++;

}
#if defined (S_or_R)
static void errmsg(char *string) {
  LASSO2_PROBLEM "%s\n", string LASSO2_RECOVER(LASSO2_NULL_ENTRY);
}
#else
static void errmsg(char *where, char *string) {
  fprintf(stderr, "Error in %s: %s\n", where, string);
  exit(EXIT_FAILURE);
}
#endif
