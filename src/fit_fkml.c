/* Evaluate the objective function of the following methods:
 
  1. Maximum Likelihood (ML)
  2. Maximum Spacings Product (MSP)
  3. Titterington's Method (TM)
  4. Starship Method (SM)
  5. Method of TL-moments (TL)
 
Written by Benjamin Dean (21/09/2010)
Changes, Robert King (23/10/2013)

Acknowledgements: The "do_Fx" function is a modified version of gld.fmkl.fx.c, 
which was written by Robert King. This file is included in the "gld" package,
which can be downloaded from http://www.r-project.org. It was released under
the GNU General Public License. */

// TODO: Take inverse.eps as an argument to do_Fx and fit_fkml from R

#include <R.h>
#include <math.h>
#include <stdlib.h>
#include <Rmath.h>

/*##############################################################################
#                             Function prototypes                              #
##############################################################################*/

double do_Fx (const double *, const double *, const double *, const double *, 
  const double *, const double *);

void do_diffn (const double *, const double *, double *, double *, 
  const double *, const double *, const double *, const double *);

double do_QF (const double *,const double *,const double *,const double *,
  const double *);

double sum (const double *, int);

double do_fx (const double *, const double *, const double *, const double *, 
  const double *, const double *);
	
void sort (double *, const int *); 

int sort_compare (const void *, const void *);

int conditions (double *, double *, double *, int *, double *, 
  double *, double *, double *);

double TL_moments (int, double *, double *, double, double *, double *);

double int_term (double, double, double *, double *);

double TL_sample_moments (int, double *, double *, double *, int *);

void compute_l1_l2 (double *, double *, double *, int *, double *, 
  double *, double *, double *);


/*##############################################################################
#                              Objective function                              #
##############################################################################*/

void fit_fkml (int *method, double *l1, double *l2, double *l3, double *l4, 
  double x[], int *n, double *penalty_value, double *t1, double *t2,
  double *l1_tl, double *l2_tl, int *valid_tl, double *obj_result, double *p_to_inverse_eps) {

if (*l2 <= 0) {
				
  *obj_result = *penalty_value;	
	
} else {
		
  /* ================================== ML ===================================*/ 
		
  if (*method == 1) {
			    			
    /* Create an integer to use in the loop */
    int i;
		
    /* Create an array to store log densities */
    double log_density [*n];
			
    for (i = 0; i < *n; i++) {
								
      /* f(x) is computed in two stages.
      1. Given x, compute u = F(x)
      2. Use this u value to compute the density. 
      (This density then corresponds to f(x)) */				
			
      /* 1. Compute u = F(x) */
      double u = do_Fx (&x[i], l1, l2, l3, l4, p_to_inverse_eps);
			
      /* 2. Using this u value, compute the density */
      double density = do_fx (&u, &x[i], l1, l2, l3, l4);
			
      /* Compute log(density) */
      if (density == 0) {
        log_density[i] = -*penalty_value;
      } else {
        log_density[i] =  log(density);
      }
			 				
    }
		
    /* Store value of objective function */
    /* (maximising likelihood = minimising -1*likelihood) */
    *obj_result = -sum (log_density, *n);
			
  }

  /* ================================== MSP ==================================*/ 
				
  if (*method == 2) {      
				
    /* Create an integer to use in the loop */	
    int i;
			
    /* Sort x (rewrites 'x' array) */
    sort (x, n);					
			
    /* Create a vector to store the values of the distribution function */
    /* (The length is n+2) */
    double F [*n+2];
	
    /* Set first element of array to 0, last element of array to 1 */
    /* (since F(support min) = 0 and F(support max) = 1) */
    F[0] = 0;
    F[*n+1] = 1;
	
    /* Compute the distribution function at all other elements */	
    for (i = 1; i <= *n; i++) {
      F[i] = do_Fx (&x[i-1], l1, l2, l3, l4, p_to_inverse_eps);
    }
	
    /* Create a vector to store the log of the pairwise differences */
    /* (There are n+1 pairwise differences) */
    double log_differences_all [*n+1];
		
    /* To compute the differences, we start at 1 and run up to n+1. */		
    for (i = 1; i <= (*n+1); i++) {			
      double difference = F[i] - F[i-1];
      if (difference <= 0) {
		
        if (i != (*n+1)) {
				
          /* Repeated (ordered) observations produce zero differences.
          When this happens, replace the difference Di = F(x[i]) - F(x[i-1]) 
          by the density f(x[i]) (see Cheng and Amin, 1983) */
          
          /* Care is required with the indexing. F[i] - F[i-1] corresponds 
          to F(x[i-1]) - F(x[i-2]). Thus, if the difference is zero, replace
          this term by f(x[i-1]). */			
			
          /* Compute u = F(x[i-1]) */
          double u = do_Fx (&x[i-1], l1, l2, l3, l4, p_to_inverse_eps);
							
          /* Using this u value, compute the density */
          double density = do_fx (&u, &x[i-1], l1, l2, l3, l4);
			
          /* Compute log(density) */
          if (density == 0) {
            log_differences_all[i-1] = -*penalty_value;
          } else {
            log_differences_all[i-1] =  log(density);
          }
			
        } else {
					
          /* i = n+1 */
          
          /* When i = n+1, F[i] - F[i-1] corresponds to F[n+1] - F[n] 
          = 1 - F(x[n-1]). In this case, if the difference is zero, 
          replace this term by f(x[n-1]). */

          /* Compute u = F(x[n-1]) */					
          double u = do_Fx (&x[*n-1], l1, l2, l3, l4, p_to_inverse_eps);
					 
          /* Using this u value, compute the density */
          double density = do_fx (&u, &x[*n-1], l1, l2, l3, l4);
					
          /* Compute log(density) */
          if (density == 0) {
            log_differences_all[i-1] = -*penalty_value;
          } else {
            log_differences_all[i-1] =  log(density);
          }				
        }
				
      } else {
        log_differences_all[i-1] = log(difference);
      }
			
    }
		
    /* Store value of objective function */
    *obj_result = -sum (log_differences_all, (*n+1)) / (*n+1);
				
  }

  /* ================================== TM ===================================*/
		
  if (*method == 3) {    

    /* Create an integer to use in the loop */	
    int i;
			
    /* Sort x (rewrites 'x' array) */
    sort (x, n);			
			
    /* Create a vector to store the values of the distribution function */
    /* (The length is n+1) */
    double F [*n+1];
	
    /* Set first element of array to 0, last element of array to 1 */
    /* (since F(support min) = 0 and F(support max) = 1) */
    F[0] = 0;
    F[*n] = 1;
	
    /* Compute the distribution function at all other elements
    (The ith value of the distribution function (ie, F[i])
    is evaluated at 0.5*(x[i-1] + x[i]).) */	
    for (i = 1; i < *n; i++) {
      double avg_i = 0.5*(x[i-1]+x[i]);
      F[i] = do_Fx (&avg_i, l1, l2, l3, l4, p_to_inverse_eps);
    }
	
    /* Create a vector to store the log of the pairwise differences. */
    /* (There are n pairwise differences) */
    double log_differences_all [*n];
		
    /* To compute the differences, we start at 1 and run up to n. */
    for (i = 1; i <= *n; i++) {
      double difference = F[i] - F[i-1];
      if (difference <= 0) { 
		
        if (i != *n) {
					
          /* F[i] - F[i-1] corresponds to F{(x[i-1]+x[i])/2} - 
          F{(x[i-2]+x[i-1])/2}. Thus, if the difference is zero,
          replace this term by f{(x[i-1]+x[i])/2}. */	
			
          /* Compute u = F{(x[i-1]+x[i])/2} */
          double avg1 = 0.5*(x[i-1]+x[i]);
          double u = do_Fx (&avg1, l1, l2, l3, l4, p_to_inverse_eps);
			
          /* Using this u value, compute the density */
          double density = do_fx (&u, &avg1, l1, l2, l3, l4);
			
          /* Compute log(density) */
          if (density == 0) {
            log_differences_all[i-1] = -*penalty_value;
          } else {
            log_differences_all[i-1] =  log(density);
          }	 	
			
        } else {

          /* i = n */ 
          
          /* When i = n, F[i] - F[i-1] = F[n] - F[n-1]
          = 1 - F{(x[n-2]+x[n-1])/2}. In this case, if the difference
          is zero, replace this term by f{(x[n-2]+x[n-1])/2}. */

          /* Compute u = F{(x[n-2]+x[n-1])/2} */
          double avg2 = 0.5*(x[*n-2]+x[*n-1]);
          double u = do_Fx (&avg2, l1, l2, l3, l4, p_to_inverse_eps);

          /* Using this u value, compute the density */
          double density = do_fx (&u, &avg2, l1, l2, l3, l4);
					
          /* Compute log(density) */
          if (density == 0) {
            log_differences_all[i-1] = -*penalty_value;
          } else {
            log_differences_all[i-1] =  log(density);
          }	 										
        }
											
      } else {
        log_differences_all[i-1] = log(difference);
      }
			
    }
    
    /* Store value of objective function */
    *obj_result = -sum (log_differences_all, *n) / *n;	
  
  }

  /* ================================== SM ===================================*/ 
		
  if (*method == 4) {  
			
    /* Create an integer to use in the loop */	
    int i;
		
    /* Create empty array */
    double AD_sum_part [*n];
			
    /* Sort x (rewrites 'x' array) */
    sort (x, n);					
				
    /* Calculate the summation part of the AD statistic: */
    for (i = 0; i < *n; i++) {
      double u_i = do_Fx (&x[i], l1, l2, l3, l4, p_to_inverse_eps);
      AD_sum_part[i] = (2*i + 1) * log(u_i) + (2 * (*n - i) - 1) * log(1 - u_i);
    } 
	
    /* Complete computation of AD: */
    double AD = -*n - (sum(AD_sum_part, *n) / *n);
		
    /* optim doesn't like Inf, so remove it if it occurs */
    if (AD == INFINITY) {
      AD = *penalty_value;
    }
		
    /* Store value of objective function */
    *obj_result = AD;		
  
  }    				

  /* ================================== TL ===================================*/ 
		
  if (*method == 5) {  
			
  /* Sort x (rewrites x) */
  sort(x, n);

  /* Determine if conditions are satisfied */
  *valid_tl = conditions(l3,l4,x,n,t1,t2,l1_tl,l2_tl);

  /* Check conditions are satisfied */
  if (*valid_tl == 0) { /* not satisfied */
    *obj_result = *penalty_value;
  } else { /* satisfied */
    double tau3 = TL_moments(3,t1,t2,1,l3,l4) / TL_moments(2,t1,t2,1,l3,l4);
    double tau4 = TL_moments(4,t1,t2,1,l3,l4) / TL_moments(2,t1,t2,1,l3,l4);
    double tau3_sample = TL_sample_moments(3,t1,t2,x,n) / 
      TL_sample_moments(2,t1,t2,x,n);
    double tau4_sample = TL_sample_moments(4,t1,t2,x,n) / 
      TL_sample_moments(2,t1,t2,x,n); 
    *obj_result = pow((tau3-tau3_sample),2) + pow((tau4-tau4_sample),2);
  }		
  
  }
}
}

/*##############################################################################
#                        Compute F(x) for the FKML gld                         #
##############################################################################*/

double do_Fx (const double *x, const double *l1, const double *l2, 
  const double *l3, const double *l4, const double *p_to_inverse_eps) {
	
int j = 0;
double df = 0;
double dx = 0;
double dxold = 0;
double f = 0;
double fh = 0;
double fl = 0;
double temp = 0;
double xh = 0;
double xl = 0;
double rts = 0;
double u1 = 0;
double u2 = 1;
int iterations = 500;
double inverse_eps = *p_to_inverse_eps;
double inverse_eps_comp = 1 - inverse_eps;
	
double extreme1 = do_QF (&inverse_eps, l1, l2, l3, l4);	
double extreme2 = do_QF (&inverse_eps_comp, l1, l2, l3, l4);
	
if (*x < extreme1) {
  return 0;
} else {
  if (*x > extreme2) {
    return 1;
  } else {
    if (*l3 <= 0) {
      u1 = inverse_eps;
    }
    if (*l4 <= 0) {
      u2 = 1 - inverse_eps;
    }

    do_diffn (&u1, x, &fl, &df, l1, l2, l3, l4);
    do_diffn (&u2, x, &fh, &df, l1, l2, l3, l4);

    if (fl < 0) {
      xl = u1;
      xh = u2;
    } else {
      xh = u1;
      xl = u2;
    }

    rts = 0.5*(u1+u2);
    dxold = fabs(u2-u1);
    dx = dxold;
    do_diffn (&rts, x, &f, &df, l1, l2, l3, l4);

    for (j = 1; j <= iterations; j++) {
		
      if ((((rts - xh)*df - f)* ( (rts-xl)*df - f) >= 0 ) || 
        ( fabs(2*f) > fabs (dxold*df))) {
		
        dxold = dx;
        dx = 0.5*(xh - xl);
        rts = xl +dx;
        if (xl == rts ) { 
          return rts; 
        }
		
      } else {
		
        dxold = dx;
        dx = f/df;
        temp = rts;
        rts -= dx;
        if (temp == rts) { 
          return rts; 
        }
      }
	
      if (fabs(dx) < inverse_eps) { 
        return rts; 
      }
		
      do_diffn (&rts, x, &f, &df, l1, l2, l3, l4);
      if (f < 0) {
        xl = rts;
      } else { 
        xh = rts;
      }
    }
    return rts; 
  }
}
}

/*##############################################################################
#                         Differentiate the FKML gld                           #
##############################################################################*/

void do_diffn (const double *u, const double *x, double *F, double *dFdu, 
  const double *l1, const double *l2, const double *l3, const double *l4) {

if ( *l3 == 0 ) {
  if ( *l4 == 0 ) {
    /* l3 and l4 zero */
    *F = *l1 + (log(*u) - log(1 - *u))/ *l2 - *x;
    *dFdu = (1 /(*u * (1 - *u))) / *l2 ;
  } else {
    /* l3 zero, l4 nonzero */
    *F = *l1 + ( log(*u) - (( pow((1 - *u), *l4) - 1)/ *l4 ) ) / *l2 - *x;
    *dFdu = (1/ *u + pow((1 - *u), *l4 - 1) ) / *l2 ;
  }
} else {
  if ( *l4 == 0 ) {
    /* l3 nonzero, l4 zero */
    *F = *l1 + ((( pow(*u, *l3) - 1)/ *l3) - log(1 - *u) ) / *l2 - *x;
    *dFdu = (pow(*u, *l3 - 1) + 1/(1 - *u)  ) / *l2;
  } else {
    /* l3 and l4 nonzero */
    *F = ( (pow(*u, *l3) - 1 )/ *l3  - (pow((1 - *u), *l4) - 1 )/ *l4 ) / 
      *l2 + *l1 - *x;
    *dFdu = ( (pow(*u, (*l3 - 1))) + (pow((1 - *u), (*l4 - 1))) ) / *l2 ;
  }
}
}

/*##############################################################################
#               Evaluate the quantile function of the FKML gld                 #
##############################################################################*/

double do_QF (const double *u, const double *l1, const double *l2, 
  const double *l3, const double *l4) {
	
if (*l3 == 0) {
  if (*l4 == 0) {
    return ( *l1 + (log(*u) - log(1 - *u)) / *l2 );
  } else {
    return ( *l1 + (log(*u) - (pow((1 - *u), *l4) - 1)/ *l4)/ *l2 );
  }
} else {
  if (*l4 == 0) {
    return ( *l1 + ((pow(*u, *l3) - 1)/ *l3 - log(1 - *u))/ *l2 );
  } else {
    return ( *l1 + ((pow(*u, *l3) - 1) / *l3 - 
      (pow((1 - *u), *l4) - 1) / *l4 ) / *l2 );
  }
}	
}

/*##############################################################################
#                    Evaluate the density of the FKML gld                      #
##############################################################################*/

double do_fx (const double *u, const double *x, const double *l1, 
  const double *l2, const double *l3, const double *l4) {	
	
/* The density must be zero if the data point lies outside the support of the
distribution. So the first step is to compute the support. */
		  
double u_min;
double u_max;
		  		  
if (*l3 > 0) {
  u_min = 0;
} else {
  u_min = 1e-8;
}
	
if (*l4 > 0) {
  u_max = 1;
} else {
  u_max = 1 - 1e-8;
} 
	
double extreme1 = do_QF (&u_min, l1, l2, l3, l4);	
double extreme2 = do_QF (&u_max, l1, l2, l3, l4);
	
/* Set the density to zero if the data point lies outside the support. */
	
if (*x < extreme1 || *x > extreme2) {		
  return (0);
} else {
		
  /* Before we compute the density, we must check for cases where the
  density equation runs into problems. Ie, when u = 0 and l3 < 1, and 
  when u = 1 and l4 < 1. In these cases, the density equation contains
  a term where zero is raised to a negative power. This produces an error.
		  
  To avoid these situations, we make use of the following limits. 
  f(Q(u)) -> 0 as u -> 0 when l3 < 1
  f(Q(u)) -> 0 as u -> 1 when l4 < 1 */		

  if ( (*u == 0 && *l3 < 1) || (*u == 1 && *l4 < 1) ) {
    return (0);
  } else {				
    /* Compute density */			 
    return (*l2 / ( pow(*u, (*l3 - 1)) + pow((1 - *u), (*l4 - 1)) ) );
  }
}
}

/*##############################################################################
#                         Sum the elements of an array                         #
##############################################################################*/

double sum (const double array[], int n) {
	
int i;
double count = 0;

for (i = 0; i < n; i++) {
  count += array[i];
}
return(count);
}

/*##############################################################################
#                      Sort an array into ascending order                      #
##############################################################################*/

void sort (double array[], const int *n) {
qsort (array, *n, sizeof (double), sort_compare);	
}

int sort_compare (const void *a, const void *b) {
const double *da = (const double *) a;
const double *db = (const double *) b;
return (*da > *db) - (*da < *db);
}

/*##############################################################################
#                        Check the TL-moment conditions                        #
##############################################################################*/

int conditions (double *l3, double *l4, double x[], int *n, double *t1, 
  double *t2, double *l1_tl, double *l2_tl) {
  
int valid = 0; 

/* Check that parameter values, moments and support are valid */

if (*l3 > -(1+ *t1) && *l4 > -(1+ *t2)) { /* l3 and l4 valid */ 
  compute_l1_l2 (l3,l4,x,n,t1,t2,l1_tl,l2_tl); /* Compute l1 and l2 */  
  if (TL_moments(2,t1,t2,*l2_tl,l3,l4) > 0) { /* L2 > 0 */
    double tau3 = TL_moments(3,t1,t2,1,l3,l4) / TL_moments(2,t1,t2,1,l3,l4);
    if (-1 < tau3 && tau3 < 1) { /* tau3 valid */ 
      double tau4 = TL_moments(4,t1,t2,1,l3,l4) / TL_moments(2,t1,t2,1,l3,l4);
      if ((0.25)*(5 * pow(tau3,2) - 1) <= tau4 && tau4 < 1) { /* tau4 valid */
        double u_lower = 0;
        double u_upper = 1;
        if (do_QF(&u_lower,l1_tl,l2_tl,l3,l4) <= x[0] && 
          do_QF(&u_upper,l1_tl,l2_tl,l3,l4) >= x[*n-1]) {
          valid = 1; /* All conditions satisfied */  
        }
      }
    }
  }
}

return(valid);
}

/*##############################################################################
#                            Compute the TL-moments                            #
##############################################################################*/

double TL_moments (int r, double *t1, double *t2, double l2, double *l3,
  double *l4) {

int k;
double a, b;
double sum_total = 0;

for (k = 0; k < r; k++) {
  a = r-k+ *t1-1;
  b = k+ *t2;
  sum_total = sum_total + pow(-1,k) * choose(r-1,k) *
    choose(r+ *t1+ *t2-1, k+ *t2) * int_term(a,b,l3,l4);
}

return( (r+ *t1+ *t2) / (r * l2) * sum_total );
}

/*##############################################################################
#                 Compute the integration term for TL-moments                  #
##############################################################################*/

double int_term (double a, double b, double *l3, double *l4) {

/* This function computes the integration term. Limit cases are used
when l3 and l4 are approximately zero (within a certain tolerance).
This prevents int_term going to infinity. */

double tol = 1e-5;

if (fabs(*l3) < tol) {
  if (fabs(*l4) < tol) { /* l3 and l4 zero */
    return( beta(a+1,b+1)*(digamma(a+1) - digamma(b+1)) );
  } else { /* l3 zero, l4 nonzero */
    return( beta(a+1,b+1)*(digamma(a+1) - digamma(a+b+2)) -
      beta(a+1,*l4+b+1)/ *l4 );
  }
} else {
  if (fabs(*l4) < tol) { /* l3 nonzero, l4 zero */
    return( beta(*l3+a+1,b+1)/ *l3 -
      beta(a+1,b+1)*(digamma(b+1) - digamma(a+b+2)) );
  } else { /* l3 and l4 nonzero */
    return( beta(*l3+a+1,b+1)/ *l3 - beta(a+1,*l4+b+1)/ *l4 );
  }
}
}

/*##############################################################################
#                        Compute the sample TL-moments                         #
##############################################################################*/

double TL_sample_moments (int r, double *t1, double *t2, double x[], int *n) {

int i, k;
double sum_total = 0;

for (i = 0; i < *n; i++) {
  for (k = 0; k < r; k++) {
    sum_total = sum_total + pow(-1,k) * choose(r-1,k) * choose(i,r-k+ *t1-1) *
      choose(*n-i-1,k+ *t2) * x[i];
  }
}

return( sum_total / (r * choose(*n,r+ *t1+ *t2)) );
}

/*##############################################################################
#                       Compute l1 and l2 for TL-moments                       #
##############################################################################*/

void compute_l1_l2 (double *l3, double *l4, double x[], int *n, double *t1, 
  double *t2, double *l1_tl, double *l2_tl) {

/* Compute l2 */
*l2_tl = TL_moments(2,t1,t2,1,l3,l4) / TL_sample_moments(2,t1,t2,x,n); 

double shift;
double tol = 1e-5;

/* Compute appropriate shift for l1 */
if (fabs(*l3) < tol) {
  if (fabs(*l4) < tol) { /* l3 and l4 zero */
    shift = 0;
  } else { /* l3 zero, l4 nonzero */
    shift = -1/(*l2_tl * *l4);
  }
} else {
  if (fabs(*l4) < tol) { /* l3 nonzero, l4 zero */
    shift = 1/(*l2_tl * *l3);
  } else { /* l3 and l4 nonzero */
    shift = 1/(*l2_tl * *l3) - 1/(*l2_tl * *l4);
  }
}

/* Compute l1 */
*l1_tl = shift + TL_sample_moments(1,t1,t2,x,n) - 
  TL_moments(1,t1,t2,*l2_tl,l3,l4);

}
