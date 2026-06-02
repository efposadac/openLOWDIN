#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>

#include "string_f.h"

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* A. First, the interface to the gsl function that calculates the minimum
   of an N-dimensional function, with the knowledge of the function itself
   and its gradient. */

/* This is a type used to communicate with Fortran; func_d is the type of the
   interface to a Fortran subroutine that calculates the function and its
   gradient. */
typedef void (*func_d)(const int*, const double*, double*, const int*, double*);
typedef struct{
  func_d func;
} param_fdf_t;

/* Compute both f and df together. */
static void
my_fdf (const gsl_vector *v, void *params, double *f, gsl_vector *df)
{
  double *x, *gradient, ff2;
  int i, dim, getgrad;
  param_fdf_t * p;

  p = (param_fdf_t *) params;

  dim = v->size;
  x = (double *)malloc(dim*sizeof(double));
  gradient = (double *)malloc(dim*sizeof(double));

  for(i=0; i<dim; i++) x[i] = gsl_vector_get(v, i);
  getgrad = (df == NULL) ? 0 : 1;

  /* printf("Estoy en fdf C\n"); */

  p->func(&dim, x, &ff2, &getgrad, gradient);

  if(f != NULL)
    *f = ff2;
      
  if(df != NULL){
    for(i=0; i<dim; i++)
      gsl_vector_set(df, i, gradient[i]);
  }

  free(x); 
  free(gradient);
}

static double
my_f (const gsl_vector *v, void *params)
{
  double val;

  /* printf("Estoy en f C\n"); */
  my_fdf(v, params, &val, NULL);
  return val;
}

/* The gradient of f, df = (df/dx, df/dy). */
static void
my_df (const gsl_vector *v, void *params, gsl_vector *df)
{
  /* printf("Estoy en df C\n"); */
  my_fdf(v, params, NULL, df);
}

typedef void (*print_f_ptr)(const int*, const int*, const double*, const double*, const double*, const double*);

int FC_FUNC_(low_minimize, LOW_MINIMIZE)
     (const int *method, const int *dim, double *point, const double *step, 
      const double *tolgrad, const double *toldr, const int *maxiter, func_d f, 
      const print_f_ptr write_info, double *minimum)
{
  int iter = 0;
  int status;
  double maxgrad, maxdr, tol;
  int i;
  double * oldpoint;
  double * grad;

  const gsl_multimin_fdfminimizer_type *T = NULL;
  gsl_multimin_fdfminimizer *s;
  gsl_vector *x;
  gsl_vector *absgrad, *absdr;
  gsl_multimin_function_fdf my_func;

  param_fdf_t p;

  p.func = f;

  /* printf("Dimension: %d\n", *dim); */

  oldpoint = (double *) malloc(*dim * sizeof(double));
  grad     = (double *) malloc(*dim * sizeof(double));

  my_func.f = &my_f;
  my_func.df = &my_df;
  my_func.fdf = &my_fdf;
  my_func.n = *dim;
  my_func.params = (void *) &p;

  /* Starting point */
  x = gsl_vector_alloc (*dim);
  for(i=0; i<*dim; i++) gsl_vector_set (x, i, point[i]);

  /* Allocate space for the gradient */
  absgrad = gsl_vector_alloc (*dim);
  absdr = gsl_vector_alloc (*dim);

  tol = 0.1;
  switch(*method){
  case 1: 
    T = gsl_multimin_fdfminimizer_steepest_descent;
    break;
  case 2: 
    T = gsl_multimin_fdfminimizer_conjugate_fr;
    break;
  case 3: 
    T = gsl_multimin_fdfminimizer_conjugate_pr;
    break;
  case 4: 
    T = gsl_multimin_fdfminimizer_vector_bfgs;
    break;
  case 5: 
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    break;
  }

  s = gsl_multimin_fdfminimizer_alloc (T, *dim);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, *step, tol);

  do
    {
      iter++;
      /* printf("Mauricio Rodas\n"); */
      /* printf("Estoy en C iterador: %d\n", iter); */
      for(i=0; i<*dim; i++) oldpoint[i] = point[i];

      /* Iterate */
      status = gsl_multimin_fdfminimizer_iterate (s);
      /* printf("===========================================================\n"); */
      /* Get current minimum, point and gradient */
      *minimum = gsl_multimin_fdfminimizer_minimum(s);
      for(i=0; i<*dim; i++) point[i] = gsl_vector_get(gsl_multimin_fdfminimizer_x(s), i);
      for(i=0; i<*dim; i++) grad[i] = gsl_vector_get(gsl_multimin_fdfminimizer_gradient(s), i);

      /* Compute convergence criteria */
      for(i=0; i<*dim; i++) gsl_vector_set(absdr, i, fabs(point[i]-oldpoint[i]));
      maxdr = gsl_vector_max(absdr);
      for(i=0; i<*dim; i++) gsl_vector_set(absgrad, i, fabs(grad[i]));
      maxgrad = gsl_vector_max(absgrad);

      /* Print information */
      write_info(&iter, dim, minimum, &maxdr, &maxgrad, point);
      
      /* Store infomation for next iteration */
      /* for(i=0; i<*dim; i++) oldpoint[i] = point[i]; */

      if (status)
        break;

      if ( (maxgrad <= *tolgrad) || (maxdr <= *toldr) ) status = GSL_SUCCESS;
      else status = GSL_CONTINUE;
    }
  while (status == GSL_CONTINUE && iter <= *maxiter);

  if(status == GSL_CONTINUE) status = 1025;

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x); gsl_vector_free(absgrad); gsl_vector_free(absdr);

  free(oldpoint);
  free(grad);

  return status;
}




/* /\* B. Second, the interface to the gsl function that calculates the minimum */
/*    of a one-dimensional function, with the knowledge of the function itself, */
/*    but not its gradient. *\/ */

/* /\* This is a type used to communicate with Fortran; func1 is the type of the */
/*    interface to a Fortran subroutine that calculates the function and its */
/*    gradient. *\/ */
/* typedef void (*func1)(const double*, double*); */
/* typedef struct{ */
/*   func1 func; */
/* } param_f1_t; */


/* double fn1(double x, void * params) */
/* { */
/*   param_f1_t * p = (param_f1_t* ) params; */
/*   double fx; */
/*   p->func(&x, &fx); */
/*   return fx; */
/* } */

/* void FC_FUNC_(oct_1dminimize, OCT_1DMINIMIZE)(double *a, double *b, double *m, func1 f, int *status) */
/* { */
/*   int iter = 0; */
/*   int max_iter = 100; */
/*   const gsl_min_fminimizer_type *T; */
/*   gsl_min_fminimizer *s; */
/*   gsl_function F; */
/*   int ierr; */
/*   param_f1_t p; */

/*   p.func = f; */

/*   F.function = &fn1; */
/*   F.params = (void *) &p; */

/*   T = gsl_min_fminimizer_brent; */
/*   s = gsl_min_fminimizer_alloc (T); */

/*   ierr = gsl_min_fminimizer_set (s, &F, *m, *a, *b); */

/*   gsl_set_error_handler_off(); */

/*   do */
/*     { */
/*       iter++; */
/*       *status = gsl_min_fminimizer_iterate (s); */

/*       *m = gsl_min_fminimizer_x_minimum (s); */
/*       *a = gsl_min_fminimizer_x_lower (s); */
/*       *b = gsl_min_fminimizer_x_upper (s); */

/*       *status = gsl_min_test_interval (*a, *b, 0.00001, 0.0); */

/*       /\*if (*status == GSL_SUCCESS) printf ("Converged:\n");*\/ */
/*       /\*printf ("%5d [%.7f, %.7f] %.7f \n", iter, *a, *b,*m);*\/ */
/*     } */
/*   while (*status == GSL_CONTINUE && iter < max_iter); */
/*   gsl_min_fminimizer_free(s); */

/* } */




/* /\* C. Third, the interface to the gsl function that calculates the minimum */
/*    of an N-dimensional function, with the knowledge of the function itself, */
/*    but not its gradient. *\/ */

/* /\* This is a type used to communicate with Fortran; funcn is the type of the */
/*    interface to a Fortran subroutine that calculates the function and its */
/*    gradient. *\/ */
/* typedef void (*funcn)(int*, double*, double*); */
/* typedef struct{ */
/*   funcn func; */
/* } param_fn_t; */

/* double fn(const gsl_vector *v, void * params) */
/* { */
/*   double val; */
/*   double *x; */
/*   int i, dim; */
/*   param_fn_t * p; */

/*   p = (param_fn_t *) params; */
/*   dim = v->size; */
/*   x = (double *)malloc(dim*sizeof(double)); */

/*   for(i=0; i<dim; i++) x[i] = gsl_vector_get(v, i); */
/*   p->func(&dim, x, &val); */

/*   free(x); */
/*   return val; */
/* } */

/* typedef void (*print_f_fn_ptr)(const int*, const int*, const double*, const double*, const double*); */

/* int FC_FUNC_(oct_minimize_direct, OCT_MINIMIZE_DIRECT) */
/*      (const int *method, const int *dim, double *point, const double *step,  */
/*       const double *toldr, const int *maxiter, funcn f,  */
/*       const print_f_fn_ptr write_info, double *minimum) */
/* { */
/*   int iter = 0, status, i; */
/*   double size; */

/*   const gsl_multimin_fminimizer_type *T = NULL; */
/*   gsl_multimin_fminimizer *s = NULL; */
/*   gsl_vector *x, *ss; */
/*   gsl_multimin_function my_func; */

/*   param_fn_t p; */
/*   p.func = f; */

/*   my_func.f = &fn; */
/*   my_func.n = *dim; */
/*   my_func.params = (void *) &p; */

/*   /\* Set the initial vertex size vector *\/ */
/*   ss = gsl_vector_alloc (*dim); */
/*   gsl_vector_set_all (ss, *step); */

/*   /\* Starting point *\/ */
/*   x = gsl_vector_alloc (*dim); */
/*   for(i=0; i<*dim; i++) gsl_vector_set (x, i, point[i]); */

/*   switch(*method){ */
/*   case 6: */
/*     T = gsl_multimin_fminimizer_nmsimplex; */
/*     break; */
/*   } */

/*   s = gsl_multimin_fminimizer_alloc (T, *dim); */
/*   gsl_multimin_fminimizer_set (s, &my_func, x, ss); */

/*   do */
/*     { */
/*       iter++; */
/*       status = gsl_multimin_fminimizer_iterate(s); */

/*       if(status) break; */

/*       *minimum = gsl_multimin_fminimizer_minimum(s); */
/*       for(i=0; i<*dim; i++) point[i] = gsl_vector_get(gsl_multimin_fminimizer_x(s), i); */

/*       size = gsl_multimin_fminimizer_size (s); */
/*       status = gsl_multimin_test_size (size, *toldr); */

/*       write_info(&iter, dim, minimum, &size, point); */

/*     } */
/*   while (status == GSL_CONTINUE && iter < *maxiter); */

/*   if(status == GSL_CONTINUE) status = 1025; */

/*   gsl_vector_free(x);  */
/*   gsl_vector_free(ss); */
/*   gsl_multimin_fminimizer_free(s); */
/*   return status; */
/* } */
