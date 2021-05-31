#ifndef SEGCLUST_SEGCLUST_TYPES_H
#define SEGCLUST_SEGCLUST_TYPES_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct {
  gsl_vector *phi;
  gsl_matrix *rupt;
    
} param_struct;

#endif
