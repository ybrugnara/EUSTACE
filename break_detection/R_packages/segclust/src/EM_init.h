#ifndef SEGCLUST_EM_INIT_H
#define SEGCLUST_EM_INIT_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace segclust {

gsl_vector * EM_init(gsl_vector *x, gsl_matrix *rupt, long K, long P,bool vh);

}

#endif
