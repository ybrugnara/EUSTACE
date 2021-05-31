#ifndef SEGCLUST_REPMAT_H
#define SEGCLUST_REPMAT_H

#include <gsl/gsl_matrix.h>

namespace segclust {

gsl_matrix *
repmat(gsl_matrix *, int, int);

}
#endif
