#ifndef SEGCLUST_CONV_HULL_H
#define SEGCLUST_CONV_HULL_H

#include <gsl/gsl_vector.h>

namespace segclust {

gsl_vector *conv_hull(gsl_vector *J, gsl_vector *pen);

}

#endif
