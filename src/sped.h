
#ifndef SPED_SPED_H
#define SPED_SPED_H

#include <R.h>
#include <Rinternals.h>

SEXP descent(SEXP pa, SEXP ma, SEXP args, SEXP genes);
SEXP kinship(SEXP pa, SEXP ma);

#endif /* SPED_SPED_H */

