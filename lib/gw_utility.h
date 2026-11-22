/*
Nomenclature:
m mpz
w woltman
g giants
*/

#include <gmp.h>
#include "gwnum.h"

#define GW_THRESHOLD 8192

typedef struct cm_gw_ctx_t cm_gw_ctx_t;

void small_base_exp ( gwhandle *gwdata, mpz_ptr mr, double db,
   mpz_srcptr me, cm_gw_ctx_t *ctx );
void medium_base_exp ( gwhandle *gwdata, mpz_ptr mr, double db,
   mpz_srcptr me, cm_gw_ctx_t *ctx );
void large_base_exp ( gwhandle *gwdata, mpz_ptr mr, mpz_srcptr mb,
   mpz_srcptr me, cm_gw_ctx_t *ctx );
int gw_prp ( cm_gw_ctx_t *ctx );
void gw_powm ( cm_gw_ctx_t *ctx, mpz_ptr mr, mpz_srcptr mb_inp,
   mpz_srcptr me );

/*****************************************************************************/
/*****************************************************************************/
