/*
Nomenclature:
m mpz
w woltman
g giants
*/

#include <gmp.h>
#include "gwnum.h"

#define GW_THRESHOLD 8192

void small_base_exp ( gwhandle *gwdata, mpz_ptr mr, double db,
   mpz_srcptr me, mpz_srcptr mn );
void medium_base_exp ( gwhandle *gwdata, mpz_ptr mr, double db,
   mpz_srcptr me, mpz_srcptr mn );
void large_base_exp ( gwhandle *gwdata, mpz_ptr mr, mpz_srcptr mb,
   mpz_srcptr me, mpz_srcptr mn );
int gw_prp ( gwhandle *gwdata, mpz_srcptr mn );
void gw_powm ( gwhandle *gwdata, mpz_ptr mr, mpz_srcptr mb_inp,
   mpz_srcptr me, mpz_srcptr mn );

/*****************************************************************************/
/*****************************************************************************/
