/*

pari.c - functions using pari; for factoring polynomials and for computing
generators of class groups

Copyright (C) 2010, 2015, 2018, 2021, 2022, 2023, 2024 Andreas Enge

This file is part of CM.

CM is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the license, or (at your
option) any later version.

CM is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with CM; see the file COPYING. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/
#define sigjmp_buf jmp_buf
#include <pari/pari.h>
#include "cm-impl.h"

static GEN mpz_get_Z (mpz_srcptr z);
static void Z_get_mpz (mpz_ptr z, GEN x);
static GEN icl_get_Z (int_cl_t z);
static int_cl_t Z_get_icl (GEN x);
static GEN mpzx_get_FpX (mpzx_srcptr f, mpz_srcptr p);
static void FpX_get_mpzx (mpzx_ptr f, GEN x);
static GEN cert_get_V (mpz_t **cert, int depth);
static bool primecertentryisvalid (GEN c, int i, bool type, bool verbose);
static int primecertpartialisvalid (GEN c, bool type, bool verbose);
GEN cm_pari_primecertentryisvalid_worker (GEN i, GEN c, GEN type,
   GEN verbose);

/*****************************************************************************/
/*                                                                           */
/* Conversion functions between PARI, GMP and internal types.                */
/*                                                                           */
/*****************************************************************************/

static GEN mpz_get_Z (mpz_srcptr z)
   /* returns the GEN of type t_INT corresponding to z */

{
   const long l = z->_mp_size;
   const long lz = labs (l) + 2;
   int i;
   GEN x = cgeti (lz);

   x [1] = evalsigne ((l > 0 ? 1 : (l < 0 ? -1 : 0))) | evallgefint (lz);
   for (i = 0; i < labs (l); i++)
      *int_W (x, i) = (z->_mp_d) [i];

   return x;
}

/*****************************************************************************/

static void Z_get_mpz (mpz_ptr z, GEN x)
   /* returns via z the gmp integer corresponding to the t_INT x */

{
   const long l = lgefint (x) - 2;
   int i;

   _mpz_realloc (z, l);
   z->_mp_size = (signe (x) > 0 ? l : -l);
   for (i = 0; i < l; i++)
      (z->_mp_d) [i] = *int_W (x, i);
}

/*****************************************************************************/

static GEN icl_get_Z (int_cl_t z)
   /* returns the GEN of type t_INT corresponding to z; for the time being,
      we use a quick and dirty implementation int_cl_t -> mpz_t -> GEN,
      while it should be easier to move the 64 bits around... */

{
   mpz_t zz;
   GEN x;

   mpz_init (zz);
   cm_classgroup_mpz_set_icl (zz, z);
   x = mpz_get_Z (zz);
   mpz_clear (zz);

   return x;
}

/*****************************************************************************/

static int_cl_t Z_get_icl (GEN x)
   /* returns the int_cl_t correspondong to x, again using a quick and
      dirty implementation */

{
   mpz_t zz;
   int_cl_t i;

   mpz_init (zz);
   Z_get_mpz (zz, x);
   i = cm_classgroup_mpz_get_icl (zz);
   mpz_clear (zz);

   return i;
}

/*****************************************************************************/

static GEN mpzx_get_FpX (mpzx_srcptr f, mpz_srcptr p)
   /* Return a GEN of type t_POL over t_INT corresponding to the
      polynomial f, reduced modulo p. */

{
   int i;
   GEN res;
   mpz_t tmp;

   mpz_init (tmp);

   res = cgetg (f->deg + 3, t_POL);
   setvarn (res, 0);
   if (f->deg == -1)
      setsigne (res, 0);
   else {
      setsigne (res, 1);
      for (i = 0; i <= f->deg; i++) {
         mpz_mod (tmp, f->coeff [i], p);
         gel (res, i+2) = mpz_get_Z (tmp);
      }
   }

   mpz_clear (tmp);

   return res;
}

/*****************************************************************************/

static void FpX_get_mpzx (mpzx_ptr f, GEN x)
   /* Return the lift of x from FpX to a polynomial over the integers
      in f. */
{
   int deg, i;

   deg = (signe (x) == 0 ? -1 : (int) poldegree (x, 0));
   mpzx_set_deg (f, deg);
   for (i = 0; i <= deg; i++)
      Z_get_mpz (f->coeff [i], gel (x, i+2));
}

/*****************************************************************************/

static GEN cert_get_V (mpz_t **cert, int depth)
   /* Given a (partial or complete) ECPP certificate of the given depth
      in cert, return a corresponding PARI object. */
{
   GEN c, ci;
   int i, j;

   c = cgetg (depth + 1, t_VEC);
   for (i = 0; i < depth; i++) {
      ci = cgetg (6, t_VEC);
      gel (c, i+1) = ci;
      for (j = 0; j < 4; j++)
         gel (ci, j+1) = mpz_get_Z (cert [i][j]);
      gel (ci, 5) = mkvec2 (mpz_get_Z (cert [i][4]),
                            mpz_get_Z (cert [i][5]));
   }

   return c;
}

/*****************************************************************************/
/*                                                                           */
/* Various simple functions.                                                 */
/*                                                                           */
/*****************************************************************************/

void cm_pari_init ()
{
   static entree functions_gp []={
      {"_cm_pari_primecertentryisvalid_worker", 0,
       (void*) cm_pari_primecertentryisvalid_worker, 11, "GGGG",
       "worker for primecertentryisvalid",
       NULL, 0, 0, NULL},
      {NULL, 0, NULL, 0, NULL, NULL, NULL, 0, 0, NULL}};

#ifdef HAVE_FLINT
   cm_flint_init ();
#endif
   pari_init_opts (1ul<<23, 0, INIT_JMPm | INIT_DFTm | INIT_noIMTm);
      /* Do not capture SIGSEGV. */
   paristack_setsize (1ul<<23, 1ul<<31);
   if (sizeof (long int) * CHAR_BIT >= 64)
     sd_threadsizemax ("1G", d_SILENT);
   pari_add_module (functions_gp);
   pari_mt_init ();
}

/*****************************************************************************/

void cm_pari_clear ()

{
   pari_close ();
#ifdef HAVE_FLINT
   cm_flint_clear ();
#endif
}

/*****************************************************************************/

void cm_pari_print_library ()
{
   pari_sp av;
   GEN v;

   av = avma;

   v = pari_version ();
   printf ("PARI: include %i.%li.%li, lib %li.%li.%li %s, %lu threads\n",
         PARI_VERSION_CODE >> 16, (PARI_VERSION_CODE >> 8) & 255ul,
         PARI_VERSION_CODE & 255ul,
         itos (gel (v, 1)), itos (gel (v, 2)), itos (gel (v, 3)),
         paricfg_mt_engine,
         pari_mt_nbthreads);

   avma = av;
}

/*****************************************************************************/

bool cm_pari_eval_int (mpz_ptr n, char *e)
   /* If the PARI expression e evaluates to a t_INT, set n to this integer
      and return true; otherwise keep n unchanged and return false. */
{
   pari_sp av;
   GEN z;
   bool ok = true;

   av = avma;

   z = gp_read_str (e);
   if (typ (z) == t_INT)
      Z_get_mpz (n, z);
   else
      ok = false;

   avma = av;

   return ok;
}

/*****************************************************************************/

char* cm_pari_sprintf_hfactor (int_cl_t d)
   /* Return a newly allocated string containing the factorisation of the
      class number of d. */
{
   char *h, *res;
   pari_sp av = avma;

   default0 ("output", "0");
   h = pari_sprintf ("%Ps", factorint (qfbclassno0 (icl_get_Z (d), 0), 0));
   res = (char *) malloc ((strlen (h) + 1) * sizeof (char));
   strcpy (res, h);

   avma = av;

   return res;
}

/*****************************************************************************/
/*                                                                           */
/* Functions for mpzx modulo p relying on PARI.                              */
/*                                                                           */
/*****************************************************************************/

int cm_pari_good_root_of_unity (mpz_ptr zeta, mpz_srcptr p, const int deg)
   /* Compute in zeta a root of unity in F_p with p odd that is suitable
      for finding a root of a totally split polynomial of degree deg > 1;
      its order n is returned. A good choice seems to be n close to deg;
      we decrement it until it divides p-1. */

{
   GEN pp, pm1, factn, e, base, zetap;
   int n;

   pari_sp av = avma;

   pp = mpz_get_Z (p);
   pm1 = subis (pp, 1ul);
   for (n = deg + 1; !dvdiu (pm1, n); n--);

   factn = Z_factor (stoi (n));
   e = diviuexact (pm1, n);
   base = gen_1;
   do {
      base = addis (base, 1l);
      zetap = Fp_pow (base, e, pp);
   }
   while (!equaliu (Fp_order (zetap, factn, pp), n));

   Z_get_mpz (zeta, zetap);

   avma = av;

   return n;
}

/*****************************************************************************/

void cm_pari_mpzx_xplusa_pow_modmod (mpzx_ptr g, unsigned long int a,
   mpz_srcptr e, mpzx_srcptr m, mpz_srcptr p)
   /* Compute g = (X+a)^e modulo m and p. */
{
   GEN pp, ep, fp, mp, gp;
   mpzx_t f;

   pari_sp av = avma;

   pp = mpz_get_Z (p);
   ep = mpz_get_Z (e);
   mpzx_init (f, 1);
   mpz_set_ui (f->coeff [1], 1);
   mpz_set_ui (f->coeff [0], a);
   fp = mpzx_get_FpX (f, p);
   mpzx_clear (f);
   mp = mpzx_get_FpX (m, p);

   gp = FpXQ_pow (fp, ep, mp, pp);

   FpX_get_mpzx (g, gp);

   avma = av;
}

/*****************************************************************************/

void cm_pari_mpzx_gcd_mod (mpzx_ptr h, mpzx_srcptr f, mpzx_srcptr g,
   mpz_srcptr p)
   /* Compute h = gcd (f, g) modulo p without imposing a normalisation of the
      leading coefficient of h. */
{
   GEN pp, fp, gp, hp;

   pari_sp av = avma;

   pp = mpz_get_Z (p);
   fp = mpzx_get_FpX (f, p);
   gp = mpzx_get_FpX (g, p);

   hp = FpX_gcd (fp, gp, pp);

   FpX_get_mpzx (h, hp);

   avma = av;
}

/*****************************************************************************/

void cm_pari_mpzx_divexact_mod (mpzx_ptr h, mpzx_srcptr f, mpzx_srcptr g,
   mpz_srcptr p)
   /* Assuming that g divides f, compute the quotient in h. */
{
   GEN fp, gp, pp, hp;

   pari_sp av = avma;

   fp = mpzx_get_FpX (f, p);
   gp = mpzx_get_FpX (g, p);
   pp = mpz_get_Z (p);

   hp = FpX_div (fp, gp, pp);

   FpX_get_mpzx (h, hp);

   avma = av;
}

/*****************************************************************************/
/*                                                                           */
/* Functions for finding roots of polynomials.                               */
/*                                                                           */
/*****************************************************************************/

mpz_t* cm_pari_find_roots (int *no, mpzx_srcptr f, mpz_srcptr p)
   /* Computes all the roots (without multiplicities) of the polynomial f
      modulo p. The number of found roots is returned in no. */

{
   pari_sp av;
   mpz_t *res;
   GEN fp, pp, rootsp;
   int i;

   av = avma;

   pp = mpz_get_Z (p);
   fp = mpzx_get_FpX (f, p);
   rootsp = FpX_roots (fp, pp);
   *no = lg (rootsp) - 1;
   res = (mpz_t*) malloc ((*no) * sizeof (mpz_t));
   for (i = 0; i < *no; i++) {
      mpz_init (res [i]);
      Z_get_mpz (res [i], gel (rootsp, i+1));
   }

   avma = av;

   return res;
}

/*****************************************************************************/
/*                                                                           */
/* Functions for computing class groups.                                     */
/*                                                                           */
/*****************************************************************************/

int cm_pari_classgroup (int_cl_t disc, int_cl_t *ord, cm_form_t *gen)
   /* Given a negative discriminant, compute its ideal class group as a
      product of cyclic groups with their orders and generators. The orders
      are returned in ord and the generators in gen; their number is the
      return value. ord and gen must have been initialised with sufficient
      space to hold the result. */

{
   pari_sp av;
   int length;
   GEN d, cl, orders, gens, qfb;
   int i;

   av = avma;

   d = icl_get_Z (disc);
   cl = quadclassunit0 (d, 0, NULL, 0);
   orders = gel (cl, 2);
   gens = gel (cl, 3);
   length = glength (orders);
   for (i = 0; i < length; i++) {
      ord [i] = Z_get_icl (gel (orders, i+1));
      qfb = gel (gens, i+1);
      gen [i].a = Z_get_icl (gel (qfb, 1));
      gen [i].b = Z_get_icl (gel (qfb, 2));
   }

   avma = av;

   return length;
}

/*****************************************************************************/

void cm_pari_oneroot (mpz_ptr root, mpzx_srcptr f, mpz_srcptr p)
   /* Find a root of the polynomial f over the prime field of
      characteristic p, assuming that f splits completely, and return it
      in the variable of the same name. */
{
   GEN fp, pp, rootp;

   pari_sp av = avma;

   pp = mpz_get_Z (p);
   fp = mpzx_get_FpX (f, p);

   rootp = FpX_oneroot_split (fp, pp);
   Z_get_mpz (root, rootp);

   avma = av;
}


int cm_pari_classgroup_2quotient (int_cl_t disc, const int *p,
   int_cl_t *ord, cm_form_t *gen)
   /* Given a negative discriminant disc and a 0-terminated list of primes
      p [0], p [1], ..., p [n-1] dividing the fundamental part of disc
      with n>=2, compute the quotient of the ideal class group by the
      subgroup of order 2^(n-1) generated by the forms of norm
      p [0] * p [1], ..., p [0] * p [n-1], as a product of cyclic groups
      with their orders and generators. The orders are returned in ord and
      the generators in gen; their number is the return value. ord and gen
      must have been initialised with sufficient space to hold the result.
      The function could be combined with the previous one, but is kept
      apart for the sake of the clarity of the previous function. */

{
   pari_sp av;
   int length;
   GEN d, cl, orders, gens, qfb;
   int length2, size2;
   GEN gens2, group2, dlog2, halforder;
   GEN M, Uinv, ordersq, gensq;
   int lengthq;
   cm_form_t p0pi;
   int i, j;

   av = avma;

   /* Compute the form class group. */
   d = icl_get_Z (disc);
   cl = quadclassunit0 (d, 0, NULL, 0);
   orders = gel (cl, 2);
   gens = gel (cl, 3);
   length = glength (orders);

   /* Compute the size of the 2-elementary subgroup and its generators
      as quadratic forms. */
   for (length2 = 0;
        length2 < length && !mpodd (gel (orders, length2 + 1));
        length2++);
   gens2 = cgetg (length2 + 1, t_VEC);
   for (i = 1; i <= length2; i++)
      gel (gens2, i) = gpowgs (gel (gens, i), itou (gel (orders, i)) / 2);

   /* Enumerate the subgroup elements and their discrete logarithms. */
   size2 = 1 << length2;
   group2 = cgetg (size2 + 1, t_VEC);
   dlog2 = cgetg (size2 + 1, t_VEC);
   size2 = 1;
#if PARI_VERSION_CODE < PARI_VERSION (2, 14, 0)
   gel (group2, 1) = qfi_1 (gel (gens, 1));
#else
   gel (group2, 1) = qfb_1 (gel (gens, 1));
#endif
   gel (dlog2, 1) = zerocol (length);
   for (i = 1; i <= length2; i++) {
      halforder = shifti (gel (orders, i), -1);
      for (j = 1; j <= size2; j++) {
#if PARI_VERSION_CODE < PARI_VERSION (2, 14, 0)
         gel (group2, size2 + j)
            = qficomp (gel (group2, j), gel (gens2, i));
#else
         gel (group2, size2 + j)
            = qfbcomp (gel (group2, j), gel (gens2, i));
#endif
         gel (dlog2, size2 + j) = shallowcopy (gel (dlog2, j));
         gmael (dlog2, size2 + j, i) = halforder;
      }
      size2 <<= 1;
   }

   /* Create a matrix of relations for the quotient group. */
   M = zeromat (length, 0);
   if (disc % p [0] != 0) {
      printf ("***** Error: Calling cm_pari_classgroup_2quotient with "
              "non-ramified prime %i.\n", p [0]);
      exit (1);
   }
   for (i = 1; p [i] != 0; i++) {
      if (disc % p [i] != 0) {
         printf ("***** Error: Calling cm_pari_classgroup_2quotient with "
                 "non-ramified prime %i.\n", p [i]);
         exit (1);
      }
      /* Compute the form of norm p [0] * p [i]. */
      p0pi.a = p [0] * p [i];
      if (disc % 2 != 0)
         p0pi.b = p0pi.a;
      else if (disc % 8 == 0 || p0pi.a % 2 != 0)
         p0pi.b = 0;
      else
         p0pi.b = p0pi.a;
      cm_classgroup_reduce (&p0pi, disc);
#if PARI_VERSION_CODE < PARI_VERSION (2, 14, 0)
      qfb = qfi (icl_get_Z (p0pi.a), icl_get_Z (p0pi.b),
         icl_get_Z (cm_classgroup_compute_c (p0pi.a, p0pi.b, disc)));
#else
      qfb = mkqfb (icl_get_Z (p0pi.a), icl_get_Z (p0pi.b),
         icl_get_Z (cm_classgroup_compute_c (p0pi.a, p0pi.b, disc)),
         icl_get_Z (disc));
#endif
      /* Look up its discrete logarithm. */
      for (j = 1; !gequal (qfb, gel (group2, j)); j++);
      M = shallowconcat (M, gel (dlog2, j));
   }

   /* Compute the quotient group. */
   M = hnfmodid (M, orders);
   ordersq = ZM_snf_group (M, NULL, &Uinv);
   lengthq = glength (ordersq);
   gensq = cgetg (lengthq + 1, t_VEC);
   for (i = 1; i <= lengthq; i++)
      gel (gensq, i) = factorback2 (gens, gel (Uinv, i));

   for (i = 0; i < lengthq; i++) {
      ord [i] = Z_get_icl (gel (ordersq, i+1));
      qfb = gel (gensq, i+1);
      gen [i].a = Z_get_icl (gel (qfb, 1));
      gen [i].b = Z_get_icl (gel (qfb, 2));
   }

   avma = av;

   return lengthq;
}

/*****************************************************************************/
/*                                                                           */
/* Functions used for ECPP                                                   */
/*                                                                           */
/*****************************************************************************/

void cm_pari_prime_product (mpz_ptr prim, unsigned long int a,
   unsigned long int b)
   /* Return in prim the product of all primes p with a < p <= b.
      This function requires some memory to compute the prime list by
      a call to PARI. When calling the function, the PARI stack must be
      empty, so that it can be reset safely. */
{
   pari_sp av;
   GEN p;
   mpz_t *prod;
   int len, len2, i;

   av = avma;
   p = primes_interval_zv (a+1, b);
   len = glength (p);
   if (len == 0) {
      /* This probably never occurs... */
      len = 1;
      prod = (mpz_t *) malloc (1 * sizeof (mpz_t));
      mpz_init_set_ui (prod [0], 1);
   }
   else {
      prod = (mpz_t *) malloc (len * sizeof (mpz_t));
      for (i = 0; i < len; i++)
         mpz_init_set_si (prod [i], p [i+1]);
   }
   avma = av;
   /* The computations usually enlarge the PARI stack, but later steps
      do not need as much PARI memory, so it can be released. */
   parivstack_reset ();

   while (len > 1) {
      len2 = len / 2;
      for (i = 0; i < len2; i++)
         mpz_mul (prod [i], prod [2*i], prod [2*i+1]);
      if (len % 2 != 0) {
         mpz_set (prod [len2], prod [2*len2]);
         len2++;
      }
      for (i = len2; i < len; i++)
         mpz_clear (prod [i]);
      len = len2;
   }

   mpz_set (prim, prod [0]);
   mpz_clear (prod [0]);
   free (prod);
}

/*****************************************************************************/

/*****************************************************************************/

/*****************************************************************************/

static bool primecertentryisvalid (GEN c, int i, bool type, bool verbose)
   /* Check whether line i of the partial ECPP certificate c is valid
      (where numbering starts with 1 as conventional in PARI/GP).
      c is supposed to be of type t_VEC.
      If type is true, assume a CM style certificate, that is, the points
         are of prime order;
      if type is false, assume a PARI/GP style certificate, that is,
         the points must first be multiplied by the cofactor to reach prime
         order. */
{
   GEN ci, N, t, s, a, P, x, y, cip1, N1, m, q, PP, Q;
   int depth;
   bool res;
   pari_sp av;

   av = avma;

   /* Extract the line and check all types in passing. */
   depth = glength (c);
   if (i > depth)
      return false;
   ci = gel (c, i);
   if (glength (ci) != 5)
      return false;
   N = gel (ci, 1);
   t = gel (ci, 2);
   s = gel (ci, 3);
   a = gel (ci, 4);
   P = gel (ci, 5);
   if (typ (N) != t_INT || typ (t) != t_INT || typ (s) != t_INT
       || typ (a) != t_INT || typ (P) != t_VEC || glength (P) != 2)
      return false;
   x = gel (P, 1);
   y = gel (P, 2);
   if ( typ (x) != t_INT || typ (y) != t_INT)
      return false;

   /* Do the same for the first entry of the next line if applicable. */
   if (i < depth) {
      cip1 = gel (c, i+1);
      if (typ (cip1) != t_VEC || glength (cip1) != 5)
         return false;
      N1 = gel (cip1, 1);
      if (typ (N1) != t_INT)
         return false;
   }

   /* Check for numerical boundaries.
      From here on, a simple return is not possible any more, we need
      garbage collection. */
   res = true;
   if (cmpui (2, N) > 0
      || cmpii (sqri (t), mului (4, N)) > 0
      || cmpui (0, s) >= 0
      || cmpui (0, a) > 0 || cmpii (a, N) >= 0
      || cmpui (0, x) > 0 || cmpii (x, N) >= 0
      || cmpui (0, y) > 0 || cmpii (y, N) >= 0)
      res = false;
   else {
      m = subii (addui (1, N), t);
      if (!equalui (0, modii (m, s)))
         res = false;
      else {
         /* Check for consistency with next line if applicable. */
         q = diviiexact (m, s);
         if (i < depth && !equalii (q, N1))
            res = false;
         else {
            /* Check that q > (\sqrt [4]{N} + 1)^2; the following
               equivalent formula has been obtained by repeated sorting
               and squaring. */
            if (cmpii (sqri (subii (sqri (subiu (q, 1)), N)),
                       mului (16, mulii (q, N))) < 0)
               res = false;
            else {
               PP = cgetg (4, t_VEC);
               gel (PP, 1) = x;
               gel (PP, 2) = y;
               gel (PP, 3) = gen_1;
               /* If type is false, replace PP by s*PP and check that it
                  is different from the point at infinity modulo all prime
                  divisors of N, otherwise said that its third projective
                  coordinate is coprime to N. */
               if (!type) {
                  PP = FpJ_mul (PP, s, a, N);
                  if (!equalui (1, gcdii (gel (PP, 3), N)))
                     res = false;
               }
               /* Check that q*PP = O. */
               if (res) {
                  Q = FpJ_mul (PP, q, a, N);
                  res = (equalui (0, modii (gel (Q, 3), N)) != 0);
               }
            }
         }
      }
   }

   if (verbose)
      printf ("  ECPP line %5i: %i\n", i, res);

   avma = av;

   return res;
}

/*****************************************************************************/

GEN cm_pari_primecertentryisvalid_worker (GEN i, GEN c, GEN type,
   GEN verbose)
   /* A version of the previous function adapted to use in parfor. */
{
   int il = itos (i);
   bool typel = (itos (type) == 1);
   bool verbosel = (itos (verbose) == 1);

   bool res = primecertentryisvalid (c, il, typel, verbosel);

   return stoi ((int) res);
}

/*****************************************************************************/

static int primecertpartialisvalid (GEN c, bool type, bool verbose)
   /* Check whether c is a valid partial ECPP certificate.
      type indicates whether the certificate is of type created by CM
      or not.
      The return value is 0, 1 or 2.
      0 means the certificate is invalid.
      2 means the certificate is valid, but only partial.
      1 means the certificate is valid and its first entry is a
        proven prime. */
{
   int depth, val, res;
   GEN ci, N, t, s, q;
   pari_sp av;

   av = avma;

   if (typ (c) != t_VEC)
      return 0;

   depth = glength (c);
   val = 0;
   {
      /* Call primecertentryisvalid in parallel with early abort as soon
         as an entry is invalid. */
      parfor_t iter;
      GEN p2;
      parfor_init (&iter, gen_1, stoi (depth),
         strtoclosure("_cm_pari_primecertentryisvalid_worker", 3,
            c, stoi (type == true), stoi (verbose == true)));

      while ((p2 = parfor_next (&iter))) {
         if (equalui (0, gel (p2, 2))) {
            parfor_stop (&iter);
            break;
         }
         else
            val++;
      }
    }

   if (val != depth)
      return 0;

   ci = gel (c, depth);
   N = gel (ci, 1);
   t = gel (ci, 2);
   s = gel (ci, 3);
   q = diviiexact (subii (addui (1, N), t), s);
   if (expi (q) >= 64)
      res = 2;
   else if (equalui (1, gispseudoprime (N, 0)))
      res = 1;
   else
      res = 0;

   avma = av;

   return res;
}

/*****************************************************************************/

bool cm_pari_ecpp_check (mpz_t **cert, int depth, bool debug)
   /* Given a complete ECPP certificate (after step 2) of length depth in
      cert, use PARI to check it and return its validity.
      If debug is true, some progress information is printed. */
{
   pari_sp av;
   bool res;
   GEN c;

   av = avma;

   c = cert_get_V (cert, depth);

   res = (primecertpartialisvalid (c, true, debug) == 1);

   avma = av;

   return res;
}

/*****************************************************************************/

int cm_pari_ecpp_file_check (const char* filename, bool type, bool verbose)
   /* Given the name of a file containing a (partial) ECPP certificate
      in PARI format, call primecertpartialisvalid on the file content.
      The values of type and verbose are passed through. */
{
   pari_sp av;
   GEN c;
   int res;

   av = avma;

   c = gp_read_file (filename);

   res = primecertpartialisvalid (c, type, verbose);

   avma = av;

   return res;
}

/*****************************************************************************/
