/*

tcm.c - tests for cm

Copyright (C) 2009, 2010, 2011, 2012, 2016, 2021 Andreas Enge

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

#include <string.h>
#include "cm-impl.h"

/*****************************************************************************/

static void test_cuberoot (void)
{
   const unsigned long primes[] = {
      2, 3, 5, 7, 13, 19, 31, 37, 43, 61, 73, 109, 163, 271
   };
   mpz_t p, a, root, check, exponent, test;
   mpzx_t f;
   unsigned int i;
   unsigned long x;
   int expected;

   mpz_init (p);
   mpz_init (a);
   mpz_init (root);
   mpz_init (check);
   mpz_init (exponent);
   mpz_init (test);
   for (i = 0; i < sizeof (primes) / sizeof (primes [0]); i++) {
      mpz_set_ui (p, primes [i]);
      if (primes [i] % 3 == 1) {
         mpz_sub_ui (exponent, p, 1);
         mpz_divexact_ui (exponent, exponent, 3);
      }
      for (x = 0; x < primes [i]; x++) {
         mpz_set_ui (a, x);
         if (x == 0 || primes [i] % 3 != 1)
            expected = 1;
         else {
            mpz_powm (test, a, exponent, p);
            expected = (mpz_cmp_ui (test, 1) == 0);
         }
         if (cm_nt_mpz_cuberoot (root, a, p) != expected) {
            fprintf (stderr, "Cube root existence check failed for p=%lu, a=%lu\n",
               primes [i], x);
            exit (1);
         }
         if (expected) {
            mpz_powm_ui (check, root, 3, p);
            if (mpz_cmp (check, a) != 0) {
               fprintf (stderr, "Cube root verification failed for p=%lu, a=%lu\n",
                  primes [i], x);
               exit (1);
            }
         }
      }
   }

   for (i = 0; i < sizeof (primes) / sizeof (primes [0]); i++) {
      if (primes [i] % 3 == 1) {
         unsigned long a0 = 1, a1 = 2, a2 = 3;
         mpz_set_ui (p, primes [i]);
         mpzx_init (f, 3);
         mpz_set_si (f->coeff [0], -(long) (a0 * a1 * a2));
         mpz_mod (f->coeff [0], f->coeff [0], p);
         mpz_set_ui (f->coeff [1], a0 * a1 + a0 * a2 + a1 * a2);
         mpz_mod (f->coeff [1], f->coeff [1], p);
         mpz_set_si (f->coeff [2], -(long) (a0 + a1 + a2));
         mpz_mod (f->coeff [2], f->coeff [2], p);
         mpz_set_ui (f->coeff [3], 1);
         mpzx_oneroot_split_mod (root, f, p, NULL, false, false);
         mpz_set (check, f->coeff [3]);
         mpz_mul (check, check, root);
         mpz_add (check, check, f->coeff [2]);
         mpz_mul (check, check, root);
         mpz_add (check, check, f->coeff [1]);
         mpz_mul (check, check, root);
         mpz_add (check, check, f->coeff [0]);
         mpz_mod (check, check, p);
         if (mpz_sgn (check) != 0) {
            fprintf (stderr, "Cardano root verification failed for p=%lu\n",
               primes [i]);
            exit (1);
         }
         mpzx_clear (f);
      }
   }
   mpz_clear (p);
   mpz_clear (a);
   mpz_clear (root);
   mpz_clear (check);
   mpz_clear (exponent);
   mpz_clear (test);
}

/*****************************************************************************/

static void test_curve (int_cl_t d, char invariant, bool verbose)
{
   cm_param_t param;
   cm_class_t c;
   mpz_t a, b, x, y, p, n, l, co;
   cm_timer_t clock;

   cm_timer_start (clock);

   if (!cm_param_init (param, d, invariant,
      -1, CM_SUBFIELD_PREFERRED, verbose))
      exit (1);
   mpz_init (a);
   mpz_init (b);
   mpz_init (x);
   mpz_init (y);
   mpz_init (p);
   mpz_init (n);
   mpz_init (l);
   mpz_init (co);

   cm_curve_crypto_param (p, n, l, co, d, 200, verbose);

   /* First test with class polynomials. */
   printf ("d = %"PRIicl", inv = %c; ", d, invariant);
   printf ("class polynomial: ");
   fflush (stdout);
   cm_class_init (c, param, verbose);
   if (!cm_class_compute (c, param, true, false, verbose)) {
      printf ("class polynomial computation failed\n");
      exit (1);
   }
   cm_curve_and_point (a, b, x, y, param, c, p, l, co, CM_MODPOLDIR,
      false, verbose, false);
      /* CM_MODPOLDIR is a preprocessor variable defined by the -D
         parameter of gcc */
   cm_class_clear (c);
   printf ("ok; ");

   /* Then test with a class field tower. */
   printf ("class field tower: ");
   fflush (stdout);
   cm_class_init (c, param, verbose);
   if (!cm_class_compute (c, param, false, true, verbose)) {
      printf ("class field tower computation failed\n");
      exit (1);
   }
   cm_curve_and_point (a, b, x, y, param, c, p, l, co, CM_MODPOLDIR,
      false, verbose, false);
   cm_class_clear (c);
   printf ("ok\n");

   mpz_clear (a);
   mpz_clear (b);
   mpz_clear (x);
   mpz_clear (y);
   mpz_clear (p);
   mpz_clear (n);
   mpz_clear (l);
   mpz_clear (co);

   cm_timer_stop (clock);
   if (verbose)
      printf ("\n--- Elapsed time: %.1f\n", cm_timer_get (clock));
}


/*****************************************************************************/

static void small_test (void)
{
   test_curve ((int_cl_t) (-108715), CM_INVARIANT_DOUBLEETA, false);
}

/*****************************************************************************/

static void big_test (void)
{
   /* Weber: d divisible by 4, not by 32, not by 3 */
   test_curve ((int_cl_t) (-108740), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108712), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108716), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108752), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108724), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108728), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108764), CM_INVARIANT_WEBER, false);
   /* Weber: d divisible by 12, not by 32 */
   test_curve ((int_cl_t) (-108708), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108720), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108732), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108744), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108756), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108780), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108792), CM_INVARIANT_WEBER, false);

   /* d==-16, corresponds to fixed bug */
   test_curve ((int_cl_t) (-16), CM_INVARIANT_J, false);

   test_curve ((int_cl_t) (-108708), CM_INVARIANT_J, false);
   test_curve ((int_cl_t) (-108707), CM_INVARIANT_J, false);
   test_curve ((int_cl_t) (-108703), CM_INVARIANT_J, false);
   test_curve ((int_cl_t) (-108707), CM_INVARIANT_GAMMA2, false);
   test_curve ((int_cl_t) (-108703), CM_INVARIANT_GAMMA2, false);
   test_curve ((int_cl_t) (-108711), CM_INVARIANT_GAMMA3, false);
   test_curve ((int_cl_t) (-108708), CM_INVARIANT_ATKIN, false); /* p=71 */
   test_curve ((int_cl_t) (-108707), CM_INVARIANT_ATKIN, false); /* p=47 */
   test_curve ((int_cl_t) (-108711), CM_INVARIANT_ATKIN, false); /* p=59 */
   test_curve ((int_cl_t) (-58767),  CM_INVARIANT_ATKIN, false); /* p=131 */
   test_curve ((int_cl_t) (-299),    CM_INVARIANT_SIMPLEETA, false);

   test_curve ((int_cl_t) (-1032), CM_INVARIANT_DOUBLEETA, false);
      /* N=3*43, s=e=2, both primes ramified */
   test_curve ((int_cl_t) (-1043), CM_INVARIANT_DOUBLEETA, false);
      /* N=3*7, s=2, e=1 */

   test_curve ((int_cl_t) (-455), CM_INVARIANT_MULTIETA, false);
      /* N=2*5*5, s=e=1, 5 and 7 ramified */
   test_curve ((int_cl_t) (-1919), CM_INVARIANT_MULTIETA, false);
      /* N=2*3*5, s=e=3, complex */
}

/*****************************************************************************/

int main (void)
{
   cm_pari_init ();

   test_cuberoot ();
   small_test ();
   big_test ();
   
   cm_pari_clear ();

   return 0;
}
