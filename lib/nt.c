/*

nt.c - number theoretic helper functions

Copyright (C) 2009, 2010, 2015, 2021, 2022, 2023 Andreas Enge

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

#include <ecm.h>
#include "cm-impl.h"
#include "gw_utility.h"

#ifndef HAVE_FLINT3
static void tree_gcd (mpz_t *gcd, mpz_srcptr n, mpz_t *m, int no_m);
#endif
void cm_nt_powmod (mpz_ptr res, mpz_srcptr base, mpz_srcptr power, mpz_srcptr p);

/*****************************************************************************/

long int cm_nt_gcd (long int a, long int b)
   /* returns the positive gcd of a and b */

{
   long int r, a_local, b_local;

   if (a == 0)
      return b;
   else if (b == 0)
      return a;
   else
   {
      if (a < 0)
         a_local = -a;
      else
         a_local = a;
      if (b < 0)
         b_local = -b;
      else
         b_local = b;
      r = a_local % b_local;
      while (r > 0)
      {
         a_local = b_local;
         b_local = r;
         r = a_local % b_local;
      }
      return b_local;
   }
}

/*****************************************************************************/

int cm_nt_kronecker (int_cl_t a, int_cl_t b)
   /* Compute the Kronecker symbol (a / b) following Algorithm 1.4.12
      in Cohen93 (by the binary algorithm). */
{
   const int tab [8] = {0, 1, 0, -1, 0, -1, 0, 1};
   int k;
   int_cl_t r;

   /* step 1 */
   if (b == 0) {
      if (a == 1 || a == -1)
         return 1;
      else
         return 0;
   }

   /* step 2*/
   if (a % 2 == 0 && b % 2 == 0)
      return 0;

   while (b % 4 == 0)
      b >>= 2;
   if (b % 2 == 0) {
      b >>= 1;
      k = tab [cm_classgroup_mod (a, (uint_cl_t) 8)];
   }
   else
      k = 1;

   if (b < 0) {
      b = -b;
      if (a < 0)
         k = -k;
   }

   /* step 3 */
   a = cm_classgroup_mod (a, b);

   /* steps 4 to 6 */
   while (a != 0) {
      while (a % 4 == 0)
         a >>= 2;
      if (a % 2 == 0) {
         a >>= 1;
         k *= tab [cm_classgroup_mod (b, (uint_cl_t) 8)];
      }
      if (b > a) {
         r = b - a;
         if (   cm_classgroup_mod (a, (uint_cl_t) 4) == 3
             && cm_classgroup_mod (b, (uint_cl_t) 4) == 3)
            k = -k;
         b = a;
         a = r;
      }
      else
         a -= b;
   }

   if (b > 1)
      return 0;
   else
      return k;
}

/*****************************************************************************/

/*****************************************************************************/

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define psieve_max_primes ((160000))
#define psieve_sieve_span ((1 << 22))
#define psieve_sieve_words ((psieve_sieve_span >> 7))
#define psieve_wheel_size ((3 * 5 * 7 * 11 * 13))
typedef long long ll;
typedef unsigned long long ull;
const ull psieve_mask[64] = {
  0x0000000000000001ull,0x0000000000000002ull,0x0000000000000004ull,0x0000000000000008ull,
  0x0000000000000010ull,0x0000000000000020ull,0x0000000000000040ull,0x0000000000000080ull,
  0x0000000000000100ull,0x0000000000000200ull,0x0000000000000400ull,0x0000000000000800ull,
  0x0000000000001000ull,0x0000000000002000ull,0x0000000000004000ull,0x0000000000008000ull,
  0x0000000000010000ull,0x0000000000020000ull,0x0000000000040000ull,0x0000000000080000ull,
  0x0000000000100000ull,0x0000000000200000ull,0x0000000000400000ull,0x0000000000800000ull,
  0x0000000001000000ull,0x0000000002000000ull,0x0000000004000000ull,0x0000000008000000ull,
  0x0000000010000000ull,0x0000000020000000ull,0x0000000040000000ull,0x0000000080000000ull,
  0x0000000100000000ull,0x0000000200000000ull,0x0000000400000000ull,0x0000000800000000ull,
  0x0000001000000000ull,0x0000002000000000ull,0x0000004000000000ull,0x0000008000000000ull,
  0x0000010000000000ull,0x0000020000000000ull,0x0000040000000000ull,0x0000080000000000ull,
  0x0000100000000000ull,0x0000200000000000ull,0x0000400000000000ull,0x0000800000000000ull,
  0x0001000000000000ull,0x0002000000000000ull,0x0004000000000000ull,0x0008000000000000ull,
  0x0010000000000000ull,0x0020000000000000ull,0x0040000000000000ull,0x0080000000000000ull,
  0x0100000000000000ull,0x0200000000000000ull,0x0400000000000000ull,0x0800000000000000ull,
  0x1000000000000000ull,0x2000000000000000ull,0x4000000000000000ull,0x8000000000000000ull
};
int psieve_primes[psieve_max_primes], psieve_mcnt;
ull psieve_sieve[psieve_sieve_words];
ull psieve_pattern[psieve_wheel_size];
void psieve_mark(ull* s, int o) { s[o >> 6] |= psieve_mask[o & 63]; }
void psieve_unmark(ull* s, int o) { s[o >> 6] &= ~psieve_mask[o & 63]; }
int psieve_test(ull* s, int o) { return !(s[o >> 6] & psieve_mask[o & 63]); }
void psieve_update_sieve(int base) {
    int o = base % psieve_wheel_size;
    o = (o + ((o * 105) & 127) * psieve_wheel_size) >> 7;
    for (int i = 0, k; i < psieve_sieve_words; i += k, o = 0) {
        k = min(psieve_wheel_size - o, psieve_sieve_words - i);
        memcpy(psieve_sieve + i, psieve_pattern + o, sizeof(*psieve_pattern) * k);
    }
    if (base == 0) {
        psieve_sieve[0] |= psieve_mask[0];
        psieve_sieve[0] &= ~(psieve_mask[1] | psieve_mask[2] | psieve_mask[3] | psieve_mask[5] | psieve_mask[6]);
    }
    for (int i = 0; i < psieve_mcnt; ++i) {
        ll j = psieve_primes[i] * psieve_primes[i];
        if (j > base + psieve_sieve_span - 1) break;
        if (j > base) j = (j - base) >> 1;
        else {
            j = psieve_primes[i] - base % psieve_primes[i];
            if (!(j & 1)) j += psieve_primes[i];
            j >>= 1;
        }
        while (j < psieve_sieve_span >> 1) {
            psieve_mark(psieve_sieve, j);
            j += psieve_primes[i];
        }
    }
}
void pre_sieve() {
    for (int i = 0; i < (1048576 >> 7); ++i) psieve_sieve[i] = 0;
    for (int i = 3; i < 1024; i += 2) {
        if (psieve_test(psieve_sieve, i >> 1)) {
            for (int j = (i * i) >> 1; j < 1048576; j += i) psieve_mark(psieve_sieve, j);
        }
    }
    psieve_mcnt = 0;
    for (int i = 8; i < 1048576; ++i) {
        if (psieve_test(psieve_sieve, i)) psieve_primes[psieve_mcnt++] = (i << 1) + 1;
    }
    for (int i = 0; i < psieve_wheel_size; ++i) psieve_pattern[i] = 0;
    for (int i = 1; i < psieve_wheel_size * 64; i += 3) psieve_mark(psieve_pattern, i);
    for (int i = 2; i < psieve_wheel_size * 64; i += 5) psieve_mark(psieve_pattern, i);
    for (int i = 3; i < psieve_wheel_size * 64; i += 7) psieve_mark(psieve_pattern, i);
    for (int i = 5; i < psieve_wheel_size * 64; i += 11) psieve_mark(psieve_pattern, i);
    for (int i = 6; i < psieve_wheel_size * 64; i += 13) psieve_mark(psieve_pattern, i);
}
int segment_sieve(int base, int lim, int pcnt_lim, unsigned long int* all, int* pcnt) {
    psieve_update_sieve(base);
    for (int i = 0, p, u; i < psieve_sieve_words; ++i) {
        ull o = ~psieve_sieve[i];
        while (o) {
            p = __builtin_ctzll(o);
            u = base + (i << 7) + (p << 1) | 1;
            if (*pcnt >= pcnt_lim || u >= lim) return 1;
            all[(*pcnt)++] = u;
            o -= o & ((~o) + 1);
        }
    }
    return 0;
}
void fast_sieve(int lim, int pcnt_lim, unsigned long int* all, int* pcnt) {
    pre_sieve();
    all[(*pcnt)++] = 2;
    for (int base = 0;; base += psieve_sieve_span)
        if (segment_sieve(base, lim, pcnt_lim, all, pcnt) == 1) return;
}
#ifdef WITH_MPI
static unsigned long int psieve_P [1200000];
      /* primes up to 10^7 */
#else
static unsigned long int psieve_P [100000];
      /* primes up to 10^5 */
#endif
const unsigned int psieve_size = sizeof (psieve_P) / sizeof (unsigned long int);

unsigned long int cm_nt_next_prime (const unsigned long int n)
   /* returns the prime following n */

{
   static bool init = true;

   if (init) {
      unsigned int i = 0;
      fast_sieve (INT_MAX, psieve_size, psieve_P, &i);
      init = 0;
   }

   if (n < psieve_P [psieve_size - 1]) {
      /* search in prime list; loop invariants:
         left <= right
         P [0], ..., P [left - 1] <= n
         P [right], ..., P [size - 1] > n      */
      unsigned int left = 0, right = psieve_size - 1, middle;
      while (left != right) {
         middle = (left + right) / 2; /* < right */
         if (psieve_P [middle] > n)
            right = middle; /* becomes smaller */
         else
            left = middle + 1; /* becomes bigger, remains <= right */
      }
      return psieve_P [left];
   }
   else {
      printf ("***** Error: cm_nt_next_prime called with an argument\n"
         "that is too large for the precomputed list.\n");
      exit (1);
   }
}

/*****************************************************************************/

void cm_nt_factor (uint_cl_t d, uint_cl_t *factors, unsigned int *exponents)
   /* Factor d by trial division. The prime factors are stored in factors,
      their multiplicities in exponents, which must provide sufficiently
      much space. The list of prime factors is terminated by 0, so that
      17 entries suffice for the 64 bit type. */

{
   uint_cl_t p, p2;
   int i;

   i = 0;
   p = 2;
   p2 = 4;
   do {
      if (d % p == 0) {
         factors [i] = p;
         exponents [i] = 1;
         d /= p;
         while (d % p == 0) {
            exponents [i]++;
            d /= p;
         }
         i++;
      }
      /* We may wish to use cm_nt_next_prime, but its implementation
         is slow beyond the precomputed table. */
      if (p == 2) {
         p = 3;
         p2 = 9;
      }
      else {
         p2 += 4*(p+1);
         p += 2;
      }
   } while (p2 <= d);

   if (d != 1) {
     /* There is a prime factor left. */
     factors [i] = d;
     exponents [i] = 1;
     i++;
   }

   factors [i] = 0;
}

/*****************************************************************************/

uint_cl_t cm_nt_largest_factor (uint_cl_t n)
   /* Return the largest prime factor of n or 1 if n==1. */
{
   uint_cl_t factors [17];
   unsigned int exponents [17];
   int i;

   cm_nt_factor (n, factors, exponents);
   for (i = 0; factors [i] != 0; i++);

   if (i == 0)
      return 1;
   else
      return factors [i - 1];
}


void cm_nt_powmod (mpz_ptr res, mpz_srcptr base, mpz_srcptr power, mpz_srcptr p)
{
   mpz_powm (res, base, power, p);
}

/*****************************************************************************/

unsigned int cm_nt_mpz_tonelli_generator (mpz_ptr q, mpz_ptr z, mpz_srcptr p)
   /* For p an odd prime, compute and return e such that p-1 = 2^e * q. */
{
   unsigned int e = mpz_scan1 (p, 1);
   mpz_tdiv_q_2exp (q, p, e);
   if (e > 1) {
      for (mpz_set_ui (z, 2ul); mpz_legendre (z, p) != -1;
         mpz_add_ui (z, z, 1ul));
      cm_nt_powmod (z, z, q, p);
   }
   return e;
}

/*****************************************************************************/

void cm_nt_mpz_tonelli_with_generator_unnormalized (mpz_ptr root, mpz_srcptr a,
   mpz_srcptr p, unsigned int e, mpz_srcptr q, mpz_srcptr z)
{
   mpz_t a_local, y, x, b, tmp;
   unsigned long int r, m;

   if (mpz_sgn (a) == 0) {
      mpz_set_ui (root, 0ul);
      return;
   }
   mpz_init_set (a_local, a);
   mpz_mod (a_local, a_local, p);
   mpz_init (y); mpz_init (x); mpz_init (b); mpz_init (tmp);
   if (e == 1) {
      mpz_add_ui (tmp, p, 1ul);
      mpz_tdiv_q_2exp (tmp, tmp, 2ul);
      cm_nt_powmod (x, a_local, tmp, p);
   }
   else {
      mpz_set (y, z); r = e;
      mpz_sub_ui (tmp, q, 1ul);
      mpz_tdiv_q_2exp (tmp, tmp, 1ul);
      cm_nt_powmod (x, a_local, tmp, p);
      mpz_mul (b, x, x); mpz_mod (b, b, p);
      mpz_mul (b, b, a_local); mpz_mod (b, b, p);
      mpz_mul (x, x, a_local); mpz_mod (x, x, p);
      while (mpz_cmp_ui (b, 1) != 0) {
         m = 1;
         mpz_mul (tmp, b, b); mpz_mod (tmp, tmp, p);
         while (m < r && mpz_cmp_ui (tmp, 1) != 0) {
            m++; mpz_mul (tmp, tmp, tmp); mpz_mod (tmp, tmp, p);
         }
         if (m == r) {
            mpz_set_ui (root, 0ul);
            mpz_clear (a_local); mpz_clear (y); mpz_clear (x);
            mpz_clear (b); mpz_clear (tmp);
            return;
         }
         mpz_set (tmp, y);
         for (unsigned long i = 0; i < r - m - 1; i++) {
            mpz_mul (tmp, tmp, tmp); mpz_mod (tmp, tmp, p);
         }
         mpz_mul (y, tmp, tmp); mpz_mod (y, y, p);
         r = m;
         mpz_mul (x, x, tmp); mpz_mod (x, x, p);
         mpz_mul (b, b, y); mpz_mod (b, b, p);
      }
   }
   mpz_set (root, x);
   mpz_clear (a_local); mpz_clear (y); mpz_clear (x);
   mpz_clear (b); mpz_clear (tmp);
}

static void cm_nt_mpz_tonelli_with_generator_internal (mpz_ptr root,
   mpz_srcptr a, mpz_srcptr p, unsigned int e, mpz_srcptr q, mpz_srcptr z,
   gwhandle *gwdata)
   /* Compute a square root of a modulo p by the Tonelli-Shanks algorithm. */
{
   mpz_t a_local, y, x, b, tmp;
   unsigned long int r, m;
   mpz_init_set (a_local, a); mpz_mod (a_local, a_local, p);
   if (mpz_sgn (a_local) == 0) {
      mpz_set_ui (root, 0ul); mpz_clear (a_local); return;
   }
   mpz_init (y); mpz_init (x); mpz_init (b); mpz_init (tmp);
   if (e == 1) {
      mpz_add_ui (tmp, p, 1ul);
      mpz_tdiv_q_2exp (tmp, tmp, 2ul);
      if (gwdata != NULL)
         cm_gw_powm (gwdata, x, a_local, tmp, p);
      else
         cm_nt_powmod (x, a_local, tmp, p);
   }
   else {
      mpz_set (y, z); r = e;
      mpz_sub_ui (tmp, q, 1ul);
      mpz_tdiv_q_2exp (tmp, tmp, 1ul);
      if (gwdata != NULL)
         cm_gw_powm (gwdata, x, a_local, tmp, p);
      else
         cm_nt_powmod (x, a_local, tmp, p);
      mpz_mul (b, x, x); mpz_mod (b, b, p);
      mpz_mul (b, b, a_local); mpz_mod (b, b, p);
      mpz_mul (x, x, a_local); mpz_mod (x, x, p);
      while (mpz_cmp_ui (b, 1) != 0) {
         m = 1;
         mpz_mul (tmp, b, b); mpz_mod (tmp, tmp, p);
         while (m < r && mpz_cmp_ui (tmp, 1) != 0) {
            m++; mpz_mul (tmp, tmp, tmp); mpz_mod (tmp, tmp, p);
         }
         if (m == r) {
            printf ("*** mpz_tonelli called with a = "); mpz_out_str (stdout, 10, a);
            printf (" and p = "); mpz_out_str (stdout, 10, p);
            printf (",\nbut a is not a square modulo p!\n"); exit (1);
         }
         mpz_set (tmp, y);
         for (unsigned long i = 0; i < r - m - 1; i++) {
            mpz_mul (tmp, tmp, tmp); mpz_mod (tmp, tmp, p);
         }
         mpz_mul (y, tmp, tmp); mpz_mod (y, y, p);
         r = m;
         mpz_mul (x, x, tmp); mpz_mod (x, x, p);
         mpz_mul (b, b, y); mpz_mod (b, b, p);
      }
   }
   mpz_set (root, x);
   mpz_clear (a_local); mpz_clear (y); mpz_clear (x);
   mpz_clear (b); mpz_clear (tmp);
}

/*****************************************************************************/

void cm_nt_mpz_tonelli_with_generator (mpz_ptr root, mpz_srcptr a,
   gwhandle *gwdata, unsigned int e, mpz_srcptr q, mpz_srcptr z)
{
   mpz_t p;

   mpz_init (p);
   mpz_mul_2exp (p, q, e);
   mpz_add_ui (p, p, 1);
   cm_nt_mpz_tonelli_with_generator_internal (root, a, p, e, q, z,
      gwdata);
   mpz_clear (p);
}

/*****************************************************************************/

void cm_nt_mpz_tonelli (mpz_ptr root, mpz_srcptr a, mpz_srcptr p)
{
   unsigned int e;
   mpz_t q, z;
   mpz_init (q); mpz_init (z);
   e = cm_nt_mpz_tonelli_generator (q, z, p);
   cm_nt_mpz_tonelli_with_generator_internal (root, a, p, e, q, z, NULL);
   mpz_clear (q); mpz_clear (z);
}

/*****************************************************************************/

bool cm_nt_mpz_cuberoot (mpz_ptr root, mpz_srcptr a, mpz_srcptr p)
   /* Compute a cube root modulo the prime p using the
      Adleman-Manders-Miller algorithm. Return false if a is not a cube. */
{
   mpz_t a_local, t, qpow, exponent, inverse, a_t, a_q;
   mpz_t candidate, test, generator, omega, omega2;
   mpz_t residual, correction, correction_power, power_j;
   mpz_t power_q, logarithm, root_q, root_t, check;
   unsigned long int s, j, d;
   bool success = false;

   if (mpz_cmp_ui (p, 2) < 0)
      return false;

   mpz_inits (a_local, t, qpow, exponent, inverse, a_t, a_q,
      candidate, test, generator, omega, omega2, residual, correction,
      correction_power, power_j, power_q, logarithm, root_q,
      root_t, check, NULL);
   mpz_mod (a_local, a, p);

   if (mpz_cmp_ui (p, 3) == 0) {
      mpz_set (root, a_local);
      success = true;
      goto cleanup;
   }
   if (mpz_fdiv_ui (p, 3) == 2) {
      mpz_mul_2exp (exponent, p, 1);
      mpz_sub_ui (exponent, exponent, 1);
      mpz_divexact_ui (exponent, exponent, 3);
      mpz_powm (root, a_local, exponent, p);
      success = true;
      goto cleanup;
   }
   if (mpz_sgn (a_local) == 0) {
      mpz_set_ui (root, 0);
      success = true;
      goto cleanup;
   }
   if (mpz_fdiv_ui (p, 9) == 4) {
      mpz_mul_2exp (exponent, p, 1);
      mpz_add_ui (exponent, exponent, 1);
      mpz_divexact_ui (exponent, exponent, 9);
      mpz_powm (root, a_local, exponent, p);
      mpz_powm_ui (check, root, 3, p);
      success = (mpz_cmp (check, a_local) == 0);
      goto cleanup;
   }
   if (mpz_fdiv_ui (p, 9) == 7) {
      mpz_add_ui (exponent, p, 2);
      mpz_divexact_ui (exponent, exponent, 9);
      mpz_powm (root, a_local, exponent, p);
      mpz_powm_ui (check, root, 3, p);
      success = (mpz_cmp (check, a_local) == 0);
      goto cleanup;
   }

   /* Write p - 1 = 3^s * t with 3 not dividing t. */
   mpz_sub_ui (t, p, 1);
   for (s = 0; mpz_divisible_ui_p (t, 3); s++)
      mpz_divexact_ui (t, t, 3);
   mpz_ui_pow_ui (qpow, 3, s);

   /* Project a onto the subgroup of order t and take its cube root there. */
   mpz_mod (exponent, qpow, t);
   if (mpz_invert (inverse, exponent, t) == 0)
      goto cleanup;
   mpz_mul (exponent, qpow, inverse);
   mpz_powm (a_t, a_local, exponent, p);
   if (mpz_invert (inverse, a_t, p) == 0)
      goto cleanup;
   mpz_mul (a_q, a_local, inverse);
   mpz_mod (a_q, a_q, p);

   if (mpz_cmp_ui (t, 1) == 0)
      mpz_set_ui (root_t, 1);
   else {
      mpz_set_ui (exponent, 3);
      if (mpz_invert (inverse, exponent, t) == 0)
         goto cleanup;
      mpz_powm (root_t, a_t, inverse, p);
   }

   /* A non-cube generates the full 3-Sylow subgroup after raising to t. */
   mpz_sub_ui (exponent, p, 1);
   mpz_divexact_ui (exponent, exponent, 3);
   for (mpz_set_ui (candidate, 2); ; mpz_add_ui (candidate, candidate, 1)) {
      if (mpz_cmp (candidate, p) >= 0)
         mpz_set_ui (candidate, 2);
      mpz_powm (test, candidate, exponent, p);
      if (mpz_cmp_ui (test, 1) != 0)
         break;
   }
   mpz_powm (generator, candidate, t, p);

   /* Lift the base-3 digits of the discrete logarithm in the 3-Sylow
      subgroup, then divide that logarithm by 3 to obtain the root. */
   mpz_ui_pow_ui (power_q, 3, s - 1);
   mpz_powm (omega, generator, power_q, p);
   mpz_mul (omega2, omega, omega);
   mpz_mod (omega2, omega2, p);
   mpz_set (residual, a_q);
   if (mpz_invert (correction, generator, p) == 0)
      goto cleanup;
   mpz_set_ui (power_j, 1);
   mpz_set_ui (logarithm, 0);

   for (j = 0; j < s; j++) {
      mpz_powm (test, residual, power_q, p);
      if (mpz_cmp_ui (test, 1) == 0)
         d = 0;
      else if (mpz_cmp (test, omega) == 0)
         d = 1;
      else if (mpz_cmp (test, omega2) == 0)
         d = 2;
      else
         goto cleanup;

      if (j == 0 && d != 0)
         goto cleanup;
      mpz_addmul_ui (logarithm, power_j, d);
      mpz_pow_ui (correction_power, correction, d);
      mpz_mul (residual, residual, correction_power);
      mpz_mod (residual, residual, p);
      if (j + 1 < s) {
         mpz_pow_ui (correction, correction, 3);
         mpz_mod (correction, correction, p);
         mpz_mul_ui (power_j, power_j, 3);
         mpz_divexact_ui (power_q, power_q, 3);
      }
   }

   if (!mpz_divisible_ui_p (logarithm, 3))
      goto cleanup;
   mpz_divexact_ui (exponent, logarithm, 3);
   mpz_powm (root_q, generator, exponent, p);
   mpz_mul (root, root_t, root_q);
   mpz_mod (root, root, p);
   mpz_powm_ui (check, root, 3, p);
   success = (mpz_cmp (check, a_local) == 0);

cleanup:
   mpz_clears (a_local, t, qpow, exponent, inverse, a_t, a_q,
      candidate, test, generator, omega, omega2, residual, correction,
      correction_power, power_j, power_q, logarithm, root_q,
      root_t, check, NULL);
   return success;
}

/*****************************************************************************/

void cm_nt_mpz_tonelli_si_with_generator (mpz_ptr root,
   const long int a, gwhandle *gwdata, unsigned int e, mpz_srcptr q,
   mpz_srcptr z)
{
   mpz_t tmp_a;

   mpz_init_set_si (tmp_a, a);
   cm_nt_mpz_tonelli_with_generator (root, tmp_a, gwdata, e, q, z);
   mpz_clear (tmp_a);
}

/*****************************************************************************/

void cm_nt_mpz_tonelli_si (mpz_ptr root, const long int a, mpz_srcptr p)
   /* computes a square root of a modulo p */

{
   mpz_t tmp_a;

   mpz_init_set_si (tmp_a, a);
   cm_nt_mpz_tonelli (root, tmp_a, p);
   mpz_clear (tmp_a);
}

/*****************************************************************************/

static void tree_gcd (mpz_t *gcd, mpz_srcptr n, mpz_t *m, int no_m, bool ecm)
   /* This is the workhorse behind the cm_nt_mpz_tree_gcd function; it
      assumes that the total size of the no_m entries of m is less than
      that of n and directly computes a tree. */
{
   mpz_t **tree;
   mpz_t *ecm_numbers;
   ecm_params *param;
   int *width;
   int levels;
   int i, j;

   /* Compute the height of the subproduct tree of the m in levels. */
   for (i = no_m, levels = 1; i > 1; i = (i+1) / 2, levels++);

   /* Compute bottom-up a subproduct tree with m on the leaves. */
   tree = (mpz_t **) malloc (levels * sizeof (mpz_t *));
   width = (int *) malloc (levels * sizeof (int));
   width [0] = no_m;
   tree [0] = (mpz_t *) malloc (no_m * sizeof (mpz_t));
#ifndef WITH_MPI
     _Pragma("omp parallel for num_threads(8)")
#endif
   for (j = 0; j < no_m; j++)
      mpz_init_set (tree [0][j], m [j]);
   for (i = 1; i < levels; i++) {
      width [i] = (width [i-1] + 1) / 2;
      tree [i] = (mpz_t *) malloc (width [i] * sizeof (mpz_t));
#ifndef WITH_MPI
      _Pragma("omp parallel for num_threads(8)")
#endif
      for (j = 0; j < width [i-1] / 2; j++) {
         mpz_init (tree [i][j]);
         mpz_mul (tree [i][j], tree [i-1][2*j], tree [i-1][2*j+1]);
      }
      if (width [i-1] % 2 != 0) {
         j = width [i-1] / 2;
         mpz_init (tree [i][j]);
         mpz_set (tree [i][j], tree [i-1][2*j]);
      }
   }

   /* Replace the tree tops by n modulo the entry. */
   for (j = 0; j < width [levels-1]; j++)
      mpz_mod (tree [levels-1][j], n, tree [levels-1][j]);

   /* Replace top-down the tree entries by n modulo the entry. */
   for (i = levels - 2; i >= 0; i--) {
#ifndef WITH_MPI
      _Pragma("omp parallel for num_threads(8)")
#endif
      for (j = 0; j < (width [i] / 2) * 2; j++)
         mpz_mod (tree [i][j], tree [i+1][j/2], tree [i][j]);
      if (width [i] % 2 != 0) {
         j = width [i] - 1;
         mpz_set (tree [i][j], tree [i+1][j/2]);
      }
   }
   /* Compute the gcd of n mod m [j] and m [j]. */
   if (ecm) {
#ifndef WITH_MPI
      _Pragma("omp parallel for num_threads(8)")
#endif
      for (j = 0; j < no_m; j++) {
         mpz_gcd (tree [0][j], tree [0][j], m [j]);
         mpz_set (gcd [j], tree [0][j]);
      }
      ecm_numbers = (mpz_t *) malloc (no_m * sizeof (mpz_t));
      param = (ecm_params *) malloc (no_m * sizeof (ecm_params));
#ifndef WITH_MPI
      _Pragma("omp parallel for num_threads(8)")
#endif
      for (j = 0; j < no_m; j++) {
         mpz_init (ecm_numbers [j]);
         mpz_div (ecm_numbers [j], m [j], tree [0][j]);
         mpz_gcd (tree [0][j], ecm_numbers [j], tree [0][j]);
         while (mpz_cmp_ui (tree [0][j], 1) != 0) {
            mpz_div (ecm_numbers [j], ecm_numbers [j], tree [0][j]);
            mpz_gcd (tree [0][j], ecm_numbers [j], tree [0][j]);
         }
      }
      for (i = 0; i < 2; i++) {
#ifndef WITH_MPI
         _Pragma("omp parallel for num_threads(8)")
#endif
         for (j = 0; j < no_m; j++) {
            ecm_init (param [j]);
            double B1 = 10000;
            bool small_n = mpz_sizeinbase (m [0], 2) < 1024;
            if (i == 0) {
               B1 = 50000;
               if (small_n) B1 = 500000;
               param [j]->method = ECM_PM1;
            }
            else {
               if (small_n) B1 = 100000;
               param [j]->method = ECM_ECM;
               param [j]->param = ECM_PARAM_SUYAMA;
            }
            int ret = ecm_factor (tree [0][j], ecm_numbers [j], B1, param [j]);
            if (ret > 0 && mpz_cmp (tree [0][j], ecm_numbers [j]) != 0) {
               mpz_mul (gcd [j], gcd [j], tree [0][j]);
               mpz_div (ecm_numbers [j], ecm_numbers [j], tree [0][j]);
               mpz_gcd (tree [0][j], ecm_numbers [j], tree [0][j]);
               while (mpz_cmp_ui (tree [0][j], 1) != 0) {
                  mpz_div (ecm_numbers [j], ecm_numbers [j], tree [0][j]);
                  mpz_gcd (tree [0][j], ecm_numbers [j], tree [0][j]);
               }
            }
            ecm_clear (param [j]);
         }
         fprintf(stderr, "ECM round %d/5 done\r", i + 1);
      }
      for (j = 0; j < no_m; j++)
         mpz_clear (ecm_numbers [j]);
      free (ecm_numbers);
      free (param);
   }
   else {
#ifndef WITH_MPI
      _Pragma("omp parallel for num_threads(8)")
#endif
      for (j = 0; j < no_m; j++)
         mpz_gcd (gcd [j], tree [0][j], m [j]);
   }

   /* Clear the tree. */
   for (i = 0; i < levels; i++) {
      for (j = 0; j < width [i]; j++)
         mpz_clear (tree [i][j]);
      free (tree [i]);
   }
   free (tree);
   free (width);
}

/*****************************************************************************/

void cm_nt_mpz_tree_gcd (mpz_t *gcd, mpz_srcptr n, mpz_t *m, int no_m, bool ecm)
   /* Given a (large) positive integer n and an array of (smaller) no_m >= 1
      positive integers in m, compute the gcd of n with all the m and return
      them in gcd, which needs to provide sufficient space and initialised
      entries. */
{
#if 1
   const unsigned long int size_n = mpz_sizeinbase (n, 2);
   unsigned long int size_m;
   int offset, no_batch;
   for (offset = 0; offset < no_m; offset += no_batch) {
      /* Treat the entries m [offset], ..., m [offset + no_batch - 1],
         where no_batch is chosen such that their product is at most
         size_n / 2; this keeps the subproduct tree of the m, which
         can require a considerable amount of memory, of a size such
         that reduction of n modulo the product at the root still
         requires some work. */
      for (size_m = mpz_sizeinbase (m [offset], 2), no_batch = 1;
           offset + no_batch < no_m
              && size_m + mpz_sizeinbase (m [offset + no_batch], 2)
                 < size_n / 2;
           size_m += mpz_sizeinbase (m [offset + no_batch], 2), no_batch++);
      tree_gcd (gcd + offset, n, m + offset, no_batch, ecm);
   }
#else
   tree_gcd (gcd, n, m, no_m, ecm);
#endif
}

/*****************************************************************************/
/*****************************************************************************/

bool cm_nt_fget_z (mpz_t out, ftype in)
   /* Tries to round the real value in to an integer value out. The return
      value reflects the success of the operation. */

{
   ftype rounded, diff;
   bool ok;
   mp_exp_t expo;

   finit (rounded, fget_prec (in));
   finit (diff, fget_prec (in));

   fround (rounded, in);
   fsub (diff, in, rounded);
   if (fsgn (diff) == 0 || (- fget_exp (diff)) >= 10) {
      expo = fget_z_exp (out, rounded);
      if (expo > 0)
         ok = false;
      else if (expo < 0) {
         mpz_tdiv_q_2exp (out, out, (unsigned long int) (-expo));
         ok = true;
      }
      else
         ok = true;
   }
   else
      ok = false;

   if (!ok) {
      printf ("***** Error in rounding:\n");
      fout_str (stdout, 10, 0ul, in);
      printf ("\n");
      fout_str (stdout, 10, 0ul, rounded);
      printf ("\n");
   }

   fclear (rounded);
   fclear (diff);

   return ok;
}

/*****************************************************************************/

bool cm_nt_cget_zz (mpz_ptr out1, mpz_ptr out2, ctype in, ctype omega)
   /* Tries to round the complex value to an imaginary-quadratic integer
      out1+out2*omega, where omega is the second element of an integral
      basis (or more generally, a non-real complex number). The return
      value reflects the success of the operation. */
{
   ftype tmp;
   bool ok;

   finit (tmp, cget_prec (in));

   fdiv (tmp, cimagref (in), cimagref (omega));
   ok = cm_nt_fget_z (out2, tmp);

   if (ok) {
      fmul_z (tmp, crealref (omega), out2);
      fsub (tmp, crealref (in), tmp);
      ok = cm_nt_fget_z (out1, tmp);
   }

   fclear (tmp);

   return ok;
}

/*****************************************************************************/
/*****************************************************************************/
