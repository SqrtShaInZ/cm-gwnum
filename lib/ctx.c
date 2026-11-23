/*

ctx.c - cached gwnum contexts for number theoretic operations

This file is part of CM.

CM is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your
option) any later version.

*/

#include "cm-impl.h"
#include "gw_utility.h"

#define ZQ 17586324677945679473ULL

void cm_gw_ctx_init (cm_gw_ctx_t *ctx)
{
   mpz_init (ctx->cached_p);
   ctx->initialized = 0;
}

void cm_gw_ctx_clear (cm_gw_ctx_t *ctx)
{
   if (ctx->initialized) {
      gwdone (&ctx->gwdata);
      ctx->initialized = 0;
   }
   mpz_clear (ctx->cached_p);
}

gwhandle *cm_gw_ctx_get (cm_gw_ctx_t *ctx, mpz_srcptr p)
{
   mpz_t modulus;
   giant giant_modulus;
   size_t string_capacity;
   char *string;
   int giant_length;
   int setup_error;

   if (ctx->initialized && mpz_cmp (ctx->cached_p, p) == 0)
      return &ctx->gwdata;

   if (ctx->initialized) {
      gwdone (&ctx->gwdata);
      ctx->initialized = 0;
   }

   mpz_init (modulus);
   mpz_mul_ui (modulus, p, (unsigned long) ZQ);
   string_capacity = mpz_sizeinbase (modulus, 10) + 2;
   string = (char *) malloc (string_capacity);
   if (string == NULL) {
      fprintf (stderr, "cm_gw_ctx_get: unable to allocate modulus string\n");
      abort ();
   }
   mpz_get_str (string, 10, modulus);
   giant_length = (int) (mpz_sizeinbase (modulus, 2) / 32) + 64;
   giant_modulus = allocgiant (giant_length);
   ctog (string, giant_modulus);

   gwinit (&ctx->gwdata);
   gwset_maxmulbyconst (&ctx->gwdata, 255);
   gwset_larger_fftlen_count (&ctx->gwdata, 0);
   setup_error = gwsetup_general_mod_giant (&ctx->gwdata, giant_modulus);
   free (giant_modulus);
   free (string);
   mpz_clear (modulus);
   if (setup_error != 0) {
      char error_text [256];
      gwerror_text (&ctx->gwdata, setup_error, error_text,
         sizeof (error_text));
      fprintf (stderr, "cm_gw_ctx_get: gwnum setup failed: %s\n",
         error_text);
      gwdone (&ctx->gwdata);
      abort ();
   }

   mpz_set (ctx->cached_p, p);
   ctx->initialized = 1;
   return &ctx->gwdata;
}

void cm_mpz_powm (cm_gw_ctx_t *ctx, mpz_ptr z, mpz_srcptr a,
   mpz_srcptr e, mpz_srcptr p)
{
   if (mpz_sizeinbase (e, 2) < GW_THRESHOLD) {
      mpz_powm (z, a, e, p);
      return;
   }
   cm_gw_ctx_get (ctx, p);
   gw_powm (&ctx->gwdata, z, a, e, p);
}

void cm_gw_powm (gwhandle *gwdata, mpz_ptr z, mpz_srcptr a,
   mpz_srcptr e, mpz_srcptr p)
{
   gw_powm (gwdata, z, a, e, p);
}

int cm_nt_is_prime (cm_gw_ctx_t *ctx, mpz_srcptr n)
{
   if (mpz_sizeinbase (n, 2) < GW_THRESHOLD) {
      int res = 0;
      mpz_t r;
      mpz_t e;
      mpz_init (r);
      mpz_init (e);
      mpz_sub_ui (e, n, 1);
      mpz_set_ui (r, 3);
      mpz_tdiv_q_2exp (e, e, 1);
      mpz_powm (r, r, e, n);
      if (mpz_cmp_ui (r, 1) == 0) res = 1;
      else{
         mpz_add_ui (r, r, 1);
         if (mpz_cmp (r, n) == 0) res = 1;
      }
      return res;
   }
   cm_gw_ctx_get (ctx, n);
   return gw_prp (ctx);
}
