/*

tecpp.c - test of ECPP

Copyright (C) 2021, 2022 Andreas Enge

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

#include "cm-impl.h"
#include <sys/stat.h>
#include <unistd.h>

/*****************************************************************************/

static void test_ecpp_class_numbers (void)
{
   unsigned short h [50];
   unsigned short cached_h [32];
   char tmpdir [64], filename [80];

   cm_ecpp_compute_h_chunk (h, 0, 100);
   if (h [0] != 1 || h [1] != 1 || h [2] != 1 || h [3] != 1
       || h [4] != 1 || h [6] != 2) {
      printf ("Incorrect class numbers for discriminants -3, -4, -7, "
         "-8, -11 or -15.\n");
      exit (1);
   }

   snprintf (tmpdir, sizeof (tmpdir), "/tmp/cm-ecpp-test-%lu",
      (unsigned long) getpid ());
   if (mkdir (tmpdir, 0700) != 0
       || !cm_file_write_h (tmpdir, h, 5)
       || !cm_file_read_h (tmpdir, cached_h, 5)
       || memcmp (h, cached_h, sizeof (cached_h)) != 0) {
      printf ("Class-number cache write/read failed.\n");
      exit (1);
   }

   snprintf (filename, sizeof (filename), "%s/cm_h.dat", tmpdir);
   if (remove (filename) != 0 || rmdir (tmpdir) != 0) {
      printf ("Could not remove temporary class-number cache.\n");
      exit (1);
   }
}

/*****************************************************************************/

static void test_ecpp (mpz_srcptr n)
   /* Try to compute an ECPP certificate for n. */
{
   bool res;

   res =    cm_ecpp (n, CM_MODPOLDIR,
      NULL /* filename */,
      NULL /* tmpdir */,
      false /* print */,
      true /* trust */,
      true /* check */,
      0, /* phases */
      32, /* delta */
      false /* verbose */,
      false /* debug */);

   if (!res) {
      printf ("ECPP certificate verification failed for\n");
      mpz_out_str (stdout, 10, n);
      printf ("\n");
      exit (1);
   }
}

/*****************************************************************************/

int main (int argc, char **argv)
{
   mpz_t n;

   test_ecpp_class_numbers ();
   if (argc == 2 && !strcmp (argv [1], "--class-numbers"))
      return 0;

   mpz_init (n);
   mpz_set_ui (n, 10);
   mpz_pow_ui (n, n, 200);
   mpz_nextprime (n, n);

   cm_pari_init ();
   test_ecpp (n);
   cm_pari_clear ();

   mpz_clear (n);

   return 0;
}
