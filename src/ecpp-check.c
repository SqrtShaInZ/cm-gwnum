/*

ecpp-check.c - executable for verifying ECPP certificates

Copyright (C) 2024 Andreas Enge

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

#include "params.h"

int main (int argc, char* argv [])
{
   char *filename;
   int res, res_pari;

   cm_pari_init ();

   evaluate_parameters_ecpp_check (argc, argv, &filename);

   res = cm_pari_ecpp_file_check (filename, true, true);

   if (res == 1)
      printf ("The file contains a valid ECPP certificate "
         "in the CM format.\n");
   else if (res == 2)
      printf ("The file contains a valid partial ECPP certificate "
         "in the CM format.\n");
   else {
      res_pari = cm_pari_ecpp_file_check (filename, false, true);
      printf ("The certificate in the file is invalid according to the "
         "CM specification.\n");
      if (res_pari == 1)
         printf ("But it is valid ");
      else if (res_pari == 2)
         printf ("But it is a valid partial certificate ");
      else
         printf ("It is also invalid ");
      printf ("according to the PARI/GP specification.\n");
   }

   cm_pari_clear ();

   return 0;
}
