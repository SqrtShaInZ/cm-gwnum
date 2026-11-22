# cm-gwnum
A fork of [Enge's CM](https://gitlab.inria.fr/enge/cm) that injects George Woltman's `gwnum`.

Feel free to contribute through pull requests and issues. AI coding is allowed, provided enough debug is performed.

## Install instructions
Before compiling, make sure to replace the following headers with the ones used to build gwnum.a:

* [lib/giants.h](lib/giants.h)
* [lib/gwbench.h](lib/gwbench.h)
* [lib/gwcommon.h](lib/gwcommon.h)
* [lib/gwdbldbl.h](lib/gwdbldbl.h)
* [lib/gwini.h](lib/gwini.h)
* [lib/gwnum.h](lib/gwnum.h)
* [lib/gwtables.h](lib/gwtables.h)
* [lib/gwthread.h](lib/gwthread.h)
* [lib/gwutil.h](lib/gwutil.h)

You also need to have a copy or symlink to gwnum.a named libgwnum.a, so `$CC $CFLAGS -lgwnum` is able to detect it.

If you built from the source tree, you will need to run `autoreconf`.

After that, type `./configure && make` to build the package, and (`sudo`) `make install` to install the package at `$PREFIX`.

## Original description

The CM software implements the construction of ring class fields of
imaginary quadratic number fields and of elliptic curves with complex
multiplication via floating point approximations. It consists of libraries
that can be called from within a C program and of executable command
line applications.

For more details, please see the texinfo documentation provided with
the software or available at the website
   http://www.multiprecision.org/cm/

