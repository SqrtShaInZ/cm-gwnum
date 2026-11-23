AM_CFLAGS=-std=c99

include_HEADERS = cm.h cm-arith.h

include_HEADERS += gwnum.h gwcommon.h gwthread.h giants.h gw_utility.h

lib_LTLIBRARIES = libcm.la
if MPI
lib_LTLIBRARIES += libcm_mpi.la
endif

libcm_la_LDFLAGS = -version-info 2:0:1

libcm_la_SOURCES = cm.h cm-arith.h cm-impl.h get_version.c \
   timer.c file.c nt.c ctx.c qdev.c modular.c modpol.c mpzx.c classgroup.c \
   modclass.c curve.c param.c class.c jmodp.c pari.c corn.c \
   ecpp.c flint.c gw_utility.c

if MPI
libcm_mpi_la_SOURCES = $(libcm_la_SOURCES) mpi.c

libcm_mpi_la_CFLAGS = $(libcm_CFLAGS) $(MPI_CFLAGS) -DWITH_MPI

libcm_mpi_la_LDFLAGS = $(libcm_la_LDFLAGS) $(MPI_CLDFLAGS)
endif
