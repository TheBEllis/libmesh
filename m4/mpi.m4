dnl ---------------------------------------------------------------------------
dnl check for the required MPI library
dnl ---------------------------------------------------------------------------
AC_DEFUN([ACX_MPI], [

dnl if MPIHOME is empty, set it to /usr
AS_IF([test "x$MPIHOME" = "x"], [
      dnl also check $MPI_HOME; it is a common name as well.
      AS_IF([test "x$MPI_HOME" = "x"],
      [
      MPIHOME="/usr"
      ],
      [MPIHOME=$MPI_HOME])
    ])

AC_ARG_WITH([mpi],
            AS_HELP_STRING([--with-mpi=PATH],
                           [Prefix where MPI is installed (MPIHOME)]),
            [MPI="$withval"],
            [
              AS_ECHO(["note: MPI library path not given... trying prefix=$MPIHOME"])
              MPI=$MPIHOME
            ])

AS_IF([test -z "$MPI"], [MPI="/usr"])

MPI_LIBS_PATH="$MPI/lib"
MPI_INCLUDES_PATH="$MPI/include"

dnl Check that the compiler uses the library we specified...

dnl look for LAM or other MPI implementation
AS_IF([test -e $MPI_LIBS_PATH/libmpi.a || test -e $MPI_LIBS_PATH/libmpi.so],
      [
        AS_ECHO(["note: using $MPI_LIBS_PATH/libmpi(.a/.so)"])

        dnl Ensure the compiler finds the library...
        tmpLIBS=$LIBS
        AC_LANG_SAVE
        AC_LANG_CPLUSPLUS

        LIBS="-L$MPI_LIBS_PATH $LIBS"

        dnl look for lam_version_show in liblam.(a/so)
        dnl (this is needed in addition to libmpi.(a/so) for
        dnl LAM MPI
        AC_CHECK_LIB([lam],
                     [lam_show_version],
                     [
                       LIBS="-llam $LIBS"
                       MPI_LIBS="-llam $MPI_LIBS"
                     ],
                     [])

        dnl Quadrics MPI requires the elan library to be included too
        if (nm $MPI_LIBS_PATH/libmpi.* | grep elan > /dev/null); then
          AS_ECHO(["note: MPI found to use Quadrics switch, looking for elan library"])
          AC_CHECK_LIB([elan],
                       [elan_init],
                       [
                         LIBS="-lelan $LIBS"
                         MPI_LIBS="-lelan $MPI_LIBS"
                       ],
                       [AC_MSG_ERROR([Could not find elan library... exiting])])
        fi

        AC_CHECK_LIB([mpi],
                     [MPI_Init],
                     [
                       MPI_LIBS="-lmpi $MPI_LIBS"
                       MPI_LIBS_PATHS="-L$MPI_LIBS_PATH"
                       MPI_IMPL="mpi"
                       AC_MSG_RESULT([Found valid MPI installation...])
                     ],
                     [AC_MSG_RESULT([Could not link in the MPI library...]); enablempi=no])

        AC_LANG_RESTORE
        LIBS=$tmpLIBS
      ])

AS_IF([test -e $MPI_LIBS_PATH/libmpich.a || test -e $MPI_LIBS_PATH/libmpich.so],
      [
        AS_ECHO(["note: using $MPI_LIBS_PATH/libmpich(.a/.so)"])

        dnl Ensure the compiler finds the library...
        tmpLIBS=$LIBS
        AC_LANG_SAVE
        AC_LANG_CPLUSPLUS
        LIBS="-L$MPI_LIBS_PATH $LIBS"

        dnl Myricomm MPICH requires the gm library to be included too
        if (nm $MPI_LIBS_PATH/libmpich.* | grep gm_open > /dev/null); then
          AS_ECHO(["note: MPICH found to use Myricomm's Myrinet, looking for gm library"])

                AS_IF([test "x$GMHOME" = "x"], [GMHOME="/usr"])
                AC_ARG_WITH([gm],
                            AS_HELP_STRING([--with-gm=PATH],
                                           [Prefix where GM is installed (GMHOME)]),
                            [GM="$withval"],
                            [
                              AS_ECHO(["note: GM library path not given... trying prefix=$MPIHOME"])
                              GM=$GMHOME
                            ])

                LIBS="-L$GM/lib $LIBS"
                MPI_LIBS="-L$GM/lib $MPI_LIBS"

                AC_CHECK_LIB([gm],
                             [gm_open],
                             [
                               LIBS="$LIBS -lgm"
                               MPI_LIBS="$MPI_LIBS -lgm"
                             ],
                             [AC_MSG_ERROR( [Could not find gm library... exiting])])
        fi

        dnl look for MPI_Init in libmich.(a/so)
        dnl try adding libmpl if we see it there; some MPICH2 versions
        dnl require it.
        AS_IF([test -e $MPI_LIBS_PATH/libmpl.a || test -e $MPI_LIBS_PATH/libmpl.so],
              [
                LIBS="-L$MPI_LIBS_PATH -lmpl $tmpLIBS"
                AC_CHECK_LIB([mpich],
                             [MPI_Init],
                             [
                               MPI_LIBS="-lmpich -lmpl $MPI_LIBS"
                               MPI_LIBS_PATHS="-L$MPI_LIBS_PATH"
                               MPI_IMPL="mpich"
                               AC_MSG_RESULT([Found valid MPICH installation with libmpl...])
                             ],
                             [AC_MSG_RESULT([Could not link in the MPI library...]); enablempi=no])
              ],
              [
                AC_CHECK_LIB([mpich],
                             [MPI_Init],
                             [
                               MPI_LIBS="-lmpich $MPI_LIBS"
                               MPI_LIBS_PATHS="-L$MPI_LIBS_PATH"
                               MPI_IMPL="mpich"
                               AC_MSG_RESULT([Found valid MPICH installation...])
                             ],
                             [AC_MSG_RESULT([Could not link in the MPI library...]); enablempi=no])
              ])

        AC_LANG_RESTORE
        LIBS=$tmpLIBS
      ])

AS_IF([test "x$MPI_IMPL" != x],
      [
        dnl Ensure the compiler finds the header file...
        AS_IF([test -e $MPI_INCLUDES_PATH/mpi.h],
              [
                AS_ECHO(["note: using $MPI_INCLUDES_PATH/mpi.h"])
                tmpCPPFLAGS=$CPPFLAGS
                AC_LANG_SAVE
                AC_LANG_CPLUSPLUS
                CPPFLAGS="-I$MPI_INCLUDES_PATH $CPPFLAGS"
                AC_CHECK_HEADER([mpi.h],
                                [AC_MSG_RESULT([MPI-headers are working as expected.])],
                                [AC_MSG_RESULT([Could not compile in the MPI headers...]); enablempi=no])
                MPI_INCLUDES_PATHS="-I$MPI_INCLUDES_PATH"
                AC_LANG_RESTORE
                CPPFLAGS=$tmpCPPFLAGS
              ],
              [test -e $MPI_INCLUDES_PATH/mpi/mpi.h],
              [
                MPI_INCLUDES_PATH=$MPI_INCLUDES_PATH/mpi
                AS_ECHO(["note: using $MPI_INCLUDES_PATH/mpi.h"])
                tmpCPPFLAGS=$CPPFLAGS
                AC_LANG_SAVE
                AC_LANG_CPLUSPLUS
                CPPFLAGS="-I$MPI_INCLUDES_PATH $CPPFLAGS"
                AC_CHECK_HEADER([mpi.h],
                                [AC_MSG_RESULT([MPI-headers are working as expected.])],
                                [AC_MSG_RESULT([Could not compile in the MPI headers...]); enablempi=no] )
                MPI_INCLUDES_PATHS="-I$MPI_INCLUDES_PATH"
                AC_LANG_RESTORE
                CPPFLAGS=$tmpCPPFLAGS
              ],
              [
                AC_MSG_RESULT([Could not find MPI header <mpi.h>...])
                enablempi=no
              ])
      ],
      [
        dnl no MPI install found, see if the compiler "natively" supports it by
        dnl attempting to link a test application without any special flags.
        AC_TRY_LINK([@%:@include <mpi.h>],
                    [int np; MPI_Comm_size (MPI_COMM_WORLD, &np);],
                    [
                       MPI_IMPL="built-in"
                       AC_MSG_RESULT( [$CXX Compiler Supports MPI] )
                    ],
                    [AC_MSG_RESULT([$CXX Compiler Does NOT Support MPI...]); enablempi=no])
      ])
      VERSION_MPI
])


AC_DEFUN([VERSION_MPI], [

  AS_IF([test "x$enablempi" != xno],
  [
                CPPFLAGS="-I$MPI_INCLUDES_PATH $CPPFLAGS"
        dnl check that MPI_VERSION and MPI_SUBVERSION are defined.
        AC_TRY_LINK([@%:@include <mpi.h>],  dnl better now?
                    [#ifndef MPI_VERSION
                     #error "MPI-Version undefined."
                    #endif
                    ],
                    [ dnl true: Check which version it is:
                    dnl check that MPI_VERSION is not 1; else disable MPI
                    AC_TRY_LINK([@%:@include <mpi.h>],
                                [
                                #if MPI_VERSION == 1
                                #error "Bad MPI-version"
                                #endif
                                int np; MPI_Comm_size (MPI_COMM_WORLD, &np);
                                ],
                                [ dnl it gives true if MPI_VERSION !=1...
                                 dnl check that MPI_VERSION is not 2; else through a warning
                                 AC_TRY_LINK([@%:@include <mpi.h>],
                                             [
                                             #if MPI_VERSION == 2
                                             #error "deprecated MPI-version"
                                             #endif
                                             int np; MPI_Comm_size (MPI_COMM_WORLD, &np);
                                             ],
                                             [
                                               dnl MPI_VERSION >2: Be happy and configure with MPI.
                                              AC_MSG_RESULT([ The MPI found has version >= 3.0. ]);
                                             ],
                                             [
                                              AC_MSG_WARN([MPI-version is 2.X. This feature is deprecated. ]);
                                             ])
                                dnl set necessary variables etc: Here, all enablempi-cases are in one place.
                                AC_DEFINE(HAVE_MPI, 1, [Flag indicating whether or not MPI is available])
                                ],
                                [
                                 AC_MSG_RESULT(["ERROR: MPI-version is 1.X. Currently we support only MPI version >= 2.0.  Disable MPI now..."]); enablempi=no
                                 mpi_major=`grep "MPI_VERSION" $MPI_INCLUDES_PATH/mpi.h`
                                 dnl AS_ECHO([" FOO: $mpi_major"])
                                ])
                    ],
                    [ dnl MPI_VERSION is not defined.
                     AC_MSG_ERROR(["ERROR: MPI-version seems to be too low: MPI < 1.5. Disable MPI now..."]); enablempi=no
                    ])
  ])
dnl Save variables...
AC_SUBST(MPI)
AC_SUBST(MPI_IMPL)
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_LIBS_PATH)
AC_SUBST(MPI_LIBS_PATHS)
AC_SUBST(MPI_INCLUDES_PATH)
AC_SUBST(MPI_INCLUDES_PATHS)
])
