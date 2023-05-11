# Initialize
AC_INIT([slaterlib], [1.0])
AM_PROG_AR
LT_INIT
AC_CONFIG_MACRO_DIR([m4])

# Generate the Makefile
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

AC_LANG(C)

AC_PROG_FC(FC)

PKG_PROG_PKG_CONFIG([])
PKG_INSTALLDIR([])

# Generate config
AC_CONFIG_HEADERS([config.h])

AC_SUBST(CFLAGS)
AC_SUBST(AM_CFLAGS)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_LDFLAGS)

# Generate Makefile from Makefile.in
AC_CONFIG_FILES([Makefile
                 pkgconfig/slaterlib.pc])

# Done
AC_OUTPUT

echo \
"-------------------------------------------------

${PACKAGE_NAME} Version ${PACKAGE_VERSION} 

Prefix ........:  ${prefix}
CC ............:  ${CC}
FC ............:  ${FC}
CPPFLAGS ......:  ${CPPFLAGS}
CFLAGS ........:  ${CFLAGS}
FCFLAGS........:  ${FCFLAGS}
FC ............:  ${FC}
FCLAGS ........:  ${FCFLAGS}
LDFLAGS .......:  ${LDFLAGS}
LIBS ..........:  ${LIBS}


Now type 'make @<:@<target>@:>@'
where the optional <target> is:
  all           - build C library and Fortran binding
  check	        - run unit tests
  install       - install the library 
--------------------------------------------------"
