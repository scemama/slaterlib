# Initialize
AC_INIT([slaterlib], [1.0])
AM_PROG_AR
LT_INIT
AC_CONFIG_MACRO_DIR([m4])

# Generate the Makefile
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

AC_LANG(C)

AC_PROG_FC(FC)

# Generate config
AC_CONFIG_HEADERS([config.h])

VERSION_MAJOR=`echo ${PACKAGE_VERSION} | cut -d. -f1`
VERSION_MINOR=`echo ${PACKAGE_VERSION} | cut -d. -f2`
VERSION_PATCH=`echo ${PACKAGE_VERSION} | cut -d. -f3 | cut -d- -f1`
AC_DEFINE_UNQUOTED([QMCKL_VERSION_MAJOR], [$VERSION_MAJOR], [major version])
AC_DEFINE_UNQUOTED([QMCKL_VERSION_MINOR], [$VERSION_MINOR], [minor version])
AC_DEFINE_UNQUOTED([QMCKL_VERSION_PATCH], [$VERSION_PATCH], [patch version])
AC_SUBST([VERSION_MAJOR])
AC_SUBST([VERSION_MINOR])
AC_SUBST([VERSION_PATCH])

AC_SUBST(CFLAGS)
AC_SUBST(AM_CFLAGS)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_LDFLAGS)

# Generate Makefile from Makefile.in
AC_CONFIG_FILES([Makefile])

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
