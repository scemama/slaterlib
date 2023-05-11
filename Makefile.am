AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 --install

VERSION_MAJOR   = @VERSION_MAJOR@
VERSION_MINOR   = @VERSION_MINOR@
VERSION_PATCH   = @VERSION_PATCH@

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = pkgconfig/slaterlib.pc

src_main_c = src/core.c src/debug.c src/simple.c src/test.c
src_main_h = src/slater_condon.h src/slater_condon_simple.h src/debug.h

lib_LTLIBRARIES = lib/libslater.la
include_HEADERS = $(src_main_h)
lib_libslater_la_SOURCES = $(src_main_c) $(src_main_h)
lib_libslater_la_CFLAGS = $(AM_CFLAGS)
lib_libslater_la_CPPFLAGS = "-I${top_srcdir}/src"
lib_libslater_la_CPPFLAGS += "-I${top_srcdir}/include"
lib_libslater_la_LDFLAGS = -version-info 0:0:0
