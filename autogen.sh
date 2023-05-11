#!/bin/sh

aclocal
autoconf
autoheader
automake --add-missing
autoreconf -i -Wall --no-recursive
