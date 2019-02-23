AC_INIT([haystack], 1.0)
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

AC_PROG_CXX
AC_PROG_LIBTOOL

AX_BOOST_BASE([1.60.0])
AX_BOOST_PROGRAM_OPTIONS

AC_ARG_ENABLE(debug,
AS_HELP_STRING([--enable-debug],
               [enable debugging, default: no]),
[case "${enableval}" in
             yes) debug=true ;;
             no)  debug=false ;;
             *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
esac],
[debug=false])
AM_CONDITIONAL(DEBUG, test x"$debug" = x"true")

AX_SUBMODULE(isl,build|bundled|system,bundled)
AM_CONDITIONAL(BUNDLED_ISL, test $with_isl = bundled)
AM_CONDITIONAL(BUILD_ISL, test $with_isl = build)

AC_SUBST(ISL_CFLAGS)
AC_SUBST(ISL_LIBS)
AC_SUBST(ISL_SRCDIR)
AC_SUBST(ISL_BUILDDIR)
case "$with_isl" in
bundled)
	ISL_CFLAGS="-I\$(top_srcdir)/isl/include -I\$(top_builddir)/isl/include"
	ISL_CFLAGS="$ISL_CFLAGS"
	ISL_SRCDIR="$srcdir/isl"
	ISL_BUILDDIR=isl
	bv_configure_args="$bv_configure_args --with-isl-builddir=../isl"
	PACKAGE_CFLAGS_ISL='-I${prefix}/include'
	;;
build)
	ISL_SRCDIR="$isl_srcdir"
	ISL_BUILDDIR=`echo @abs_builddir@ | $with_isl_builddir/config.status --file=-`
	ISL_CFLAGS="-I$isl_srcdir/include -I$ISL_BUILDDIR/include"
	ISL_LIBS="$with_isl_builddir/libisl.la"
	PACKAGE_CFLAGS_ISL='-I${prefix}/include'
	;;
system)
	PKG_CHECK_MODULES([ISL], [isl])
	PACKAGE_CFLAGS_ISL="$ISL_CFLAGS"
esac
AM_CONDITIONAL(HAVE_ISL_BUILDDIR, test "x$ISL_BUILDDIR" != "x")

AX_SUBMODULE(barvinok,build|bundled|system,bundled)
AM_CONDITIONAL(BUNDLED_BARVINOK, test $with_barvinok = bundled)
AM_CONDITIONAL(BUILD_BARVINOK, test $with_barvinok = build)

AC_SUBST(BARVINOK_CFLAGS)
AC_SUBST(BARVINOK_LIBS)
AC_SUBST(BARVINOK_SRCDIR)
AC_SUBST(BARVINOK_BUILDDIR)
case "$with_barvinok" in
bundled)
	BARVINOK_CFLAGS="-I\$(top_srcdir)/barvinok -I\$(top_builddir)/barvinok"
	BARVINOK_CFLAGS="$BARVINOK_CFLAGS"
	BARVINOK_SRCDIR="$srcdir/barvinok"
	BARVINOK_BUILDDIR=barvinok
	;;
build)
	BARVINOK_SRCDIR="$barvinok_srcdir"
	BARVINOK_CFLAGS="$BARVINOK_CFLAGS"
	BARVINOK_BUILDDIR=`echo @abs_builddir@ | $with_BARVINOK_builddir/config.status --file=-`
	BARVINOK_CFLAGS="-I$barvinok_srcdir/ -I$BARVINOK_BUILDDIR/"
	BARVINOK_LIBS="$with_barvinok_builddir/libisl.la"
	;;
system)
	PKG_CHECK_MODULES([BARVINOK], [barvinok])
	PACKAGE_CFLAGS_BARVINOK="$BARVINOK_CFLAGS"
esac
AM_CONDITIONAL(HAVE_BARVINOK_BUILDDIR, test "x$BARVINOK_BUILDDIR" != "x")

AX_SUBMODULE(pet,bundled|build|system,bundled)
AC_SUBST(PET_CFLAGS)
AC_SUBST(PET_LIBS)
AC_SUBST(PET_BUILDDIR)
case "$with_pet" in
bundled)
	PET_CFLAGS="-I\$(top_srcdir)/pet/include"
	;;
build)
	PET_BUILDDIR=`echo @abs_builddir@ | $with_pet_builddir/config.status --file=-`
	PET_CFLAGS="-I$pet_srcdir/include"
	;;
system)
	PKG_CHECK_MODULES([PET], [pet])
	PACKAGE_CFLAGS_PET="$PET_CFLAGS"
	;;
esac
AM_CONDITIONAL(BUNDLED_PET, test $with_pet = bundled)
AM_CONDITIONAL(BUILD_PET, test $with_pet = build)

PACKAGE_CFLAGS="$PACKAGE_CFLAGS_ISL $PACKAGE_CFLAGS_BARVINOK $PACKAGE_CFLAGS_PET"
PACKAGE_LIBS="-lisl -lbarvinok -lpet"
AX_CREATE_PKGCONFIG_INFO

AC_CONFIG_FILES(Makefile src/Makefile)
if test $with_isl = bundled; then
	AC_CONFIG_SUBDIRS(isl)
fi
if test $with_barvinok = bundled; then
	AC_CONFIG_SUBDIRS(barvinok)
fi
if test $with_pet = bundled; then
	AC_CONFIG_SUBDIRS(pet)
fi
AC_CONFIG_COMMANDS_POST([
	ac_configure_args="$ac_configure_args $bv_configure_args"
])
AC_OUTPUT