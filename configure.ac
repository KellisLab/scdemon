AC_INIT([scdemon], [0.0.2])

AC_MSG_NOTICE([Creating VERSION file])
echo m4_defn([AC_PACKAGE_VERSION]) > VERSION

AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET

AC_LANG([C++])
AC_PROG_CXX

PKG_PROG_PKG_CONFIG
PKG_CHECK_MODULES([EIGEN], [eigen3 >= 3.3.9], [], [AC_MSG_WARN([Eigen flags may not be set for python])])
PKG_CHECK_MODULES([HDF5], [hdf5], [], [AC_MSG_WARN([HDF5 may not be set for testing])])
PKG_CHECK_MODULES([GSL], [gsl], [], [AC_MSG_ERROR([GSL is required, please install])])
AC_OPENMP

AC_CHECK_PROG([R_FOUND], [R], [yes], [no])
if test "x$R_FOUND" = xno; then
  AC_MSG_WARN([R is not installed.])
fi
AC_PROG_MAKE_SET

AM_INIT_AUTOMAKE([foreign subdir-objects])
AM_MAINTAINER_MODE


AC_ARG_ENABLE([cxxtest],
    [AS_HELP_STRING([--enable-cxxtest], [Enable the compilation of c++ binary])],
    [cxxtest=$enableval],
    [cxxtest=no])

AM_CONDITIONAL([ENABLE_CXX_TEST], [test "x$cxxtest" = "xyes"])

AC_SUBST([PKG_CXXFLAGS], ["$CXXFLAGS $OPENMP_CXXFLAGS $GSL_CFLAGS"])
AC_SUBST([PKG_LIBS], ["$LDFLAGS $OPENMP_CXXFLAGS $GSL_LIBS"]) 
AC_SUBST([EIGEN_CXXFLAGS], ["$EIGEN_CFLAGS"])
AC_SUBST([EIGEN_LIBS], ["$EIGEN_LIBS"])
AC_SUBST([OPENMP_CXXFLAGS], ["$OPENMP_CXXFLAGS"])
AC_SUBST([OPENMP_LIBS], ["$OPENMP_CXXFLAGS"])
AC_SUBST([HDF5_CXXFLAGS], ["$HDF5_CXXFLAGS"])
AC_SUBST([HDF5_LIBS], ["$HDF5_LIBS"])
AC_SUBST([GSL_CXXFLAGS], ["$GSL_CFLAGS"])
AC_SUBST([GSL_LIBS], ["$GSL_LIBS"])
AC_CONFIG_FILES([src/Makevars Makefile DESCRIPTION docs/sphinx/conf.py])
AC_OUTPUT
