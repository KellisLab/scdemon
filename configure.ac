AC_INIT([scdemon], [0.0.3])

dnl AC_MSG_NOTICE([Creating VERSION file])
dnl echo m4_defn([AC_PACKAGE_VERSION]) > VERSION

: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
   AC_MSG_ERROR([Could not determine R_HOME.])
fi

CXX=`"${R_HOME}/bin/R" CMD config CXX`
if test -z "${CXX}"; then
    AC_MSG_ERROR([No C++ compiler available])
fi
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_LANG([C++])
AC_REQUIRE_CPP
AC_PROG_CXX

AC_OPENMP
AC_SUBST([PKG_CXXFLAGS], ["$CXXFLAGS $OPENMP_CXXFLAGS"])
AC_SUBST([PKG_LIBS], ["$LDFLAGS $OPENMP_CXXFLAGS"]) 
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
