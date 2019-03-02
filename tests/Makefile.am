if BUNDLED_ISL
    ISL_LA = $(top_builddir)/isl/libisl.la
endif
if BUILD_ISL
    ISL_LA = $(ISL_BUILDDIR)/libisl.la
endif

if BUNDLED_BARVINOK
    BARVINOK_LA = $(top_builddir)/barvinok/libbarvinok.la
endif
if BUILD_BARVINOK
    BARVINOK_LA = $(BARVINOK_BUILDDIR)/libbarvinok.la
endif

if BUNDLED_PET
    PET_LA = $(top_builddir)/pet/libpet.la
endif
if BUILD_PET
    PET_LA = $(PET_BUILDDIR)/libpet.la
endif

if DEBUG
    AM_CFLAGS = -g3 -O0
    AM_CXXFLAGS = -g3 -O0
else
    AM_CFLAGS = -O3 -DNDEBUG
    AM_CXXFLAGS = -O3 -DNDEBUG
endif

FORCE:
../src/libhaystack.la: FORCE
	cd ../src; $(MAKE) $(AM_MAKEFLAGS) libhaystack.la

AM_CPPFLAGS = @ISL_CFLAGS@ @BARVINOK_CFLAGS@ @PET_CFLAGS@ 
LIB_ISL = $(ISL_LA) @ISL_LIBS@
LIB_BARVINOK = $(BARVINOK_LA) @BARVINOK_LIBS@
LIB_PET = $(PET_LA) @PET_LIBS@

check_LTLIBRARIES = libgtest.la
libgtest_la_SOURCES = ../googletest/googletest/src/gtest-all.cc
libgtest_la_CPPFLAGS = -I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/googletest/googletest
libgtest_la_LDFLAGS = -pthread

check_PROGRAMS = \
    ToyTest \
    BlinkTest \
    CopyTest \
    CholeskyTest \
    MultiplicationTest \
    StencilTest

BlinkTest_SOURCES = \
    main.cpp \
    CacheEmulator.cpp \
    BlinkTest.cpp

BlinkTest_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) -I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/googletest/googletest -pthread
BlinkTest_LDFLAGS = -pthread
BlinkTest_LDADD = libgtest.la ../src/libhaystack.la $(LIB_ISL) $(LIB_BARVINOK) $(LIB_PET)

CopyTest_SOURCES = \
    main.cpp \
    CacheEmulator.cpp \
    CopyTest.cpp

CopyTest_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) -I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/googletest/googletest -pthread
CopyTest_LDFLAGS = -pthread
CopyTest_LDADD = libgtest.la ../src/libhaystack.la $(LIB_ISL) $(LIB_BARVINOK) $(LIB_PET)


ToyTest_SOURCES = \
    main.cpp \
    CacheEmulator.cpp \
    ToyTest.cpp

ToyTest_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) -I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/googletest/googletest -pthread
ToyTest_LDFLAGS = -pthread
ToyTest_LDADD = libgtest.la ../src/libhaystack.la $(LIB_ISL) $(LIB_BARVINOK) $(LIB_PET)

CholeskyTest_SOURCES = \
    main.cpp \
    CacheEmulator.cpp \
    CholeskyTest.cpp

CholeskyTest_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) -I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/googletest/googletest -pthread
CholeskyTest_LDFLAGS = -pthread
CholeskyTest_LDADD = libgtest.la ../src/libhaystack.la $(LIB_ISL) $(LIB_BARVINOK) $(LIB_PET)

MultiplicationTest_SOURCES = \
    main.cpp \
    CacheEmulator.cpp \
    MultiplicationTest.cpp

MultiplicationTest_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) -I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/googletest/googletest -pthread
MultiplicationTest_LDFLAGS = -pthread
MultiplicationTest_LDADD = libgtest.la ../src/libhaystack.la $(LIB_ISL) $(LIB_BARVINOK) $(LIB_PET)

StencilTest_SOURCES = \
    main.cpp \
    CacheEmulator.cpp \
    ToyTest.cpp

StencilTest_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) -I$(top_srcdir)/googletest/googletest/include -I$(top_srcdir)/googletest/googletest -pthread
StencilTest_LDFLAGS = -pthread
StencilTest_LDADD = libgtest.la ../src/libhaystack.la $(LIB_ISL) $(LIB_BARVINOK) $(LIB_PET)

TESTS = ToyTest 
TESTS += CopyTest 
TESTS += BlinkTest 
TESTS += CholeskyTest 
TESTS += MultiplicationTest 
TESTS += StencilTest