if BUNDLED_ISL
    MAYBE_ISL = isl
endif
if BUNDLED_BARVINOK
    MAYBE_BARVINOK = barvinok
endif
if BUNDLED_PET
    MAYBE_PET = pet
endif

SUBDIRS = $(MAYBE_ISL) $(MAYBE_BARVINOK) $(MAYBE_PET) src

ACLOCAL_AMFLAGS = -I m4

FORCE:
isl/libisl.la: FORCE
	cd isl; $(MAKE) $(AM_MAKEFLAGS) libisl.la
barvinok/libbarvinok.la: FORCE
	cd barvinok; $(MAKE) $(AM_MAKEFLAGS) libbarvinok.la
pet/libpet.la: FORCE
	cd pet; $(MAKE) $(AM_MAKEFLAGS) libpet.la

