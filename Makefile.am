
AUTOMAKE_OPTIONS = foreign

include doxygen.am

SUBDIRS = src

EXTRA_DIST = README COPYING autogen.sh

pkgconfig_DATA = bbcat-dsp-$(BBCAT_DSP_MAJORMINOR).pc

bbcat-dsp-$(BBCAT_DSP_MAJORMINOR).pc: bbcat-dsp.pc
	cp bbcat-dsp.pc bbcat-dsp-$(BBCAT_DSP_MAJORMINOR).pc

bbcat-dsp-$(BBCAT_DSP_MAJORMINOR)-uninstalled.pc: bbcat-dsp-uninstalled.pc
	cp bbcat-dsp-uninstalled.pc bbcat-dsp-$(BBCAT_DSP_MAJORMINOR)-uninstalled.pc

BUILT_SOURCES = bbcat-dsp-$(BBCAT_DSP_MAJORMINOR)-uninstalled.pc

CLEANFILES = bbcat-dsp-$(BBCAT_DSP_MAJORMINOR).pc bbcat-dsp-$(BBCAT_DSP_MAJORMINOR)-uninstalled.pc

ACLOCAL_FLAGS = -I m4
