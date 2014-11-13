
AUTOMAKE_OPTIONS = foreign

include doxygen.am

SUBDIRS = src test

dsp_DATA  = share/licences.txt

dspdir = $(prefix)/share/$(PACKAGE)

EXTRA_DIST = README COPYING autogen.sh $(dsp_DATA)

pkgconfig_DATA = bbcat-dsp-$(BBCAT_DSP_MAJORMINOR).pc

bbcat-dsp-$(BBCAT_DSP_MAJORMINOR).pc: bbcat-dsp.pc
	cp bbcat-dsp.pc bbcat-dsp-$(BBCAT_DSP_MAJORMINOR).pc

bbcat-dsp-$(BBCAT_DSP_MAJORMINOR)-uninstalled.pc: bbcat-dsp-uninstalled.pc
	cp bbcat-dsp-uninstalled.pc bbcat-dsp-$(BBCAT_DSP_MAJORMINOR)-uninstalled.pc

BUILT_SOURCES = bbcat-dsp-$(BBCAT_DSP_MAJORMINOR)-uninstalled.pc

CLEANFILES = bbcat-dsp-$(BBCAT_DSP_MAJORMINOR).pc bbcat-dsp-$(BBCAT_DSP_MAJORMINOR)-uninstalled.pc

ACLOCAL_FLAGS = -I m4
