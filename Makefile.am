
AUTOMAKE_OPTIONS = foreign

include doxygen.am

SUBDIRS = src

dsp_DATA = share/licences.txt

dspdir = $(prefix)/share/$(PACKAGE)

EXTRA_DIST = README COPYING autogen.sh $(dsp_DATA)

pkgconfig_DATA = bbcat-dsp-$(BBCAT_DSP_MAJORMINOR).pc

bbcat-dsp-$(BBCAT_DSP_MAJORMINOR).pc: bbcat-dsp.pc
	cp bbcat-dsp.pc bbcat-dsp-$(BBCAT_DSP_MAJORMINOR).pc

bbcat-dsp-$(BBCAT_DSP_MAJORMINOR)-uninstalled.pc: bbcat-dsp-uninstalled.pc
	cp bbcat-dsp-uninstalled.pc bbcat-dsp-$(BBCAT_DSP_MAJORMINOR)-uninstalled.pc

BUILT_SOURCES = bbcat-dsp-$(BBCAT_DSP_MAJORMINOR)-uninstalled.pc gitbranch.txt

CLEANFILES = bbcat-dsp-$(BBCAT_DSP_MAJORMINOR).pc bbcat-dsp-$(BBCAT_DSP_MAJORMINOR)-uninstalled.pc src/version.h src/version.cpp gitbranch.txt .branch

ACLOCAL_FLAGS = -I m4

.PHONY: gitbranch.txt

gitbranch.txt:
	@-git rev-parse --abbrev-ref HEAD >.branch && (diff gitbranch.txt .branch >/dev/null || cp .branch gitbranch.txt)
	@if [ ! -e "gitbranch.txt" ] ; then touch "gitbranch.txt" ; fi
	@echo "Branch: `cat gitbranch.txt`"
