
AUTOMAKE_OPTIONS = foreign

include doxygen.am

SUBDIRS = src

control_DATA  = share/licences.txt

controldir = $(prefix)/share/$(PACKAGE)

EXTRA_DIST = README COPYING autogen.sh $(control_DATA)

pkgconfig_DATA = bbcat-control-$(BBCAT_CONTROL_MAJORMINOR).pc

bbcat-control-$(BBCAT_CONTROL_MAJORMINOR).pc: bbcat-control.pc
	cp bbcat-control.pc bbcat-control-$(BBCAT_CONTROL_MAJORMINOR).pc

bbcat-control-$(BBCAT_CONTROL_MAJORMINOR)-uninstalled.pc: bbcat-control-uninstalled.pc
	cp bbcat-control-uninstalled.pc bbcat-control-$(BBCAT_CONTROL_MAJORMINOR)-uninstalled.pc

BUILT_SOURCES = bbcat-control-$(BBCAT_CONTROL_MAJORMINOR)-uninstalled.pc

CLEANFILES = bbcat-control-$(BBCAT_CONTROL_MAJORMINOR).pc bbcat-control-$(BBCAT_CONTROL_MAJORMINOR)-uninstalled.pc src/version.h

ACLOCAL_FLAGS = -I m4
