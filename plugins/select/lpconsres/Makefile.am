# Makefile for select/lpconsres plugin

AUTOMAKE_OPTIONS = foreign

PLUGIN_FLAGS = -module -avoid-version --export-dynamic

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src/common

pkglib_LTLIBRARIES = select_lpconsres.la

# Consumable resources node selection plugin.
select_lpconsres_la_SOURCES =  select_lpconsres.c select_lpconsres.h \
                              dist_tasks.c dist_tasks.h \
			      job_test.c job_test.h
select_lpconsres_la_LDFLAGS = $(SO_LDFLAGS) $(PLUGIN_FLAGS)
