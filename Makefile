SUBNAME = dynlib
LIB = smartmet-$(SUBNAME)
SPEC = smartmet-library-$(SUBNAME)
INCDIR = smartmet/$(SUBNAME)

DEFINES = -DUNIX -D_REENTRANT

# macgyver is used only for the header-only Fmi::Matrix template.
REQUIRES :=

include $(shell smartbuildcfg --prefix)/share/smartmet/devel/makefile.inc

# Fortran sources are vendored under third_party/dynlib/. The compile
# order below respects the module USE graph: kind -> const -> config ->
# derivatives -> diag / detect_rwb_contour / detect_lines / utils -> detect ->
# dynlib_wrapper (our ISO_C_BINDING shim).
FORTRAN_DIR  = third_party/dynlib
FORTRAN_SRCS = \
	kind.f90 \
	const.f90 \
	config.f90 \
	derivatives.f90 \
	diag.f90 \
	detect_rwb_contour.f90 \
	detect_lines.f90 \
	utils.f90 \
	detect.f90 \
	dynlib_wrapper.f90
FORTRAN_OBJS = $(patsubst %.f90, $(objdir)/%.o, $(FORTRAN_SRCS))

FC     := gfortran
# Upstream dynlib uses `intent(out)` on output buffers that are only
# partially written (e.g. `line_locate` fills the valid prefix of lnoff
# and leaves the rest). Under gfortran's aggressive -O1/-O2 tree-level
# optimisations the shim's pre-call NaN fill of those buffers is seen
# as a dead store and elided; unused slots then contain stack leftovers
# from neighbouring frames, and the C decoder's monotonic-offset scan
# reads them as spurious lines. This is strictly deterministic but
# depends on caller stack layout, so extra test binary symbols can
# shift results from "78 points detected" to "0 points detected". Use
# -O0 for the vendored Fortran until we port the detection code or
# refactor the shim to not depend on `intent(out)` pre-fill semantics.
FFLAGS := -O0 -g -fPIC -fno-range-check -J$(objdir) -cpp

# C++ sources
SRC_DIRS = $(SUBNAME)
vpath %.cpp $(SRC_DIRS)
vpath %.h   $(SRC_DIRS)

SRCS = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.cpp))
HDRS = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.h))
OBJS = $(patsubst %.cpp, $(objdir)/%.o, $(notdir $(SRCS)))

INCLUDES := -Idynlib $(INCLUDES)

LIBFILE = libsmartmet-$(SUBNAME).so

# Runtime deps: LAPACK for diag.f90, libgfortran for the Fortran runtime.
LIBS += \
	$(PREFIX_LDFLAGS) \
	$(REQUIRED_LIBS) \
	-llapack \
	-lgfortran \
	-lm

.PHONY: test rpm all clean format examples install debug release

all: objdir $(LIBFILE)
debug: all
release: all

# End-to-end smoke test against a real querydata file. Links against
# the just-built libsmartmet-dynlib.so plus the installed
# smartmet-library-newbase headers and .so.
examples: $(LIBFILE) examples/demo

examples/demo: examples/demo.cpp $(LIBFILE)
	$(CXX) $(CFLAGS) -Idynlib -isystem /usr/include/smartmet \
		examples/demo.cpp ./$(LIBFILE) \
		-lsmartmet-newbase -lsmartmet-macgyver -lboost_iostreams \
		-llapack -lgfortran -lm \
		-Wl,-rpath,$(CURDIR) \
		-o examples/demo

$(LIBFILE): $(FORTRAN_OBJS) $(OBJS)
	$(CXX) $(LDFLAGS) -shared -rdynamic -o $(LIBFILE) $(FORTRAN_OBJS) $(OBJS) $(LIBS)
	@echo Checking $(LIBFILE) for unresolved references
	@if ldd -r $(LIBFILE) 2>&1 | c++filt | grep ^undefined\ symbol | \
			grep -Pv ':\ __(?:(?:a|t|ub)san_|sanitizer_)'; \
		then rm -v $(LIBFILE); \
		exit 1; \
	fi

clean:
	rm -f $(LIBFILE) *~ $(SUBNAME)/*~ $(FORTRAN_DIR)/*~
	rm -rf $(objdir)
	$(MAKE) -C test clean

format:
	clang-format -i -style=file $(SUBNAME)/*.h $(SUBNAME)/*.cpp test/*.cpp

install:
	@mkdir -p $(includedir)/$(INCDIR)
	@list='$(HDRS)'; for hdr in $$list; do \
		HDR=$$(echo $$hdr | sed -e 's:^$(SUBNAME)/::'); \
		$(INSTALL_DATA) $$hdr $(includedir)/$(INCDIR)/$$HDR; \
	done
	@mkdir -p $(libdir)
	$(INSTALL_PROG) $(LIBFILE) $(libdir)/$(LIBFILE)

test test-installed:
	$(MAKE) -C test $@

rpm: clean $(SPEC).spec
	rm -f $(SPEC).tar.gz
	tar -czvf $(SPEC).tar.gz --exclude-vcs --transform "s,^,$(SPEC)/," *
	rpmbuild -tb $(SPEC).tar.gz $(RPMBUILD_OPT)
	rm -f $(SPEC).tar.gz

# C++ compilation
$(objdir)/%.o : %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) $(INCLUDES) -c -MD -MF $(patsubst $(objdir)/%.o, $(objdir)/%.d, $@) -MT $@ -o $@ $<

# Fortran compilation. Each step produces a .mod into $(objdir) (-J flag);
# compilation order is enforced via explicit dependencies.
$(objdir)/%.o : $(FORTRAN_DIR)/%.f90
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c $< -o $@

# Fortran module dependency chain (required because gfortran reads .mod
# files, so each dependent must be compiled after its dependencies).
$(objdir)/config.o:             $(objdir)/kind.o
$(objdir)/const.o:              $(objdir)/kind.o
$(objdir)/derivatives.o:        $(objdir)/kind.o $(objdir)/config.o $(objdir)/const.o
$(objdir)/diag.o:               $(objdir)/derivatives.o
$(objdir)/detect_rwb_contour.o: $(objdir)/kind.o $(objdir)/config.o
$(objdir)/detect_lines.o:       $(objdir)/derivatives.o $(objdir)/const.o
$(objdir)/utils.o:              $(objdir)/kind.o $(objdir)/config.o $(objdir)/const.o
$(objdir)/detect.o:             $(objdir)/diag.o $(objdir)/detect_rwb_contour.o \
                                $(objdir)/detect_lines.o $(objdir)/utils.o
$(objdir)/dynlib_wrapper.o:     $(objdir)/detect.o

objdir:
	@mkdir -p $(objdir)

ifneq ($(wildcard $(objdir)/*.d),)
-include $(wildcard $(objdir)/*.d)
endif
