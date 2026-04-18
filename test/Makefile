PROG = $(patsubst %.cpp,%,$(wildcard *Test.cpp))

REQUIRES =

include $(shell echo $${PREFIX-/usr})/share/smartmet/devel/makefile.inc

FLAGS = -std=$(CXX_STD) -Wall -W -fdiagnostics-color=$(GCC_DIAG_COLOR) -Wno-unused-parameter

CFLAGS = -DUNIX -O0 -g $(FLAGS)

INCLUDES += -I../dynlib

LIBS += ../libsmartmet-dynlib.so \
	$(REQUIRED_LIBS) \
	$(PREFIX_LDFLAGS) \
	-llapack -lgfortran -lm

all: $(PROG)
clean:
	rm -f $(PROG) *~

test: $(PROG)
	@echo Running tests:
	@rc=0; for prog in $(PROG); do \
		LD_LIBRARY_PATH=..:$${LD_LIBRARY_PATH} ./$$prog || rc=1; \
	done; exit $$rc

test-installed: $(PROG)
	@echo Running tests against installed library:
	@rc=0; for prog in $(PROG); do ./$$prog || rc=1; done; exit $$rc

$(PROG) : % : %.cpp ../libsmartmet-dynlib.so
	$(CXX) $(CFLAGS) -o $@ $@.cpp $(INCLUDES) $(LIBS)
