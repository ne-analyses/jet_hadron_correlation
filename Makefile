os = $(shell uname -s)

#INCFLAGS      = -I$(ROOTSYS)/include -I$(FASTJETDIR)/include -I$(PYTHIA8DIR)/include -I$(STARPICOPATH)
INCFLAGS      = -I$(ROOTSYS)/include -I$(FASTJETDIR)/include -I$(PYTHIA8DIR)/include -I$(PYTHIA8DIR)/include/Pythia8 -I$(STARPICODIR) -I/opt/local/include
INCFLAGS      += -I${BASICAJDIR}/src

ifeq ($(os),Linux)
CXXFLAGS      = 
else
CXXFLAGS      = -O -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
## for debugging:
# CXXFLAGS      = -g -O0 -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
endif

ifeq ($(os),Linux)
LDFLAGS       = -g
LDFLAGSS      = -g --shared 
else
LDFLAGS       = -O -Xlinker -bind_at_load -flat_namespace
LDFLAGSS      = -flat_namespace -undefined suppress
LDFLAGSSS     = -bundle
endif

ifeq ($(os),Linux)
CXX          = g++ 
else
CXX          = clang
endif


ROOTLIBS      = $(shell root-config --libs)

LIBPATH       = $(ROOTLIBS) -L$(FASTJETDIR)/lib -L$(PYTHIA8DIR)/lib -L$(STARPICODIR)
LIBS          = -lfastjet -lfastjettools -lpythia8 -lTStarJetPico

LIBPATH       += -L$(BASICAJDIR)/lib
LIBS	      += -lMyJetlib

# for cleanup
SDIR          = src
ODIR          = src/obj
BDIR          = bin


###############################################################################
################### Remake when these headers are touched #####################
###############################################################################


###############################################################################
# standard rules
$(ODIR)/%.o : $(SDIR)/%.cxx $(INCS)
	@echo 
	@echo COMPILING
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

$(BDIR)/%  : $(ODIR)/%.o 
	@echo 
	@echo LINKING
	$(CXX) $(LDFLAGS) $(LIBPATH) $(LIBS) $^ -o $@

###############################################################################

###############################################################################
############################# Main Targets ####################################
###############################################################################
all : $(BDIR)/test $(BDIR)/auau_correlation $(BDIR)/pp_correlation

$(SDIR)/dict.cxx                : $(SDIR)/ktTrackEff.hh
	cd ${SDIR}; rootcint -f dict.cxx -c -I. ./ktTrackEff.hh

$(ODIR)/dict.o                  : $(SDIR)/dict.cxx
$(ODIR)/ktTrackEff.o            : $(SDIR)/ktTrackEff.cxx $(SDIR)/ktTrackEff.hh
$(ODIR)/corrFunctions.o					: $(SDIR)/corrFunctions.cxx $(SDIR)/corrFunctions.hh

#$(ODIR)/qa_v1.o 		: $(SDIR)/qa_v1.cxx
$(ODIR)/test.o			: $(SDIR)/test.cxx
$(ODIR)/auau_correlation	: $(SDIR)/auau_correlation.cxx
$(ODIR)/pp_correlation		: $(SDIR)/pp_correlation.cxx

#data analysis
#$(BDIR)/qa_v1		: $(ODIR)/qa_v1.o
$(BDIR)/test			: $(ODIR)/test.o $(ODIR)/corrFunctions.o
$(BDIR)/auau_correlation		: $(ODIR)/auau_correlation.o $(ODIR)/corrFunctions.o $(ODIR)/ktTrackEff.o $(ODIR)/dict.o
$(BDIR)/pp_correlation			: $(ODIR)/pp_correlation.o	$(ODIR)/corrFunctions.o $(ODIR)/ktTrackEff.o $(ODIR)/dict.o

###############################################################################
##################################### MISC ####################################
###############################################################################

clean :
	@echo 
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -vf $(BDIR)/*
	rm -vf lib/*


