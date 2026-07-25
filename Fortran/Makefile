
# Author:
#   Yuanqing Wu, KAUST, Saudi Arabia
#
# History:
#   2015-2-9 by Yuanqing Wu
#
# Support:
#   wuyuanq@gmail.com

PROG = MulticomponentFlow
defaults: ${PROG}

FCC=gfortran

#FCFLAGS= -O2 -fdefault-real-8
FCFLAGS= -g -O3
#FCFLAGS= -O2 -fbounds-check 
#FCFLAGS= -O2 

LIBS = -L/opt/local/lib -llapack

FSRCS= RST_model.F90 RST_globalData.F90 RST_PREOS.F90 RST_viscosity.F90 RST_writeFile.F90 \
      RST_genPlotfileMatlab.F90 RST_singlePhaseFlow.F90 RST_singlePhaseFlowDriver.F90 \
      infile_RSTi_singlePhaseFlow.F90 
      
%.o : %.F90
	$(FCC) -c $(FCFLAGS) $^ -o $@

FOBJS = $(FSRCS:.F90=.o)
$(PROG): $(FOBJS)
	$(FCC) $(FCFLAGS) $(FOBJS) $(LIBS) -o $(PROG)

all: $(PROG)

.PHONY: clean

clean:
	rm -rf *.mod *.o $(PROG)

run: ${PROG} 
	./${PROG} 
