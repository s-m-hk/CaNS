override LIBS += -L$(LIBS_DIR)/2decomp-fft -ldecomp2d
override INCS += -I$(LIBS_DIR)/2decomp-fft/mod

ifeq ($(strip $(GPU)),1)
override LIBS += -L$(LIBS_DIR)/cuDecomp/build/lib -lcudecomp -lcudecomp_fort -cudalib=cufft
override INCS += -I$(LIBS_DIR)/cuDecomp/build/include
endif

ifeq ($(strip $(USE_NVTX)),1)
override LIBS += -L$(NVHPC_ROOT)/cuda/lib64 -lnvToolsExt
endif

#ifeq ($(strip $(FCOMP)),GNU)
#override LIBS += -L/lscratch/smhk2/Software/AMD/amd-fftw/lib
#endif
#ifeq ($(strip $(FCOMP)),INTEL)
#override LIBS += -L/scratch/smhk2/Software/fftw-intel/lib
#endif
#ifeq ($(strip $(FCOMP)),NVIDIA)
#override LIBS += -L$(NVHPC_ROOT)/math_libs/lib64
#endif

ifneq ($(strip $(GPU)),1)
override LIBS += -lfftw3

ifeq ($(strip $(OPENMP)),1)
override LIBS += -lfftw3_threads
endif

ifeq ($(strip $(SINGLE_PRECISION)),1)
override LIBS += -lfftw3f
ifeq ($(strip $(OPENMP)),1)
override LIBS += -lfftw3f_threads
endif
endif

ifeq ($(strip $(SINGLE_PRECISION_POISSON)),1)
override LIBS += -lfftw3f
ifeq ($(strip $(OPENMP)),1)
override LIBS += -lfftw3f_threads
endif
endif

endif