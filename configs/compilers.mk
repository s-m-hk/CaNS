FC = mpifort
ifeq ($(strip $(FCOMP)),GNU)
FC = mpifort
endif
ifeq ($(strip $(FCOMP)),INTEL)
FC = mpiifort
CPP = -fpp
endif
ifeq ($(strip $(FCOMP)),NVIDIA)
FC = mpifort
endif
ifeq ($(strip $(FCOMP)),CRAY)
FC = ftn
CPP = -eZ
endif
ifeq ($(strip $(FTN_MPI_WRAPPER)),1)
FC = ftn
endif
