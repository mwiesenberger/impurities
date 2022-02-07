device=gpu
FELTOR_PATH=../feltor

#configure machine
include $(FELTOR_PATH)/config/default.mk
include $(FELTOR_PATH)/config/version.mk
include $(FELTOR_PATH)/config/*.mk
include $(FELTOR_PATH)/config/devices/devices.mk

INCLUDE+=-I$(FELTOR_PATH)/inc/


all: impurities impurities_hpc impurities_mpi

impurities: impurities.cpp impurities.h init.h parameters.h diag.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(GLFLAGS) $(LIBS) $(JSONLIB) $(VERSION_FLAGS) -DWITH_GLFW -g

impurities_hpc: impurities.cpp impurities.h init.h parameters.h diag.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) $(VERSION_FLAGS) -g

impurities_mpi: impurities.cpp impurities.h init.h parameters.h diag.h
	$(MPICC) $(OPT) $(MPICFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) $(VERSION_FLAGS)

.PHONY: clean

clean:
	rm -f impurities impurities_hpc impurities_mpi
