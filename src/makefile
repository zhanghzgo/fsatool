CONFIG_FILE=config.h
PATCH_DIR=fsatool1_pmemd_src/amber_codes
AMBER_SRC=$(BASEDIR)/src/pmemd/src
SHELL := /bin/bash

include $(CONFIG_FILE)

.PHONY: install clean

install:
	make -C $(FSATOOL_SAMPLE_SRC)
	make -C $(AMBERHOME)/src/pmemd
	make -C fsatool2_extract/
	make -C fsatool3_msm/
	sed -i 's|\$$FSATOOLHOME|$(FSATOOLHOME)|g' fsatool

clean:
	(cd $(FSATOOL_SAMPLE_SRC) && make clean)
	(cd fsatool2_extract/ && make clean)
	(cd fsatool3_msm/ && make clean)
	sed -i 's|$(FSATOOLHOME)|\$$FSATOOLHOME|g' fsatool
