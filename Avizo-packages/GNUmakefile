
# The configuration and architecture.
MAKE_CFG ?= Debug

uname_s := $(shell sh -c 'uname -s 2>/dev/null || echo unknown')
uname_m := $(shell sh -c 'uname -m 2>/dev/null || echo unknown')

ifeq ($(uname_s),Linux)
    ifeq ($(uname_m),x86_64)
        MAKE_UNAME = LinuxAMD64
    else
        ifeq ($(uname_m),i686)
            MAKE_UNAME = Linux
        else
            $(error Unknown Linux machine type '$(uname_m)'.)
        endif
    endif
else
    ifeq ($(uname_s),Darwin)
        MAKE_UNAME = MacX
    else
        $(error Unknown operating system '$(uname_s)'.)
    endif
endif

ARCH=arch-$(MAKE_UNAME)-$(MAKE_CFG)


# Canceling implicit rules
%.o : %.c
%.o : %.cpp
%.o : %.cc
%.o : %.C
%.o : %.p
%.o : %.r
%.o : %.f
%.o : %.F
%.sym : %.def
%.o : %.mod
%.o : %.s
%.s : %.S
% : %.o
%.c : %.y
%.c : %.l
%.r : %.l
%.nl : %c
%.dvi : %.tex
%.tex : %.web
%.tex : %.w
%.p : %.web
%.c : %.w
%.dvi : %.texinfo
%.dvi : %.texi
%.dvi : %.txinfo
%.info : %.texinfo
%.info : %.texi
%.info : %.txinfo



.PHONY: default
default: fast

.PHONY: all fast
all fast: allAvizopackages
	@echo '    OK all'

.PHONY: pedantic
pedantic:
	@rm -f $(filter-out $(wildcard intermediate/$(ARCH)/*-docrun.lastbuild), $(wildcard intermediate/$(ARCH)/*.lastbuild))
	$(MAKE) all

.PHONY: tutorial
tutorial:
	$(MAKE) obj/$(ARCH)/tutorial.lastbuild

obj/$(ARCH)/tutorial.lastbuild: src/tutorial/GNUmakefile src/tutorial/CCPiAccessibleVolume.h src/tutorial/CCPiComputeThreshold.h src/tutorial/CCPiKMeansFilterITK.h src/tutorial/CCPiLabelQuantification.h src/tutorial/CCPiParticle.h src/tutorial/CCPiParticleTracker.h src/tutorial/CCPiParticleTracking.h src/tutorial/CCPiTrack.h src/tutorial/Quan3D.hpp src/tutorial/QuanWorker.hpp src/tutorial/tutorialAPI.h src/tutorial/share/resources/tutorial.rc src/tutorial/CCPiAccessibleVolume.cpp src/tutorial/CCPiComputeThreshold.cpp src/tutorial/CCPiKMeansFilterITK.cpp src/tutorial/CCPiLabelQuantification.cpp src/tutorial/CCPiParticle.cpp src/tutorial/CCPiParticleTracker.cpp src/tutorial/CCPiParticleTracking.cpp src/tutorial/CCPiTrack.cpp src/tutorial/Quan3D.cpp src/tutorial/QuanWorker.cpp src/tutorial/init.cpp src/tutorial/version.cpp
	$(MAKE) -C src/tutorial
	@mkdir -p obj/$(ARCH)
	@touch obj/$(ARCH)/tutorial.lastbuild
	@echo '    OK tutorial'

.PHONY: docAvizopackages
docAvizopackages:
	$(MAKE) obj/$(ARCH)/docAvizopackages.lastbuild

obj/$(ARCH)/docAvizopackages.lastbuild: $(wildcard  src/tutorial/doc)
	@$(MAKE) obj/$(ARCH)/docAvizopackages-docrun.lastbuild
	@mkdir -p obj/$(ARCH)
	@touch obj/$(ARCH)/docAvizopackages.lastbuild
	@echo '    OK docAvizopackages'

obj/$(ARCH)/docAvizopackages-docrun.lastbuild: $(wildcard  src/tutorial/doc)
	mkdir -p product/share/doc/avizo
	/opt/Avizo/bin/doc2html -A -skin AvizoSkin -local product -d product/share/doc/. -P src/tutorial
	@mkdir -p obj/$(ARCH)
	@touch obj/$(ARCH)/docAvizopackages-docrun.lastbuild

.PHONY: doc
doc: docAvizopackages
	@echo '    OK doc'

.PHONY: allAvizopackages
allAvizopackages:
	$(MAKE) obj/$(ARCH)/allAvizopackages.lastbuild

obj/$(ARCH)/allAvizopackages.lastbuild: obj/$(ARCH)/tutorial.lastbuild obj/$(ARCH)/docAvizopackages.lastbuild
	@mkdir -p obj/$(ARCH)
	@touch obj/$(ARCH)/allAvizopackages.lastbuild
	@echo '    OK allAvizopackages'


.PHONY: clean
clean:
	for dir in src/tutorial; do $(MAKE) -C $$dir clean ; done
	rm -rf product/share/doc

.PHONY: dep
dep: 
	for dir in src/tutorial; do $(MAKE) -C $$dir dep ; done

.PHONY: depclean
depclean:
	for dir in src/tutorial; do $(MAKE) -C $$dir depclean ; done

