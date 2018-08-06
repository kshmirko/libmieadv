#!/usr/bin/make

#main building variables
DSRC    = src
DOBJ    = build/obj/
DMOD    = build/mod/
DEXE    = build/
LIBS    =
FC      = gfortran
OPTSC   = -c -fPIC -J build/mod
OPTSL   = -shared -O3 -J build/mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)MIEV0MOD: $(MKDIRS) $(DOBJ)miev0mod.o \
	$(DOBJ)miev0.o \
	$(DOBJ)errpack.o
	@rm -f $(filter-out $(DOBJ)miev0mod.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) MIEV0MOD

#compiling rules
$(DOBJ)miev0.o: src/wiscombe/MIEV0.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)errpack.o: src/wiscombe/ErrPack.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)rdi1mach.o: src/wiscombe/RDI1MACH.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)miev0mod.o: src/wiscombe/miev0mod.f03 \
	$(DOBJ)mathutils.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mathutils.o: src/mathutils/mathutils.f03
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
