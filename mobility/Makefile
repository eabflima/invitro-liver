#================================================================================
# Directories and Libraries
#================================================================================

INC_PATHS = -I. \
-I$(QUESO_DIR)/include/ \
-I$(LIBMESH_DIR)/include/ \
-I$(BOOST_DIR)/include/ \
-I$(GSL_DIR)/include/ \
-I$(SLEPC_DIR)/include/ \
-I$(PETSC_DIR)/include/

LIBS = -lqueso -L$(QUESO_DIR)/lib/ \
	-lgsl -L$(GSL_DIR)/lib/ \
	-lboost_program_options -L$(BOOST_DIR)/lib

COMP_OPTIONS = -DNDEBUG -g -Wall

# include the library options determined by configure
include $(LIBMESH_DIR)/Make.common

# executable name
target     := ./model_code

#================================================================================
# File management.
#================================================================================

# source files
srcfiles 	:= $(wildcard *.C)

# object files
objects		:= $(patsubst %.C, %.o, $(srcfiles))

###############################################################################

.PHONY: dust clean distclean plot mat figures complete run

###############################################################################

.SUFFIXES: .o .C

# Target:
all:: $(notdir $(target))

# Production rules:  how to make the target - depends on library configuration
$(notdir $(target)): $(objects)
	@echo "Linking "$@"..."
	$(libmesh_CXX) $(libmesh_CXXFLAGS) ${objects} -o $(notdir $(target)) $(libmesh_LIBS) $(LIBS) -Wl,-rpath=$(QUESO_DIR)/lib/

%.o : %.C
	@echo "Compiling C++ (in queso mode) "$<"..."
	$(libmesh_CXX) $(COMP_OPTIONS) $(libmesh_CXXFLAGS) -c $(INC_PATHS) $<

# Useful rules.
dust:
	@echo "Deleting old output and runtime files"
	@rm -rf *~ *# *.e *.txt *.pdf *.jpeg *.png *.eps

clean: dust
	@rm -rf *.o

distclean: clean
	@rm -rf $(notdir $(target)) .libs .depend outputData

run: complete

complete: $(wildcard *.in)
	@$(MAKE) -C $(dir $(target)) $(notdir $(target))
	@echo "***************************************************************"
	@echo "* Running App " $(notdir $(target))
	@echo "***************************************************************"
	@echo " " 
	${LIBMESH_RUN} $(target) queso_options.inp ${LIBMESH_OPTIONS} 2>&1 | tee job_output.txt
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running App " $(notdir $(target))
	@echo "***************************************************************"

plot:
	@./plt.sh

figures:
	@${MATLAB} -nosplash -r "tumor_plots_ip"

mat: figures
	@cp Crop outputData/
	@cd outputData/ ; make -f Crop; pdftk *.pdf cat output allfig.pdf; mv all*.pdf ../

# include the dependency list
-include .depend

#
# Dependencies
#
.depend: $(srcfiles) $(LIBMESH_DIR)/include/libmesh/*.h
	@$(perl) $(LIBMESH_DIR)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(LIBMESH_DIR)/include $(wildcard $(LIBMESH_DIR)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend

###############################################################################
