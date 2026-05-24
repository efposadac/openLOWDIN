TOPDIR=..

include $(TOPDIR)/CONFIG

# Build folder to keep .o .mod and .a
BUILD_DIR = $(abspath $(TOPDIR)/build)

# Let Fortran knows the building dir 
FCFLAGS += -J $(BUILD_DIR)


# openlowdin library for each folder here with a configure inside
SUBDIRS = Core Output CI MM MBPT PT COSMO Ints Deriv IntegralsTransformation CalcProp DFT SCF NOCI Solver
#SUBDIRS := $(patsubst %/,%,$(dir $(wildcard */configure)))

# Reverse function
reverse = $(if $(1),$(call reverse,$(wordlist 2,$(words $(1)),$(1))) $(firstword $(1)))

# Reverse the SUBDIRS list so the linker gets the correct dependency order
LINK_ORDER = $(call reverse,$(SUBDIRS))

# Automatically generate library paths in correct order
STATIC_LIBS = $(foreach dir,$(LINK_ORDER),$(BUILD_DIR)/$(dir)/$(dir).a)

# External libraries
LIBS=$(ALL_LIBS)
LIBS += $(TOPDIR)/lib/libtransform.a 
LIBS += $(TOPDIR)/lib/libjadamilu.a
LIBS += $(TOPDIR)/lib/libmolden2aim.a

TARGET=$(BUILD_DIR)/lowdin.x

.PHONY: all subdirs $(SUBDIRS) clean setup

all: setup subdirs $(TARGET)

# Creating build folders
setup:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(foreach dir, $(SUBDIRS), $(BUILD_DIR)/$(dir))

# declare the internal dependencies here for each library
subdirs: $(SUBDIRS)

include dependencies

# Make for each library
$(SUBDIRS):
	$(MAKE) -C $@

# Linking Step
$(TARGET): $(BUILD_DIR)/lowdin.x.o
	$(FC) $(FCFLAGS) -o $@ $< $(STATIC_LIBS) $(LIBS) 
	mv $(TARGET) $(TOPDIR)/bin/

# Compiling main 
$(BUILD_DIR)/lowdin.x.o: lowdin.x.f90 subdirs 
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) 
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir clean; done

