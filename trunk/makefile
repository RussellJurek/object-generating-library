# retrieve system architecture
ARCH := $(shell uname -s)
$(info Creating busy function C/C++ library and demo programs for ARCH = $(ARCH))

# specify compilation options
CXX = c++
CC = cc
CFLAGS = -O2 -std=c++11 
LDFLAGS = -lm 
ifndef INSTALL_DIR
  INSTALL_DIR = /usr/local/lib
  $(info INSTALL_DIR not defined. Setting to default value. INSTALL_DIR = $(INSTALL_DIR).)
endif
ifndef INSTALL_INC
  INSTALL_INC = /usr/local/include
  $(info INSTALL_INC not defind. Setting to default value. INSTALL_INC = $(INSTALL_INC).)
endif

# set the library creation method
ifeq ($(ARCH),Darwin)
  LIB_METHOD = Libtool -static -o
else 
  LIB_METHOD = ar -cru
endif
$(info Using LIB_METHOD = $(LIB_METHOD))

# determine if cfitsio is installed --- try /usr/local/lib & /usr/local/include/ if nothing specified
ifeq ($(FITS_DIR), )
  THIS_FITS_DIR = /usr/local
  THIS_FITS_LIB = /lib
  THIS_FITS_INC = /include
  ifneq ($(THIS_FITS_DIR), )
    $(info FITS_DIR = $(THIS_FITS_DIR), checking libraries can be found here . . . )
    ifeq ($(shell ls $(THIS_FITS_DIR)$(THIS_FITS_LIB)/libcfitsio.a), )
      THIS_FITS_DIR = 
      $(info ERROR --- $(THIS_FITS_DIR)$(THIS_FITS_LIB)/libcfitsio.a not found.)
    endif
    ifeq ($(shell ls $(THIS_FITS_DIR)$(THIS_FITS_INC)/fitsio.h), )
      THIS_FITS_DIR = 
      $(info ERROR --- $(THIS_FITS_DIR)$(THIS_FITS_INC)/fitsio.h not found.)
    endif
  else
    THIS_FITS_DIR =
  endif
  ifeq ($(THIS_FITS_DIR), )
    $(info ERROR --- Complete cfitsio installation not found. Not using cfitsio. Change environment variable FITS_DIR to \
    a valid installation directory to use a cfitsio installation. Specify using, make FITS_DIR=path-to-installation , \
    FITS_INC=path-to-fitsio.h from FITS_DIR, FITS_LIB=path-to-libcfitsio.a from FITS_DIR .)
  else
    $(info cfitsio installation found at $(THIS_FITS_DIR). Using cfitsio.)
  endif
else
  THIS_FITS_DIR = $(FITS_DIR)
  THIS_FITS_LIB = $(FITS_LIB)
  THIS_FITS_INC = $(FITS_INC)
  ifneq ($(THIS_FITS_DIR), )
    $(info FITS_DIR = $(THIS_FITS_DIR), checking libraries can be found here . . . )
    ifeq ($(shell ls $(THIS_FITS_DIR)/$(THIS_FITS_LIB)/libwcs.a), )
      THIS_FITS_DIR = 
    $(info ERROR --- $(THIS_FITS_DIR)/$(THIS_FITS_LIB)/libwcs.a not found.)
    endif
    ifeq ($(shell ls $(THIS_FITS_DIR)/$(THIS_FITS_INC)/wcs.h), )
      THIS_FITS_DIR = 
    $(info ERROR --- $(THIS_FITS_DIR)/$(THIS_FITS_INC)/wcs.h not found.)
    endif
  else
    THIS_FITS_DIR =
  endif
  ifeq ($(THIS_FITS_DIR), )
    $(info ERROR --- Complete wcslib installation not found. Not using wcslib. Change environment variable FITS_DIR to \
    a valid installation directory to use a cfitsio installation. Specify using, make FITS_DIR=path-to-installation , \
    FITS_INC=path-to-fitsio.h from FITS_DIR, FITS_LIB=path-to-libcfitsio.a from FITS_DIR .)
  else
    $(info wcslib installation found at $(THIS_FITS_DIR)/$(THIS_FITS_LIB) & $(THIS_FITS_DIR)/$(THIS_FITS_INC). Using wcslib.)
  endif
endif

# determine if pgplot is installed
THIS_PGPLOT_DIR = $(PGPLOT_DIR)
ifneq ($(THIS_PGPLOT_DIR), )
  $(info PGPLOT_DIR = $(THIS_PGPLOT_DIR), checking libraries can be found here . . . )
  ifeq ($(shell ls $(THIS_PGPLOT_DIR)/libpgplot.a), )
    THIS_PGPLOT_DIR = 
  $(info ERROR --- libpgplot.a not found.)
  endif
  ifeq ($(shell ls $(THIS_PGPLOT_DIR)/libcpgplot.a), )
    THIS_PGPLOT_DIR = 
  $(info ERROR --- libcpgplot.a not found.)
  endif
  ifeq ($(shell ls $(THIS_PGPLOT_DIR)/cpgplot.h), )
    THIS_PGPLOT_DIR = 
  $(info ERROR --- cpgplot.h not found.)
  endif
else
  THIS_PGPLOT_DIR =
endif
ifeq ($(THIS_PGPLOT_DIR), )
  $(info ERROR --- Complete pgplot installation not found. Not using pgplot. Change environment variable PGPLOT_DIR to \
  a valid installation directory to use a pgplot installation. Specify using, make PGPLOT_DIR=path-to-installation .)
else
  $(info Pgplot installation found at $(THIS_PGPLOT_DIR). Using pgplot.)
endif

# determine if wcslib is installed --- try /usr/local/lib & /usr/local/include/wcslib if nothing specified
# use libwcs.a as the default library
ifeq ($(WCS_LIB), )
  WCS_LIB = -lwcs
endif
ifeq ($(WCS_DIR), )
  THIS_WCS_DIR = /usr/local
  THIS_WCS_INC = /include/wcslib/
  THIS_WCS_LIB = /lib
  ifneq ($(THIS_WCS_DIR), )
    $(info WCS_DIR = $(THIS_WCS_DIR), checking libraries can be found here . . . )
    ifeq ($(shell ls $(THIS_WCS_DIR)$(THIS_WCS_LIB)/libwcs.a), )
      THIS_WCS_DIR = 
    $(info ERROR --- $(THIS_WCS_DIR)$(THIS_WCS_LIB)/libwcs.a not found.)
    endif
    ifeq ($(shell ls $(THIS_WCS_DIR)$(THIS_WCS_INC)/wcs.h), )
      THIS_WCS_DIR = 
    $(info ERROR --- $(THIS_WCS_DIR)$(THIS_WCS_INC)/wcs.h not found.)
    endif
  else
    THIS_WCS_DIR =
  endif
  ifeq ($(THIS_WCS_DIR), )
    $(info ERROR --- Complete wcslib installation not found. Not using wcslib. Change environment variable WCS_DIR to \
    a valid installation directory to use a cfitsio installation. Specify using, make WCS_DIR=path-to-installation , \
    WCS_INC=path-to-wcs.h from WCS_DIR, WCS_LIB=path-to-libwcs.a from WCS_DIR .)
  else
    $(info wcslib installation found at $(THIS_WCS_DIR). Using wcslib.)
  endif
else
  THIS_WCS_DIR = $(WCS_DIR)
  THIS_WCS_LIB = $(WCS_LIB)
  THIS_WCS_INC = $(WCS_INC)
  ifneq ($(THIS_WCS_DIR), )
    $(info WCS_DIR = $(THIS_WCS_DIR), checking libraries can be found here . . . )
    ifeq ($(shell ls $(THIS_WCS_DIR)/$(THIS_WCS_LIB)/libwcs.a), )
      THIS_WCS_DIR = 
    $(info ERROR --- $(THIS_WCS_DIR)/$(THIS_WCS_LIB)/libwcs.a not found.)
    endif
    ifeq ($(shell ls $(THIS_WCS_DIR)/$(THIS_WCS_INC)/wcs.h), )
      THIS_WCS_DIR = 
    $(info ERROR --- $(THIS_WCS_DIR)/$(THIS_WCS_INC)/wcs.h not found.)
    endif
  else
    THIS_WCS_DIR =
  endif
  ifeq ($(THIS_WCS_DIR), )
    $(info ERROR --- Complete wcslib installation not found. Not using wcslib. Change environment variable WCS_DIR to \
    a valid installation directory to use a cfitsio installation. Specify using, make WCS_DIR=path-to-installation , \
    WCS_INC=path-to-wcs.h from WCS_DIR, WCS_LIB=path-to-libwcs.a from WCS_DIR .)
  else
    $(info wcslib installation found at $(THIS_WCS_DIR)/$(THIS_WCS_LIB) & $(THIS_WCS_DIR)/$(THIS_WCS_INC). Using wcslib.)
  endif
endif

# define rules used to create BusyFunction library and test programs
all: ./librjj_objgen.a ./librjj_objgen_plots.a ./librjj_objgen_wcs.a ./create_catalog_NP ./create_catalog ./ProcessSourceOnlyCube

ifneq ($(THIS_FITS_DIR), )

SOURCES = RJJ_ObjGen_DetectDefn.cpp RJJ_ObjGen_CatPrint.cpp RJJ_ObjGen_PlotGlobal.cpp RJJ_ObjGen_CreateObjs.cpp RJJ_ObjGen_ThreshObjs.cpp RJJ_ObjGen_AddObjs.cpp RJJ_ObjGen_MemManage.cpp RJJ_ObjGen_MakeMask.cpp RJJ_ObjGen_Dmetric.cpp
OBJECTS = $(SOURCES:.cpp=.o)
./librjj_objgen.a: $(SOURCES)
	@echo Creating C++ library --- INCLUDING cfitsio extensions . . . 
	$(CXX) $(CFLAGS) -c $(SOURCES) -I. -I/$(THIS_FITS_DIR)/$(THIS_FITS_INC)/ -L/$(THIS_FITS_DIR)/$(THIS_FITS_LIB)/
	$(LIB_METHOD) $@ $(OBJECTS) -L/$(THIS_FITS_DIR)/$(THIS_FITS_LIB)/ -lcfitsio

./create_catalog_NP: ./create_catalog_NP.cpp ./librjj_objgen.a
	@echo Creating terminal application create_catalog_NP . . . 
	$(CXX) $< -o $@ $(CFLAGS) -I. -I/$(THIS_FITS_DIR)/$(THIS_FITS_INC)/ -L. -lrjj_objgen -L/$(THIS_FITS_DIR)/$(THIS_FITS_LIB)/ -lcfitsio 

else

SOURCES = RJJ_ObjGen_DetectDefn.cpp RJJ_ObjGen_CatPrint.cpp RJJ_ObjGen_PlotGlobal.cpp RJJ_ObjGen_CreateObjs.cpp RJJ_ObjGen_ThreshObjs.cpp RJJ_ObjGen_AddObjs.cpp RJJ_ObjGen_MemManage.cpp RJJ_ObjGen_Dmetric.cpp
OBJECTS = $(SOURCES:.c=.o)
./librjj_objgen.a: $(SOURCES)
	@echo Creating C++ library --- EXCLUDING cfitsio extensions . . . 
	$(CXX) $(CFLAGS) -c $(SOURCES) -I. 
	$(LIB_METHOD) $@ $(OBJECTS)

./create_catalog_NP: ./create_catalog_NP.cpp ./librjj_objgen.a
	@echo NOT creating terminal application create_catalog_NP . . . Couldn\'t find cfitsio installation.

endif

ifneq ($(THIS_PGPLOT_DIR), )

./librjj_objgen_plots.a: RJJ_ObjGen_Plots.cpp 
	@echo Creating C++ plotting library . . .
	$(CXX) -c $< $(CFLAGS) -I. -I/$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_INC)/ -L/$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_LIB)/
	$(LIB_METHOD) $@ RJJ_ObjGen_Plots.o -L. -lrjj_objgen -L/$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_LIB)/ -lcpgplot -lpgplot

ifneq ($(THIS_FITS_DIR), )

./create_catalog: ./create_catalog.cpp ./librjj_objgen.a ./librjj_objgen_plots.a
	@echo Creating terminal application create_catalog . . . 
	$(CXX) $< -o $@ $(CFLAGS) -I. -I/$(THIS_FITS_DIR)/$(THIS_FITS_INC)/ -L/$(THIS_FITS_DIR)/$(THIS_FITS_LIB)/ -I/$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_INC)/ -L/$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_LIB)/ -lcpgplot -lpgplot -L. -lrjj_objgen -lrjj_objgen_plots -lcfitsio 

./ProcessSourceOnlyCube: ./ProcessSourceOnlyCube.cpp ./librjj_objgen.a ./librjj_objgen_plots.a
	@echo Creating terminal application ProcessSourceOnlyCube . . . 
	$(CXX) $< -o $@ $(CFLAGS) -I. -I/$(THIS_FITS_DIR)/$(THIS_FITS_INC)/ -I/$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_INC)/ -L/$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_LIB)/ -lcpgplot -lpgplot -L. -lrjj_objgen -lrjj_objgen_plots -L/$(THIS_FITS_DIR)/$(THIS_FITS_LIB)/ -lcfitsio 

else

./create_catalog: ./create_catalog.cpp ./librjj_objen.a
	@echo NOT creating terminal application create_catalog . . . Couldn\'t find cfitsio and pgplot installation.

./ProcessSourceOnlyCube: ./ProcessSourceOnlyCube.cpp ./librjj_objgen.a ./librjj_objgen_plots.a
	@echo NOT creating terminal application ProcessSourceOnlyCube . . . Couldn\'t find cfitsio and pgplot installation.

endif

else

./librjj_objgen_plots.a: RJJ_ObjGen_Plots.cpp
	@echo NOT creating C++ plotting library . . . Couldn\'t find pgplot installation.

endif

ifneq ($(THIS_WCS_DIR), )

./librjj_objgen_wcs.a: RJJ_ObjGen_CatPrint_WCS.cpp
	@echo Creating C++ WCS library . . .
	$(CXX) -c $< $(CFLAGS) -I. -I/$(THIS_WCS_DIR)/$(THIS_WCS_INC)/ -L/$(THIS_WCS_DIR)/$(THIS_WCS_LIB)/
	$(LIB_METHOD) $@ RJJ_ObjGen_CatPrint_WCS.o -L/. -L/$(THIS_WCS_DIR)/$(THIS_WCS_LIB)/ $(WCS_LIB) 

else

./librjj_objgen_wcs.a: RJJ_ObjGen_CatPrint_WCS.cpp
	@echo NOT creating C++ WCS library . . . Couldn\'t find WCS installation.

endif

clean:
	@echo Removing all object files . . . 
	@rm -vf *.o

distclean:
	@echo Removing all object files, libraries and compiled programs . . . 
	@rm -vf *.o
	@rm -vf librjj_objgen.a
	@rm -vf librjj_objgen_plots.a
	@rm -vf librjj_objgen_wcs.a
	@rm -vf create_catalog
	@rm -vf create_catalog_NP
	@rm -vf ProcessSourceOnly

install:
	@echo Copying libraries and include files to $(INSTALL_DIR) and $(INSTALL_INC) . . .
	@cp librjj_objgen.a $(INSTALL_DIR)
	@cp librjj_objgen_plots.a $(INSTALL_DIR)
	@cp librjj_objgen_wcs.a $(INSTALL_DIR)
	@cp RJJ_ObjGen.h $(INSTALL_INC)
	@cp RJJ_ObjGen_Plots.h $(INSTALL_INC)
	@cp RJJ_ObjGen_WCS.h $(INSTALL_INC)

