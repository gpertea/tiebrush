## add additional library path here if needed (where libhts.a was installed)
XLIB := /ccb/sw/lib
GDIR := ../gclib
PKGPATH := $(if $(XLIB),export PKG_CONFIG_PATH=$(XLIB)/pkgconfig,true)
HTSLIBS := $(shell ${PKGPATH}; pkg-config --libs --static htslib)
ifeq ($(HTSLIBS),)
 $(error ERROR: htslib installation not found. Please install htslib first)
endif

HTSINC := $(shell ${PKGPATH}; pkg-config --cflags htslib)
HLDIR := $(shell ${PKGPATH};pkg-config --libs-only-L htslib)
HLDIR := $(HLDIR:-L%=%)


#ifeq (${HLDIR}/libhts.a,$(wildcard ${HLDIR}/libhts.a))
# HTSLIBS := $(subst -lhts,${HLDIR}/libhts.a,${HTSLIBS})
#endif


### libdeflate speeds up de/compression of BAM/CRAM
# leave it blank if not available
#LIBDEFLATE :=
#LIBDEFLATE := $(if $(LIBDEFLATE),$(LIBDEFLATE),/ccb/sw/lib/libdeflate.a)

INCDIRS := -I${GDIR} ${HTSINC}

CXX   := $(if $(CXX),$(CXX),g++)

BASEFLAGS := -Wall -Wextra ${INCDIRS} -fsigned-char -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -std=c++0x -fno-strict-aliasing -fno-exceptions -fno-rtti
#for gcc 8+ add: -Wno-class-memaccess
GCCVER5 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 5)
ifeq "$(GCCVER5)" "1"
 BASEFLAGS += -Wno-implicit-fallthrough
endif

GCCVER8 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 8)
ifeq "$(GCCVER8)" "1"
  BASEFLAGS += -Wno-class-memaccess
endif

LINKER  := $(if $(LINKER),$(LINKER),g++)

LDFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)


#LIBS :=${SAMLIB}/libhts.a ${LIBDEFLATE} -llzma -lbz2 -lz -lm -lcurl -lcrypto

LIBS := ${HTSLIBS}

#ifneq (,$(findstring nothreads,$(MAKECMDGOALS)))
# NOTHREADS=1
#endif

#detect MinGW (Windows environment)
ifneq (,$(findstring mingw,$(shell ${CXX} -dumpmachine)))
 WINDOWS=1
endif

# Misc. system commands
ifdef WINDOWS
 RM = del /Q
else
 RM = rm -f
endif

# File endings
ifdef WINDOWS
 EXE = .exe
else
 EXE =
endif

DMACH := $(shell ${CXX} -dumpmachine)

ifneq (,$(filter %release %static, $(MAKECMDGOALS)))
  # -- release build
  RELEASE_BUILD=1
  CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O2)
  CXXFLAGS += -DNDEBUG $(BASEFLAGS)
else
  ifneq (,$(filter %memcheck %memdebug %tsan %tcheck %thrcheck, $(MAKECMDGOALS)))
     #use sanitizer in gcc 4.9+
     GCCVER49 := $(shell expr `${CXX} -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
     SANLIBS :=
     ifneq (,$(filter %tsan %tcheck %thrcheck, $(MAKECMDGOALS)))
        # thread sanitizer only (incompatible with address sanitizer)
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=thread -fsanitize=undefined $(BASEFLAGS)
        SANLIBS := -ltsan
     else
        # address sanitizer
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address $(BASEFLAGS)
        SANLIBS := -lasan
     endif
     ifeq "$(GCCVER5)" "1"
       CXXFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CXXFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
       #CXXFLAGS += -fcheck-pointer-bounds -mmpx
     endif
     CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG -fno-common -fstack-protector
     LIBS := ${SANLIBS} -lubsan -ldl ${LIBS}
  else
     #just plain debug build
     DEBUG_BUILD=1
     CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
     ifneq (, $(findstring darwin, $(DMACH)))
        CXXFLAGS += -gdwarf-3
     endif
     CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
  endif
endif

ifdef RELEASE_BUILD
 ifneq (,$(findstring static, $(MAKECMDGOALS)))
    STATIC_CLIB=1
 endif
endif

ifdef STATIC_CLIB
 LDFLAGS += -static-libgcc -static-libstdc++
endif

ifdef DEBUG_BUILD
  #$(warning Building DEBUG version of tiebrush.. )
  DBG_WARN=@echo
  DBG_WARN+='WARNING: built DEBUG version, use "make clean release" for a faster version of the program.'
endif

OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o ./GSam.o ./tmerge.o

#OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o ./GSam.o \
# ${GDIR}/gdna.o ${GDIR}/codons.o ${GDIR}/GFastaIndex.o ${GDIR}/GFaSeqGet.o

ifneq (,$(filter %memtrace %memusage %memuse, $(MAKECMDGOALS)))
    CXXFLAGS += -DGMEMTRACE
    OBJS += ${GDIR}/proc_mem.o
endif

#ifndef NOTHREADS
# OBJS += ${GDIR}/GThreads.o 
#endif

%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

# OBJS += rlink.o tablemaker.o tmerge.o

all release static debug: tiebrush${EXE}
memcheck memdebug tsan tcheck thrcheck: tiebrush${EXE}
memuse memusage memtrace: tiebrush${EXE}
#nothreads: tiebrush${EXE}

GSam.o : GSam.h
tiebrush.o : GSam.h tmerge.h
tmerge.o : tmerge.h
#${BAM}/libhts.a: 
#	cd ${BAM} && make lib
tiebrush: $(OBJS) tiebrush.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}
test demo tests: tiebrush${EXE}
	@./run_tests.sh
.PHONY : clean cleanall cleanAll allclean

# target for removing all object files

#	echo $(PATH)
clean:
	${RM} tiebrush${EXE} tiebrush.o* $(OBJS)
	${RM} core.*
allclean cleanAll cleanall:
	cd ${BAM} && make clean
	${RM} tiebrush${EXE} tiebrush.o* $(OBJS)
	${RM} core.*
