GDIR := ../gclib

## assumed htslib has been pulled from https://github.com/gpertea/htslib
HTSLIB := ../htslib
#my branch of htslib includes libdeflate:
LIBDEFLATE := ${HTSLIB}/xlibs/lib/libdeflate.a
LIBBZ2 := ${HTSLIB}/xlibs/lib/libbz2.a
LIBLZMA := ${HTSLIB}/xlibs/lib/liblzma.a
INCDIRS := -I. -I${GDIR} -I${HTSLIB}

BWLIB := ./bigwig

INCDIRS += -I${BWLIB}

#ifeq (${HLDIR}/libhts.a,$(wildcard ${HLDIR}/libhts.a))
# HTSLIBS := $(subst -lhts,${HLDIR}/libhts.a,${HTSLIBS})
#endif

CXX   := $(if $(CXX),$(CXX),g++)

BASEFLAGS := -Wall -Wextra ${INCDIRS} -fsigned-char -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -std=c++11 -fno-strict-aliasing -fno-exceptions -fno-rtti
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

LIBS := ${HTSLIB}/libhts.a ${LIBBZ2} ${LIBLZMA} ${LIBDEFLATE} -lz -lm -lpthread

#ifneq (,$(findstring nothreads,$(MAKECMDGOALS)))
# NOTHREADS=1
#endif

# Compiling for Windows with MinGW/MSYS2?
#ifneq ($(findstring -mingw,$(shell $(CC) -dumpmachine 2>/dev/null)),)
# LIBS += -lregex -lws2_32
#endif

#detect MinGW/MSYS (Windows environment)
ifneq (,$(findstring mingw,$(shell ${CXX} -dumpmachine)))
 WINDOWS=1
endif

# Misc. system commands
#ifdef WINDOWS
# RM = del /Q
#else
RM = rm -f
#endif

# File endings
ifdef WINDOWS
 EXE = .exe
 LIBS += -lregex -lws2_32
else
 EXE =
endif

DMACH := $(shell ${CXX} -dumpmachine)

ifneq (,$(filter %release %static %static-cpp, $(MAKECMDGOALS)))
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
 ifneq ($(findstring static,$(MAKECMDGOALS)),) 
  # static or static-cpp found
  ifneq ($(findstring static-cpp,$(MAKECMDGOALS)),) 
     #not a full static build, only c/c++ libs
     LDFLAGS := -static-libgcc -static-libstdc++ ${LDFLAGS}
  else
     #full static build
     LDFLAGS := -static -static-libgcc -static-libstdc++ ${LDFLAGS}
  endif
 endif
endif

#ifdef DEBUG_BUILD
#  DBG_WARN=@echo
#  DBG_WARN+='WARNING: built DEBUG version, use "make clean release" for a faster version of the program.'
#endif

OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o ./tmerge.o ./GSam.o
COVOBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o ./GSam.o

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

all: tiebrush tiecov
release static static-cpp debug: tiebrush tiecov
memcheck memdebug tsan tcheck thrcheck: tiebrush tiecov
memuse memusage memtrace: tiebrush tiecov

test: tiebrush tiecov
	cd test && ./run_tests.sh

valgrind: tiebrush tiecov
	cd test && ./run_valgrind.sh

GSam.o : GSam.h
tiebrush.o : GSam.h tmerge.h
tiecov.o : GSam.h
tmerge.o : tmerge.h
#${BAM}/libhts.a: 
#	cd ${BAM} && make lib

${BWLIB}/libBigWig.a:
	./build_bwlib.sh


${HTSLIB}/libhts.a: 
	cd ${HTSLIB} && ./build_lib.sh

tiebrush: ${HTSLIB}/libhts.a $(OBJS) tiebrush.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}
tiecov: ${BWLIB}/libBigWig.a ${HTSLIB}/libhts.a $(COVOBJS) tiecov.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${BWLIB}/libBigWig.a ${LIBS}
	@echo
	${DBG_WARN}

#test demo tests: tiebrush
#	@./run_tests.sh
.PHONY : clean cleanall cleanAll allclean test valgrind

# target for removing all object files

#	echo $(PATH)
clean:
	${RM} tiebrush${EXE} tiecov tiecov.o* tiebrush.o* $(OBJS)
	${RM} core.*
allclean cleanAll cleanall:
	cd ${BAM} && make clean
	${RM} tiebrush${EXE} tiebrush.o* $(OBJS)
	${RM} core.*
