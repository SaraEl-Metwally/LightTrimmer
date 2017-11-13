CXX = g++
CXXFLAGS = -g -Wall -O3 -std=c++0x 
LINKFLAGS = -lz -lpthread
DEBUG=
OBJECTS = KmerUtil.o Utility.o ReadsParsing.o LargeInt.o

all: LightTrimmer

LightTrimmer: clean main.o $(OBJECTS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)

main.o: main.cpp kseq.h KmersScanning.hpp MathUtil.hpp 
ReadsParsing.o: ReadsParsing.cpp ReadsParsing.hpp
KmerUtil.o: KmerUtil.cpp KmerUtil.hpp
Utility.o: Utility.cpp Utility.hpp
LargeInt.o: LargeInt.cpp LargeInt.hpp

k := 0$(k)
kmer_size_l32 := $(shell echo $(k)\<=32 | bc)
kmer_size_l64 := $(shell echo $(k)\<=64 | bc)
arch_stat := $(shell getconf LONG_BIT)
using_uint128 := 0
using_largeint := 0
ifeq ($(kmer_size_l32),0)
     ifeq ($(kmer_size_l64),1)
          ifeq ($(strip $(arch_stat)),64)
               CXXFLAGS += -Dkmercode_length=__uint128_t
               using_uint128 := 1
          endif
    endif
    ifeq ($(using_uint128),0)
          using_largeint := largeintlib
    endif
endif

ifeq ($(using_largeint),largeintlib)
    kmer_precision := $(shell echo \($(k)+31\)/32 | bc)
endif

ifneq ($(using_largeint),0)
    CXXFLAGS += -D$(using_largeint) -Dkmer_precision=$(kmer_precision)
endif

clean:
	rm -f *.o LightTrimmer
