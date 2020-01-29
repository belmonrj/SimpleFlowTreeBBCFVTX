RFLAGS = `root-config --cflags --libs`

PROGRAMS = cumulants cumulants_pAu flattening minbias qsimple weight

all: $(PROGRAMS)

clean:
	rm -f $(PROGRAMS)

cumulants: cumulants.cpp
	g++ -o cumulants cumulants.cpp -Wall $(RFLAGS)

cumulants_pAu: cumulants_pAu.cpp
	g++ -o cumulants_pAu cumulants_pAu.cpp -Wall $(RFLAGS)

flattening: flattening.cpp
	g++ -o flattening flattening.cpp -Wall $(RFLAGS)

minbias: minbias.cpp
	g++ -o minbias minbias.cpp -Wall $(RFLAGS)

qsimple: qsimple.cpp
	g++ -o qsimple qsimple.cpp -Wall $(RFLAGS)

weight: weight.cpp
	g++ -o weight weight.cpp -Wall $(RFLAGS)

