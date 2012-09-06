
SOURCES= SckmSwarm.c SckmParticle.c Sckm.c gillespie.c utility.c pl-utility.c pl-discrete.c resample.c
OBJECTS=$(SOURCES:.c=.o)
HEADERS=$(SOURCES:.c=.h)
SHARED_OBJECT=tlpl.so

$(SHARED_OBJECT): $(SOURCES) $(HEADERS)
	R CMD SHLIB $(SOURCES) -o $(SHARED_OBJECT)

clean: 
	rm $(OBJECTS) $(SHARED_OBJECT)


