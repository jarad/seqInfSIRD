
SOURCES=simulate.c utility.c pl-utility.c pl-discrete.c
OBJECTS=$(SOURCES:.c=.o)
SHARED_OBJECT=tlpl.so

$(SHARED_OBJECT): $(SOURCES)
	R CMD SHLIB $(SOURCES) -o $(SHARED_OBJECT)

clean: 
	rm $(OBJECTS) $(SHARED_OBJECT)


