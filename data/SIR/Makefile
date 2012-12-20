
OBJECTS = sims.RData \
          PL.RData \
          LW90.RData \
          LW95.RData \
          LW99.RData \

RCMD = R CMD BATCH --vanilla



all: $(OBJECTS) 

%.csv : %.R sims.csv
	$(RCMD) $<
	git pull origin master
	git add $@
	-git commit -m "make: updated $@"
	-git push

%.RData: %.R
	$(RCMD) $<
	git pull origin master
	git add $@
	-git commit -m "make: updated $@"
	-git push

clean:
	rm *.Rout
