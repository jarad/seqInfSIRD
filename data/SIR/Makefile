
OBJECTS = sims.RData \
          PL.RData \
          PLq.RData \
          SV.csv \
          SVq.csv \
          LW90.RData \
          LW90q.csv \
          LW95.RData \
          LW95q.csv \
          LW99.RData \
          LW99q.csv

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
