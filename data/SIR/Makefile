
OBJECTS = sims.RData \
          PL.RData \
          PLq.RData \
          SV.csv \
          SVq.csv \
          LW90.csv \
          LW90q.csv \
          LW95.csv \
          LW95q.csv \
          LW99.csv \
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
