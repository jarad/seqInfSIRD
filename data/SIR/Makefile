
OBJECTS = sims.RData \
          PL.csv \
          PLq.csv \
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
	$(RCMD) sims.R
	git pull origin master
	git add sims.RData
	-git commit -m "make: updated sims.RData"
	-git push

clean:
	rm *.Rout
