
OBJECTS = SIR-sims.csv \
          SIR-BS.csv \
          SIR-BSq.csv \
          SIR-PL.csv \
          SIR-PLq.csv \
          SIR-SV.csv \
          SIR-SVq.csv \
          SIR-LW90.csv \
          SIR-LW90q.csv \
          SIR-LW95.csv \
          SIR-LW95q.csv \
          SIR-LW99.csv \
          SIR-LW99q.csv

RCMD = R CMD BATCH --vanilla



all: $(OBJECTS) 

%.csv : %.R SIR-sims.csv
	$(RCMD) $<
	git pull origin master
	git add $@
	git commit -m "make: updated $@"
	git push

SIR-sims.csv: SIR-sims.R
	$(RCMD) SIR-sims.R
	git pull origin master
	git add SIR-sims.csv
	git commit -m "make: updated SIR-sims.csv"
	git push

clean:
	rm *.Rout
