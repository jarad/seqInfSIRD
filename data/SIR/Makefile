
OBJECTS = SIR-sims.RData \
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
	-git commit -m "make: updated $@"
	-git push

%.RData: %.R
	$(RCMD) SIR-sims.R
	git pull origin master
	git add SIR-sims.RData
	-git commit -m "make: updated SIR-sims.RData"
	-git push

clean:
	rm *.Rout
