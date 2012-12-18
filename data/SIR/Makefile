
TARGETS = SIR-sims.csv \
          SIR-PL.csv \
          SIR-PLq.csv \
          SIR-SV.csv \
          SIR-SVq.csv \
          SIR-LW99.csv \
          SIR-LW99q.csv \
          SIR-LW95.csv \
          SIR-LW95q.csv \
          SIR-LW90.csv \
          SIR-LW90q.csv \
          SIR

XX = PL SV LW99 LW95 LW90 BS

RCMD = R CMD BATCH --vanilla



SIR-$(XX)q.csv: SIR-$(XX)q.R SIR-$(XX).csv
	$(RCMD) SIR-$(XX)q.R 
	git pull origin master
	git add SIR-$(XX)q.csv
	git commit -m "make: updated SIR-$(XX)q.csv"
	git push



SIR-$(XX).csv: SIR-$(XX).R SIR-sims.csv
	$(RCMD) SIR-$(XX).R 
	git pull origin master
	git add SIR-$(XX).csv
	git commit -m "make: updated SIR-$(XX).csv"
	git push


SIR-sims.csv: SIR-sims.R
	$(RCMD) SIR-sims.R
	git pull origin master
	git add SIR-sims.csv
	git commit -m "make: updated SIR-sims.csv"
	git push


