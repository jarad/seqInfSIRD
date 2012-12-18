
TARGETS = SIR-sims.csv \
          SIR-smcLW-quantiles.csv \
          SIR-smcPL-quantiles.csv \
          SIR-smcSV-quantiles.csv

RCMD = R CMD BATCH --vanilla

method-comparison.pdf: $(TARGETS)
	$(RCMD) SIR-plots.R
	git pull origin master
	git add method-comparison.pdf
	git commit -m "make: updated method-comparison.pdf"
	git push

SIR-mcmc-density.RData: SIR-mcmc-script.RData SIR-mcmc-sequential.RData SIR-mcmc-plotData.R
	$(RCMD) SIR-mcmc-plotData.R
	git pull origin master
	git add SIR-mcmc-density.RData SIR-mcmc-quantiles.csv
	git commit -m "make: updated SIR-mcmc-density.RData and SIR-mcmc-quantiles.csv"
	git push

SIR-mcmc-sequential.RData: SIR-mcmc-script.RData SIR-mcmc-sequential.R
	$(RCMD) SIR-mcmc-sequential.R


SIR-sims.csv: SIR-sims.R
	$(RCMD) SIR-sims.R
	git pull origin master
	git add SIR-sims.csv
	git commit -m "make: updated SIR-sims.csv"
	git push


