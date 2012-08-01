
TARGETS = SIRDsims-mcmc-density.RData \
          SIRDsims.RData \
          SIRDsims-smcLW-density.RData \
          SIRDsims-smcPL-density.RData \
          SIRDsims-smcSV-density.RData \
          SIRDsims-mcmc-quantiles.csv \
          SIRDsims-smcLW-quantiles.csv \
          SIRDsims-smcPL-quantiles.csv \
          SIRDsims-smcSV-quantiles.csv

RCMD = R CMD BATCH --vanilla

method-comparison.pdf: $(TARGETS)
	$(RCMD) SIRDsims-plots.R
	git pull origin master
	git add method-comparison.pdf
	git commit -m "make: updated method-comparison.pdf"
	git push

SIRDsims-mcmc-density.RData: SIRDsims-mcmc-script.RData SIRDsis-mcmc-sequential.RData SIRDsims-mcmc-plotData.R
	$(RCMD) SIRDsims-mcmc-plotData.R
	git pull origin master
	git add SIRDsims-mcmc-density.RData SIRDsims-mcmc-quantiles.csv
	git commit -m "make: updated SIRDsims-mcmc-density.RData and SIRDsims-mcmc-quantiles.csv"
	git push

SIRDsims-mcmc-sequential.RData: SIRDsims-mcmc-script.RData SIRDsims-mcmc-sequential.R
	$(RCMD) SIRDsims-mcmc-sequential.R

SIRDsims-mcmc-script.RData: SIRDsims.RData SIRDsims-mcmc-script.R
	$(RCMD) SIRDsims-mcmc-script.R

SIRDsims.RData: SIRDsims.R
	$(RCMD) SIRDsims.R
	git pull origin master
	git add SIRDsims.RData
	git commit -m "make: updated SIRDsims.RData"
	git push


