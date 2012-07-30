
TARGETS = SIRDsims-mcmc-density.RData \
          SIRDsims.RData \
          SIRDsims-smcLW-density.RData \
          SIRDsims-smcPL-density.RData \
          SIRDsims-smcSV-density.RData \
          SIRDsims-mcmc-quantiles.csv \
          SIRDsims-smcLW-quantiles.csv \
          SIRDsims-smcPL-quantiles.csv \
          SIRDsims-smcSV-quantiles.csv



method-comparison.pdf: $(TARGETS)
	R CMD BATCH SIRDsims-plots.R

SIRDsims-mcmc-density.RData: SIRDsims-mcmc-script.RData SIRDsims-mcmc-plotData.R
	R CMD BATCH SIRDsims-mcmc-plotData.R
	git pull origin master
	git add SIRDsims-mcmc-density.RData SIRDsims-mcmc-quantiles.csv
	git commit -m "make:updated SIRDsims-mcmc-sequential.RData and SIRDsims-mcmc-quantiles.csv"
	git push

SIRDsims-mcmc-sequential.RData: SIRDsims-mcmc-script.RData SIRDsims-mcmc-sequential.R
	R CMD BATCH SIRDsims-mcmc-sequential.R
	git pull origin master
	git add SIRDsims-mcmc-sequential.RData
	git commit -m "make:updated SIRDsims-mcmc-sequential.RData"
	git push

SIRDsims-mcmc-script.RData: SIRDsims.RData SIRDsims-mcmc-script.R
	R CMD BATCH SIRDsims-mcmc-script.R
	git pull origin master
	git add SIRDsims-mcmc-script.RData
	git commit -m "make: updated SIRDsims-mcmc-script.RData"
	git push

SIRDsims.RData: SIRDsims.R
	R CMD BATCH SIRDsims.R
	git pull origin master
	git add SIRDsims.RData
	git commit -m "make: updated SIRDsims.RData"
	git push


