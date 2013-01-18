SIM_OBJECTS = sims.RData

ANALYSIS_OBJECTS = PL.RData LW90.RData LW95.RData LW99.RData

QUANT_OBJECTS = PLq.RData LW90q.RData LW95q.RData LW99q.RData



FILTER_SOURCE = filter-settings.R sims.RData

QUANT_SOURCE = quantile-settings.R 


OBJECTS = $(SIM_OBJECTS) $(ANALYSIS_OBJECTS) $(QUANT_OBJECTS)


RCMD = R CMD BATCH --vanilla

##################################################################

all: $(OBJECTS) 


sims: $(SIM_OBJECTS)

analysis: $(ANALYSIS_OBJECTS)

quantiles: $(QUANT_OBJECTS)

##################################################################


sims.RData: sims.R
	$(RCMD) $<

sims.pdf: sims-plot.R sims.RData
	$(RCMD) $<


PL.RData: PL.R $(FILTER_SOURCE)
	$(RCMD) $<

PLq.RData: PLq.R quantile-settings.R PL.RData
	$(RCMD) $<


LW90.RData: LW90.R LW.R $(FILTER_SOURCE)
	$(RCMD) $<

LW95.RData: LW95.R LW.R $(FILTER_SOURCE)
	$(RCMD) $<

LW99.RData: LW99.R LW.R $(FILTER_SOURCE)
	$(RCMD) $<


LW90q.RData: LW90q.R LWq.R quantile-settings.R LW90.RData
	$(RCMD) $<

LW95q.RData: LW95q.R LWq.R quantile-settings.R LW95.RData
	$(RCMD) $<

LW99q.RData: LW99q.R LWq.R quantile-settings.R LW95.RData
	$(RCMD) $<


clean:
	rm *.Rout
