SIM_OBJECTS = sims.RData sims.pdf

LWdelta = 90 95 99



ANALYSIS_OBJECTS = PL.RData LW90.RData LW95.RData LW99.RData

QUANT_OBJECTS = PLq.RData LW95q.RData LW99q.RData

SUM_OBJECTS = PLsum.csv LW90sum.csv LW95sum.csv LW99sum.csv example-plots.pdf

SETTINGS = settings.R

FILTER_SOURCE = $(SETTINGS) sims.RData

LW_SOURCE = $(FILTER_SOURCE) LW.R liu_west.r

LW_QUANT_SOURCE = $(SETTINGS) LWq.R 

SUM_SOURCE = summary.R

OBJECTS = $(SIM_OBJECTS) $(ANALYSIS_OBJECTS) $(QUANT_OBJECTS) $(SUM_OBJECTS)


RCMD = -R CMD BATCH --vanilla

##################################################################

all: $(OBJECTS) 


sims: $(SIM_OBJECTS)

analysis: $(ANALYSIS_OBJECTS)

quantiles: $(QUANT_OBJECTS)

summary: $(SUM_OBJECTS)

##################################################################


sims.RData: sims.R
	$(RCMD) $<

sims.pdf: sims-plot.R $(FILTER_SOURCE)
	$(RCMD) $<


PL.RData: PL.R $(FILTER_SOURCE)
	$(RCMD) $<

PLq.RData: PLq.R $(SETTINGS) PL.RData
	$(RCMD) $<

PLsum.csv: PLsum.R $(SUM_SOURCE) PLq.RData
	$(RCMD) $<



LW90.RData: LW90.R $(LW_SOURCE)
	$(RCMD) $<

LW95.RData: LW95.R $(LW_SOURCE)
	$(RCMD) $<

LW99.RData: LW99.R $(LW_SOURCE)
	$(RCMD) $<


LW90q.RData: LW90q.R $(LW_QUANT_SOURCE) LW90.RData
	$(RCMD) $<

LW95q.RData: LW95q.R $(LW_QUANT_SOURCE) LW95.RData
	$(RCMD) $<

LW99q.RData: LW99q.R $(LW_QUANT_SOURCE) LW99.RData
	$(RCMD) $<


LW90sum.csv: LW90sum.R $(SUM_SOURCE) LW90q.RData
	$(RCMD) $<

LW95sum.csv: LW95sum.R $(SUM_SOURCE) LW95q.RData
	$(RCMD) $<

LW99sum.csv: LW99sum.R $(SUM_SOURCE) LW99q.RData
	$(RCMD) $<


example-plots.pdf: example-plots.R $(QUANT_OBJECTS)
	$(RCMD) $<


clean:
	-rm *.Rout
