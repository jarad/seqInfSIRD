SIM_OBJECTS = sims.RData

ANALYSIS_OBJECTS = PL.RData LW90.RData LW95.RData LW99.RData

QUANTILE_OBJECTS = PLq.RData LW90q.RData LW95q.RData LW99q.RData

OBJECTS = $(SIM_OBJECTS) $(ANALYSIS_OBJECTS) $(QUANTILE_OBJECTS)

RCMD = R CMD BATCH --vanilla



all: $(OBJECTS) 


sims: $(SIM_OBJECTS)

analysis: $(ANALYSIS_OBJECTS)

quantiles: $(QUANTILE_OBJECTS)



%.RData: %.R
	$(RCMD) $<
	git pull origin master
	git add $@
	-git commit -m "make: updated $@"
	-git push

clean:
	rm *.Rout
