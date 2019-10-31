Task 5 Reference point comparisons (across candidate methods) 

Once reference points have been identified, their performance should be evaluated through simple management strategy evaluations. A set of appropriate stock assessment models
will be fit to the available data (Task 1 and 2). Methods used will reflect the available
data on a stock-by-stock basis (Task 2) and also methods available or developed in
Tasks 3 and 4. Performance diagnostics (e.g., residual inspection, retrospective patterns) will be run for each assessment model fit.

We will provide a review of reference points and indicators, based on a variety of assumptions, e.g. MSY based on biomass dynamic stock assessment models, and indicators such as Lopt ,L50 and Lmega for length-based methods. Then compare these using simulation, e.g. cross-testing and simulations based on FLife. This will allow the power of the methods to detect whether a stock has achieved its targets and avoid its limits. Based on this screening process a set of candidate reference points and assessment methods will be proposed for MSE.

MSE will include a Value-of-infomation analysis, where the benefits of collecting better data and new infomation will be evaluated.

Basically see the slides on next steps in the MI presentation

### **Case studies**
+ Summarise current knowledge and infomation
+ Evaluate alternative management plans and data collection options

### Evaluate **Methods**
+ Length
+ Abundance Indices
+ Catch

### Evaluate **Combined Rules** of the form $C_{t+1}=C_{t-1}rfb$

### Compare to **Category 1** rules

# Summary
## Summary
### **Management Procedures**
Consider the data along with the assumptions, estimator and management and feedback between the various elements

### **Diagnostics**
Compare across species, stocks, fisheries and advice rules

### **Value of infomation**
Models are cheap, but data is expensive

### **Risk**
An uncertainty that matters, what matters are management objectives


## R Packages

+ FLR
    + FLife
    + mydas
    + mpb
+ R Packages
    + LBSPR
    + MLZ
    + ...


+ Life histories were reviewed by Hans, so should base the OMs and priors on those.
+ Surveys do not give estimates of absolute abundance so we do as is usual practice, use relative indices.
+ Simon used the Crfb rule, in his paper where r is the relative abudance index (r) and f can use  proxy reference points for F and/or L. While b can use proxy reference points for biomass 

You need to

+ Perform lenth based assessments for all the case study stocks, using MLZ and LBSPR. This needs to include goodness of fit tests and sensitivity analyses. I.e. what if the input values are wrong?

+ Run the db-sra for the case studies

+ Run the biomass dynamic stock assessment, i.e. db-sra but with an index of abundnace if possible.

Once you have done this you need to run cross tests and set up the MSEs, for a variety of MPs for the MSEs which should be of the form Crfb, where 1st the individual rules are tested and the benefits of combnining them are evalauted.



Task 5 of the MYDAS project <https://github.com/laurieKell/mydas> requires a review of data avialble on length composition, mean size, catch and indices of abudance (commercial or survey).
Data are available from the MI  or from collective sources such as DATRAS <http://www.ices.dk/marine-data/data-portals/Pages/DATRAS.aspx> 
for survey data and the STECF <https://stecf.jrc.ec.europa.eu/data-dissemination> for commercial data. 


## Work plan

+	Condition the operating model using life histories (from table above) from brill, turbot, skate and Pollack as per the Mydas project (see the quick start). Will have to use MI estimates for sprat.  Use estimates of length (brill, sprat, skate). Obtain length frequencies for turbot and pollack from MI.

The life histories were reviewed by Hans, using data from the MI databases.

+	With further investigation look at predicted abundance estimates for brill, sprat (although on high side) and skate (abundance was for all species of skates and rays, will have to make consistent with commercial species... maybe use thornback only). Pollack and turbot catches are greater than abundance so cannot use these.

If the absolute abundance estimates are wrong relative indices can be used, as normal practice.  

+	With such good time series of survey length data for brill, sprat and skate it would good to sim test the mean size by using the mlz package <https://cran.r-project.org/web/packages/MLZ/vignettes/MLZ.html#introduction>. This would give an estimate of Z (assuming constant M) and thus look at changes in F. 

MLZ uses mean size, it is unlikely to provide good estimates of Z for sppart due to high M and high recruitment variability, the mean size is therefore determined by recruitment not F.

+	It may also be worth pursuing the LBSR for length-based composition for the above  (and compare) with life-history parameters M/K ratio and Linf to estimate F/M and F/FMSY.  With SPR being the biological reference point. <https://cran.r-project.org/web/packages/LBSPR/vignettes/LBSPR.html>

LBSPR is better at estimating F.

+	For turbot and pollack it maybe beneficially to use and empirical approach  (as survey info pretty poor) based on the commercial cpue (biomass index from MI lengths or estimated by weights) such as the ICES 2/3 rule (as per Simon Fischer).  

The surveys could be used in the 

For the observation model commercial cpues could be used.

Further: Need to discuss with MI what data we can get for Lobsters and Razors in terms of lengths.  We have 1 set of priors from Iyves.
