# Andrew C. Eagar 2025

# Random slopes modeling and GLM slope coefficient modeling

library(tidyverse) # For data manipulation & graphing. Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686. 
library(lme4) # Linear modeling package. Bates, D., Mächler, M., Bolker, B. and Walker, S., 2014. Fitting linear mixed-effects models using lme4. arXiv preprint arXiv:1406.5823
library(MuMIn) # For AIC and R-squared values. Barton, K. and Barton, M.K., 2015. Package ‘mumin’. Version, 1(18), p.439.
library(lmerTest) # For p-values. Kuznetsova, A., Brockhoff, P.B. and Christensen, R.H., 2017. lmerTest package: tests in linear mixed effects models. Journal of statistical software, 82, pp.1-26.
library(emmeans) # For variable comparisons. Lenth, R. and Lenth, M.R., 2018. Package ‘lsmeans’. The American Statistician, 34(4), pp.216-221.

# Import data  & format for random slopes modeling

RM_demography <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/rm_count_demography.csv")

Site_BA <- read.csv("C:/Users/eagar/Dropbox/IBEEM/conspecific_density_dependence/data/L1/MNF_ba_radii_allCensuses.csv")%>%
					replace(is.na(.), 0)

Site_BA_quads <- Site_BA %>%
					mutate(quadrat = as.numeric(X) + 0.5) %>%
					mutate(Transect_quad_year_unifier = paste(TransectID, quadrat, Year)) %>%
					select(-c(1, 2, 4, 5, 6, 7, 63))

## Calculate conmycorrhizal, heterospecific, and heteromycorrhizal BA for RM
# AM Species: sugar maple (sm), American elm (ae), white ash (wa), red maple (rm), black cherry (bc), serviceberry (svb)
# EM Species: American beech (ab), basswood (bs), ironwood (iw), northern red oak (nro), black oak (bo), black oak/northern red oak hybrid (bonro), white oak (wo), hemlock (hem), white pine (wp), jack pine (jp), quaking aspen (qa), bigtooth aspen (bta)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - WITH rm
{Site_BA_quads$conMyc_5m = rowSums(Site_BA_quads[, c("sm_5m", "ae_5m", "wa_5m", 
                                        "rm_5m", "bc_5m", "svb_5m")],  na.rm = TRUE)
Site_BA_quads$conMyc_10m = rowSums(Site_BA_quads[, c("sm_10m", "ae_10m", "wa_10m", 
                                        "rm_10m", "bc_10m", "svb_10m")],  na.rm = TRUE)
Site_BA_quads$conMyc_15m = rowSums(Site_BA_quads[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "rm_15m", "bc_15m", "svb_15m")],  na.rm = TRUE)
## Conmycorrhizal BA in a 5-m, 10-m, and 15-m radius - NO rm
Site_BA_quads$conMycNS_5m = rowSums(Site_BA_quads[, c("sm_5m", "ae_5m", "wa_5m", 
                                         "bc_5m", "svb_5m")],  na.rm = TRUE)
Site_BA_quads$conMycNS_10m = rowSums(Site_BA_quads[, c("sm_10m", "ae_10m", "wa_10m", 
                                          "bc_10m", "svb_10m")],  na.rm = TRUE)
Site_BA_quads$conMycNS_15m = rowSums(Site_BA_quads[, c("sm_15m", "ae_15m", "wa_15m", 
                                         "bc_15m", "svb_15m")],  na.rm = TRUE)
## Heteromycorrhizal BA in a 5-m, 10-m, and 15-m radius
Site_BA_quads$hetMyc_5m = rowSums(Site_BA_quads[, c("ab_5m", "bs_5m", "qa_5m","bta_5m", 
                                        "nro_5m","iw_5m", "bo_5m","bonro_5m",
                                        "wo_5m", "hem_5m","wp_5m","jp_5m")],  na.rm = TRUE)
Site_BA_quads$hetMyc_10m = rowSums(Site_BA_quads[, c("ab_10m", "bs_10m", "qa_10m","bta_10m", 
                                          "nro_10m","iw_10m", "bo_10m","bonro_10m",
                                          "wo_10m", "hem_10m","wp_10m","jp_10m")],  na.rm = TRUE)
Site_BA_quads$hetMyc_15m = rowSums(Site_BA_quads[, c("ab_15m", "bs_15m", "qa_15m","bta_15m", 
                                         "nro_15m","iw_15m", "bo_15m","bonro_15m",
                                         "wo_15m", "hem_15m","wp_15m","jp_15m")],  na.rm = TRUE)
}
										
## Use an ifelse statement to add a column with the closest BA census year -a B. McNichol					
{ RM_demography$baYear = ifelse(RM_demography$First.Year == "1999", "1999", 
                        ifelse(RM_demography$First.Year == "2000", "1999",
                        ifelse(RM_demography$First.Year == "2001", "1999",
                        ifelse(RM_demography$First.Year == "2002", "1999",
                        ifelse(RM_demography$First.Year == "2003", "2005",
                        ifelse(RM_demography$First.Year == "2004", "2005",
                        ifelse(RM_demography$First.Year == "2005", "2005",
                        ifelse(RM_demography$First.Year == "2006", "2005",
                        ifelse(RM_demography$First.Year == "2007", "2007",
                        ifelse(RM_demography$First.Year == "2008", "2007",
                        ifelse(RM_demography$First.Year == "2009", "2009",
                        ifelse(RM_demography$First.Year == "2010", "2009",
                        ifelse(RM_demography$First.Year == "2011", "2011",
                        ifelse(RM_demography$First.Year == "2012", "2011",
                        ifelse(RM_demography$First.Year == "2013", "2013",
                        ifelse(RM_demography$First.Year == "2014", "2013",
                        ifelse(RM_demography$First.Year == "2015", "2016",
                        ifelse(RM_demography$First.Year == "2016", "2016",
                        ifelse(RM_demography$First.Year == "2017", "2016",
                        ifelse(RM_demography$First.Year == "2018", "2019",
                        ifelse(RM_demography$First.Year == "2019", "2019",
                        ifelse(RM_demography$First.Year == "2020", "2019",
                        ifelse(RM_demography$First.Year == "2021", "2022",
                        ifelse(RM_demography$First.Year == "2022", "2022",
                        ifelse(RM_demography$First.Year == "2023", "2022",
                               "2022")))))))))))))))))))))))))
} 
							   
RM_demography_spread <- RM_demography %>%
				spread(key = "Time.Since.First.Year", value = "Count")	
					
RM_demography_spread$Starting_counts <- RM_demography_spread$'0' 

RM_demography_BA <- RM_demography_spread %>%
				gather(c(6:9), key = "Time.Since.First.Year", value = "Count") %>%
				mutate(Transect_quad_year_unifier = paste(TransectID, quadrat, baYear)) %>%
				left_join(Site_BA_quads, by = "Transect_quad_year_unifier")

# Analysis of RM density alone:

RM_rand_slopes_glm_rm <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + rm_5m + (rm_5m | SiteID + First.Year), family = binomial(link = "logit"), data = RM_demography_BA) # Model includes a random slope for tree BA

	summary(RM_rand_slopes_glm_rm)
	
RM_glm_rm <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + rm_5m + (1 | SiteID + First.Year), family = binomial(link = "logit"), data = RM_demography_BA)

	anova(RM_rand_slopes_glm_rm, RM_glm_rm)
	r.squaredGLMM(RM_rand_slopes_glm_rm)
	r.squaredGLMM(RM_glm_rm)
	
{# Output. In addition to the fixed effects, we want to pay attention to the magnitude of variance assigned to the random slope term in the model.

Generalized linear mixed model fit by maximum likelihood (Laplace
  Approximation) [glmerMod]
 Family: binomial  ( logit )
Formula: 
cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) +  
    rm_5m + (rm_5m | SiteID + First.Year)
   Data: RM_demography_BA

     AIC      BIC   logLik deviance df.resid 
168936.5 169019.3 -84459.3 168918.5    72455 

Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-14.9419  -0.7068   0.2602   0.7784  13.0079 

Random effects:
 Groups     Name        Variance Std.Dev. Corr 
 First.Year (Intercept) 0.3398   0.5829        
            rm_5m       1.0590   1.0291   -0.60
 SiteID     (Intercept) 0.2534   0.5034        
            rm_5m       8.3630   2.8919   0.37 # Large degree of site variation in slopes
Number of obs: 72464, groups:  First.Year, 20; SiteID, 12

Fixed effects:
                                  Estimate Std. Error  z value Pr(>|z|)    
(Intercept)                        2.04912    0.19418   10.553   <2e-16 ***
as.numeric(Time.Since.First.Year) -1.36928    0.00512 -267.461   <2e-16 ***
rm_5m                             -1.89625    0.78712   -2.409    0.016 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) a.(T.S
a.(T.S.F.Y) -0.038       
rm_5m        0.112  0.003

anova(RM_rand_slopes_glm_rm, RM_glm_rm) # Log liklihood

RM_glm_rm: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + rm_5m + (1 | SiteID + First.Year)
RM_rand_slopes_glm_rm: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + rm_5m + (rm_5m | SiteID + First.Year)
                      npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)
RM_glm_rm                5 169490 169536 -84740   169480                     
RM_rand_slopes_glm_rm    9 168937 169019 -84459   168919 561.74  4  < 2.2e-16

r.squaredGLMM(RM_rand_slopes_glm_rm)
r.squaredGLMM(RM_glm_rm)

                  R2m       R2c
theoretical 0.6078882 0.7890296
delta       0.5813914 0.7546372

vs.
                  R2m       R2c
theoretical 0.6267793 0.7789485
delta       0.5982123 0.7434460

}
					
# Analysis of CMDD (includes conspecifics)

RM_rand_slopes_glm_con_myc <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + conMyc_5m + (conMyc_5m | SiteID + First.Year), family = binomial(link = "logit"), data = RM_demography_BA) # Model includes a random slope for tree BA

	summary(RM_rand_slopes_glm_con_myc)
	
RM_glm_con_myc <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + conMyc_5m + (1 | SiteID + First.Year), family = binomial(link = "logit"), data = RM_demography_BA)

	anova(RM_rand_slopes_glm_con_myc, RM_glm_con_myc)
	r.squaredGLMM(RM_rand_slopes_glm_con_myc)
	r.squaredGLMM(RM_glm_con_myc)
	
{# Output. In addition to the fixed effects, we want to pay attention to the magnitude of variance assigned to the random slope term in the model.

Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: binomial  ( logit )
Formula: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) +  
    conMyc_5m + (conMyc_5m | SiteID + First.Year)
   Data: RM_demography_BA

     AIC      BIC   logLik deviance df.resid 
168502.1 168584.8 -84242.1 168484.1    72455 

Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-14.8380  -0.7041   0.2573   0.7818  13.3185 

Random effects:
 Groups     Name        Variance Std.Dev. Corr 
 First.Year (Intercept) 0.3538   0.5948        
            conMyc_5m   1.5083   1.2281   -0.52
 SiteID     (Intercept) 0.2366   0.4864        
            conMyc_5m   0.7840   0.8855   0.01 
Number of obs: 72464, groups:  First.Year, 20; SiteID, 12

Fixed effects:
                                   Estimate Std. Error  z value Pr(>|z|)    
(Intercept)                        2.060045   0.194449   10.594   <2e-16 ***
as.numeric(Time.Since.First.Year) -1.371809   0.005129 -267.458   <2e-16 ***
conMyc_5m                         -0.848231   0.416172   -2.038   0.0415 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) a.(T.S
a.(T.S.F.Y) -0.039       
conMyc_5m   -0.242  0.004


anova(RM_rand_slopes_glm_con_myc, RM_glm_con_myc) # Log liklihood

RM_glm_con_myc: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + conMyc_5m + (1 | SiteID + First.Year)
RM_rand_slopes_glm_con_myc: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + conMyc_5m + (conMyc_5m | SiteID + First.Year)
                           npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
RM_glm_con_myc                5 169554 169600 -84772   169544                         
RM_rand_slopes_glm_con_myc    9 168502 168585 -84242   168484 1059.7  4  < 2.2e-16

r.squaredGLMM(RM_rand_slopes_glm_con_myc)
r.squaredGLMM(RM_glm_con_myc)

                  R2m       R2c
theoretical 0.6278981 0.7809892
delta       0.5994603 0.7456178

vs.
                  R2m       R2c
theoretical 0.6303018 0.7775525
delta       0.6014011 0.7419001

}

# Analysis of CMDD minus conspecifics (this is heterospecific CMDD)

RM_rand_slopes_glm_con_myc_no_rm <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + conMycNS_5m + (conMycNS_5m | SiteID + First.Year), family = binomial(link = "logit"), data = RM_demography_BA) # Model includes a random slope for tree BA

	summary(RM_rand_slopes_glm_con_myc_no_rm)
	
RM_glm_con_myc_no_rm <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + conMycNS_5m + (1 | SiteID + First.Year), family = binomial(link = "logit"), data = RM_demography_BA)

	anova(RM_rand_slopes_glm_con_myc_no_rm, RM_glm_con_myc_no_rm)
	r.squaredGLMM(RM_rand_slopes_glm_con_myc_no_rm)
	r.squaredGLMM(RM_glm_con_myc_no_rm)

{# Output. In addition to the fixed effects, we want to pay attention to the magnitude of variance assigned to the random slope term in the model.

Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: binomial  ( logit )
Formula: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) +  
    conMycNS_5m + (conMycNS_5m | SiteID + First.Year)
   Data: RM_demography_BA

     AIC      BIC   logLik deviance df.resid 
168473.8 168556.5 -84227.9 168455.8    72455 

Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-14.8366  -0.7010   0.2620   0.7893  13.2126 

Random effects:
 Groups     Name        Variance Std.Dev. Corr 
 First.Year (Intercept) 0.3331   0.5771        
            conMycNS_5m 7.7377   2.7817   -0.26
 SiteID     (Intercept) 0.2387   0.4886        
            conMycNS_5m 3.9736   1.9934   -0.75
Number of obs: 72464, groups:  First.Year, 20; SiteID, 12

Fixed effects:
                                   Estimate Std. Error  z value Pr(>|z|)    
(Intercept)                        2.015638   0.185169   10.885   <2e-16 ***
as.numeric(Time.Since.First.Year) -1.371876   0.005129 -267.469   <2e-16 ***
conMycNS_5m                       -1.877477   0.794754   -2.362   0.0182 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) a.(T.S
a.(T.S.F.Y) -0.039       
conMycNS_5m -0.340  0.001

anova(RM_rand_slopes_glm_con_myc_no_rm, RM_glm_con_myc_no_rm)

RM_glm_con_myc_no_rm: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + conMycNS_5m + (1 | SiteID + First.Year)
RM_rand_slopes_glm_con_myc_no_rm: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + conMycNS_5m + (conMycNS_5m | SiteID + First.Year)
                                 npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
RM_glm_con_myc_no_rm                5 169773 169819 -84881   169763                         
RM_rand_slopes_glm_con_myc_no_rm    9 168474 168557 -84228   168456 1306.9  4  < 2.2e-16

	r.squaredGLMM(RM_rand_slopes_glm_con_myc_no_rm)
	r.squaredGLMM(RM_glm_con_myc_no_rm)

                  R2m       R2c
theoretical 0.6301049 0.7808869
delta       0.6013846 0.7452940

vs
                  R2m       R2c
theoretical 0.6253120 0.7784285
delta       0.5967479 0.7428701

}

# Analysis of HMDD.

RM_rand_slopes_glm_het_myc <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + hetMyc_5m + (hetMyc_5m | SiteID + First.Year), family = binomial(link = "logit"), data = RM_demography_BA) # Model includes a random slope for tree BA

	summary(RM_rand_slopes_glm_het_myc)
	
RM_glm_het_myc <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + hetMyc_5m + (1 | SiteID + First.Year), family = binomial(link = "logit"), data = RM_demography_BA)

	anova(RM_rand_slopes_glm_het_myc, RM_glm_het_myc)
	r.squaredGLMM(RM_rand_slopes_glm_het_myc)
	r.squaredGLMM(RM_glm_het_myc)

{# Output. In addition to the fixed effects, we want to pay attention to the magnitude of variance assigned to the random slope term in the model.

Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: binomial  ( logit )
Formula: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) +  
    hetMyc_5m + (hetMyc_5m | SiteID + First.Year)
   Data: RM_demography_BA

     AIC      BIC   logLik deviance df.resid 
169388.7 169471.4 -84685.3 169370.7    72455 

Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-14.8278  -0.7084   0.2627   0.7864  13.3545 

Random effects:
 Groups     Name        Variance Std.Dev. Corr 
 First.Year (Intercept) 0.2874   0.5361        
            hetMyc_5m   0.1441   0.3796   0.21 
 SiteID     (Intercept) 0.2630   0.5129        
            hetMyc_5m   0.8749   0.9354   -0.17
Number of obs: 72464, groups:  First.Year, 20; SiteID, 12

Fixed effects:
                                  Estimate Std. Error  z value Pr(>|z|)    
(Intercept)                        2.00509    0.19310   10.384   <2e-16 ***
as.numeric(Time.Since.First.Year) -1.36654    0.00511 -267.437   <2e-16 ***
hetMyc_5m                          0.14481    0.30258    0.479    0.632    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) a.(T.S
a.(T.S.F.Y) -0.039       
hetMyc_5m   -0.120 -0.001

anova(RM_rand_slopes_glm_het_myc, RM_glm_het_myc)

RM_glm_het_myc: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + hetMyc_5m + (1 | SiteID + First.Year)
RM_rand_slopes_glm_het_myc: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + hetMyc_5m + (hetMyc_5m | SiteID + First.Year)
                           npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
RM_glm_het_myc                5 169741 169787 -84865   169731                         
RM_rand_slopes_glm_het_myc    9 169389 169471 -84685   169371 360.04  4  < 2.2e-16

                  R2m       R2c
theoretical 0.6185625 0.7817778
delta       0.5907279 0.7465988

vs
                  R2m       R2c
theoretical 0.6257149 0.7784184
delta       0.5971312 0.7428589


}


# Analysis of HMDD interacting with CMDD - do we do this with, or without, conspecific species? Both done here, for now.

hetMyc_5m conMyc_5m or conMycNS_5m

AIC(RM_rand_slopes_glm_rm, RM_rand_slopes_glm_con_myc, RM_rand_slopes_glm_con_myc_no_rm) # Use these model comparisons to inform which CMDD term should be included in the het + con modeling?

                                 df      AIC
RM_rand_slopes_glm_rm             6 169234.7
RM_rand_slopes_glm_con_myc        9 168502.1
RM_rand_slopes_glm_con_myc_no_rm  9 168473.8

# With conspecifics

RM_rand_slopes_glm_het_con_myc <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + hetMyc_5m*conMyc_5m + (hetMyc_5m + conMyc_5m | SiteID + First.Year), family = binomial(link = "logit"), data = RM_demography_BA) # Model includes a random slope for tree BA

	summary(RM_rand_slopes_glm_het_con_myc)
	
RM_glm_het_con_myc <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + hetMyc_5m*conMyc_5m + (1 | SiteID + First.Year), family = binomial(link = "logit"), data = RM_demography_BA)

	anova(RM_rand_slopes_glm_het_con_myc, RM_glm_het_con_myc)
	r.squaredGLMM(RM_rand_slopes_glm_het_con_myc)
	r.squaredGLMM(RM_glm_het_con_myc)

{# Output. In addition to the fixed effects, we want to pay attention to the magnitude of variance assigned to the random slope term in the model.

Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: binomial  ( logit )
Formula: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) +  
    hetMyc_5m * conMyc_5m + (hetMyc_5m + conMyc_5m | SiteID +      First.Year)
   Data: RM_demography_BA

     AIC      BIC   logLik deviance df.resid 
168271.8 168428.0 -84118.9 168237.8    72447 

Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-14.8244  -0.7056   0.2558   0.7794  13.5645 

Random effects:
 Groups     Name        Variance Std.Dev. Corr       
 First.Year (Intercept) 0.3238   0.5691              
            hetMyc_5m   0.1298   0.3603    0.19      
            conMyc_5m   1.3955   1.1813   -0.52  0.13
 SiteID     (Intercept) 0.2373   0.4872              
            hetMyc_5m   0.6368   0.7980   -0.20      
            conMyc_5m   0.9025   0.9500   -0.01  0.38
Number of obs: 72464, groups:  First.Year, 20; SiteID, 12

Fixed effects:
                                   Estimate Std. Error  z value Pr(>|z|)    
(Intercept)                        2.082153   0.192625   10.809  < 2e-16 ***
as.numeric(Time.Since.First.Year) -1.373715   0.005136 -267.451  < 2e-16 ***
hetMyc_5m                          0.077212   0.266555    0.290  0.77207    
conMyc_5m                         -0.801852   0.420945   -1.905  0.05680 .  
hetMyc_5m:conMyc_5m               -0.845525   0.300862   -2.810  0.00495 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) a.(T.S htMy_5 cnMy_5
a.(T.S.F.Y) -0.040                     
hetMyc_5m   -0.131  0.000              
conMyc_5m   -0.242  0.004  0.281       
htMyc_5:M_5  0.015  0.008 -0.081 -0.113

anova(RM_rand_slopes_glm_het_con_myc, RM_glm_het_con_myc)

RM_glm_het_con_myc: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + hetMyc_5m * conMyc_5m + (1 | SiteID + First.Year)
RM_rand_slopes_glm_het_con_myc: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + hetMyc_5m * conMyc_5m + (hetMyc_5m + conMyc_5m | SiteID + First.Year)
                               npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
RM_glm_het_con_myc                7 169524 169589 -84755   169510                         
RM_rand_slopes_glm_het_con_myc   17 168272 168428 -84119   168238 1272.7 10  < 2.2e-16

                  R2m       R2c
theoretical 0.6240013 0.7830960
delta       0.5960427 0.7480091

vs
                  R2m       R2c
theoretical 0.6300341 0.7777466
delta       0.6011698 0.7421150

}

# Without conpsecifics 

RM_rand_slopes_glm_het_con_myc_no_rm <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + hetMyc_5m*conMycNS_5m + (hetMyc_5m + conMycNS_5m | SiteID + First.Year), family = binomial(link = "logit"), data = RM_demography_BA) # Model includes a random slope for tree BA

	summary(RM_rand_slopes_glm_het_con_myc_no_rm)
	
RM_glm_het_con_myc_no_rm <- glmer(cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + hetMyc_5m*conMycNS_5m + (1 | SiteID + First.Year), family = binomial(link = "logit"), data = RM_demography_BA)

	anova(RM_rand_slopes_glm_het_con_myc_no_rm, RM_glm_het_con_myc_no_rm)
	r.squaredGLMM(RM_rand_slopes_glm_het_con_myc_no_rm)
	r.squaredGLMM(RM_glm_het_con_myc_no_rm)
	
{# Output. In addition to the fixed effects, we want to pay attention to the magnitude of variance assigned to the random slope term in the model.

Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: binomial  ( logit )
Formula: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) +  
    hetMyc_5m * conMycNS_5m + (hetMyc_5m + conMycNS_5m | SiteID +      First.Year)
   Data: RM_demography_BA

     AIC      BIC   logLik deviance df.resid 
168371.3 168527.5 -84168.6 168337.3    72447 

Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-14.9220  -0.7006   0.2606   0.7862  12.8641 

Random effects:
 Groups     Name        Variance Std.Dev. Corr       
 First.Year (Intercept) 0.30560  0.5528              
            hetMyc_5m   0.08648  0.2941    0.26      
            conMycNS_5m 7.62021  2.7605   -0.24  0.02
 SiteID     (Intercept) 0.24029  0.4902              
            hetMyc_5m   0.24366  0.4936   -0.22      
            conMycNS_5m 2.67846  1.6366   -0.71  0.45
Number of obs: 72464, groups:  First.Year, 20; SiteID, 12

Fixed effects:
                                   Estimate Std. Error  z value Pr(>|z|)    
(Intercept)                        2.024819   0.187152   10.819   <2e-16 ***
as.numeric(Time.Since.First.Year) -1.372917   0.005133 -267.454   <2e-16 ***
hetMyc_5m                          0.046919   0.176558    0.266   0.7904    
conMycNS_5m                       -1.757573   0.824020   -2.133   0.0329 *  
hetMyc_5m:conMycNS_5m              0.677155   0.494362    1.370   0.1708    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) a.(T.S htMy_5 cMNS_5
a.(T.S.F.Y) -0.040                     
hetMyc_5m   -0.099  0.000              
conMycNS_5m -0.346  0.003  0.213       
htM_5:MNS_5  0.014 -0.002 -0.128 -0.060

anova(RM_rand_slopes_glm_het_con_myc, RM_glm_het_con_myc)

RM_glm_het_con_myc_no_rm: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + hetMyc_5m * conMycNS_5m + (1 | SiteID + First.Year)
RM_rand_slopes_glm_het_con_myc_no_rm: cbind(Count, Starting_counts - Count) ~ as.numeric(Time.Since.First.Year) + hetMyc_5m * conMycNS_5m + (hetMyc_5m + conMycNS_5m | SiteID + First.Year)
                                     npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
RM_glm_het_con_myc_no_rm                7 169685 169749 -84835   169671                         
RM_rand_slopes_glm_het_con_myc_no_rm   17 168371 168527 -84169   168337 1333.4 10  < 2.2e-16

	r.squaredGLMM(RM_rand_slopes_glm_het_con_myc_no_rm)
	r.squaredGLMM(RM_glm_het_con_myc_no_rm)
	
                  R2m       R2c
theoretical 0.6293791 0.7812267
delta       0.6009037 0.7458812

vs 
                 R2m       R2c
theoretical 0.6245753 0.7789949
delta       0.5961144 0.7434974

}	

AIC(RM_rand_slopes_glm_rm, RM_rand_slopes_glm_con_myc, RM_rand_slopes_glm_con_myc_no_rm, RM_rand_slopes_glm_het_myc, RM_rand_slopes_glm_het_con_myc_no_rm)

                                     df      AIC
RM_rand_slopes_glm_rm                 6 169234.7
RM_rand_slopes_glm_con_myc            9 168502.1
RM_rand_slopes_glm_con_myc_no_rm      9 168473.8 # Best model!
RM_rand_slopes_glm_het_myc            9 169388.7
RM_rand_slopes_glm_het_con_myc_no_rm 17 168371.3 # This model has different random effects and probably shouldn't be compared?

{# Plot mortality response to RM densities - separate data into two plots; 1 that includes quadrats with 0 - 0.1 RM densities (low), and 1 that includes quadrats wtih >0.1 RM densities.

hist(RM_demography_BA$rm_5m)
hist(RM_demography_BA$conMycNS_5m)
hist(RM_demography_BA$hetMyc_5m)



theme_Scatter <- theme_classic() + # Classic theme with white background.
			theme(aspect.ratio = 1, plot.title = (element_text(face = "bold", size = 25, hjust = 0.5, vjust = -2, margin = margin(b = 25))), axis.title = (element_text(face = "bold", size = 20)), legend.title = element_text(face = "bold", size = 16), axis.text = element_text(size = 15, color = "black"), legend.margin = margin(0.5,0,0,0.5, unit = "cm"), text = element_text(size = 16), legend.background = element_rect(fill='transparent'), legend.box.background = element_rect(fill='transparent', linetype = 'blank'), legend.spacing.y = unit(0.01, 'cm')) # Adjust plot and legend margins, bold titles, and make other cosmetic changes. legend.position = c(0.12, 0.28),

# Plot 1: 0 - 0.1 RM densities

(			
RM_glm_rm_plot_0_rm <- RM_demography_BA %>%
					filter(Time.Since.First.Year != 0 & rm_5m <= 0.1) %>%
					group_by(First.Year, SiteID, quadrat) %>%
					mutate(Perc_mort = Count/Starting_counts) %>%
					ggplot(aes(x = as.numeric(Time.Since.First.Year), y = Perc_mort)) +
								# geom_line(aes(color = rm_5m, group = Transect_quad_year_unifier), linewidth = 0.1) +
								geom_smooth(aes(color = rm_5m, group = Transect_quad_year_unifier), linewidth = 0.1, method = "glm", method.args = list(family = "binomial"), se = FALSE) +
								theme_Scatter +
								scale_color_gradient2(low = "grey90", mid = "firebrick1", high = "firebrick4", midpoint = 0.4)
								
								)
								
# Plot 2: > 0.1 RM densities

(			
RM_glm_rm_plot <- RM_demography_BA %>%
					filter(Time.Since.First.Year != 0 & rm_5m > 0.1) %>%
					group_by(First.Year, SiteID, quadrat) %>%
					mutate(Perc_mort = Count/Starting_counts) %>%
					ggplot(aes(x = as.numeric(Time.Since.First.Year), y = Perc_mort)) +
								# geom_line(aes(color = rm_5m, group = Transect_quad_year_unifier), linewidth = 0.1) +
								geom_smooth(aes(color = rm_5m, group = Transect_quad_year_unifier), linewidth = 0.1, method = "glm", method.args = list(family = "binomial"), se = FALSE) +
								theme_Scatter +
								scale_color_gradient2(low = "grey90", mid = "firebrick1", high = "firebrick4", midpoint = 0.4)
								
								)								
								
# Plot 3: All densities, plotted by unique densities (fewer lines!)

(			
RM_glm_rm_plot_densities <- RM_demography_BA %>%
					filter(Time.Since.First.Year != 0) %>%
					group_by(First.Year, SiteID, quadrat) %>%
					mutate(Perc_mort = Count/Starting_counts) %>%
					ggplot(aes(x = as.numeric(Time.Since.First.Year), y = Perc_mort)) +
								# geom_line(aes(color = rm_5m, group = Transect_quad_year_unifier), linewidth = 0.1) +
								geom_smooth(aes(color = rm_5m, group = rm_5m), linewidth = 0.1, method = "glm", method.args = list(family = "binomial"), se = FALSE) +
								theme_Scatter +
								scale_color_gradient2(low = "grey90", mid = "salmon", high = "firebrick4", midpoint = 0.4)
								
								)								
								
# Plot 4: All densities, plotted by quadrats

(			
RM_glm_rm_plot_densities <- RM_demography_BA %>%
					filter(Time.Since.First.Year != 0) %>%
					group_by(First.Year, SiteID, quadrat) %>%
					mutate(Perc_mort = Count/Starting_counts) %>%
					ggplot(aes(x = as.numeric(Time.Since.First.Year), y = Perc_mort)) +
								# geom_line(aes(color = rm_5m, group = Transect_quad_year_unifier), linewidth = 0.1) +
								geom_smooth(aes(color = rm_5m, group = Transect_quad_year_unifier), linewidth = 0.1, method = "glm", method.args = list(family = "binomial"), se = FALSE) +
								theme_Scatter +
								scale_color_gradient2(low = "grey90", mid = "salmon", high = "firebrick4", midpoint = 0.4)
								
								)									
								
# Plot 5: All times plotted by density

(			
RM_glm_rm_plot_densities <- RM_demography_BA %>%
					filter(Time.Since.First.Year != 0) %>%
					group_by(First.Year, SiteID, quadrat) %>%
					mutate(Perc_mort = Count/Starting_counts, First.Year.Time = paste(Time.Since.First.Year, First.Year)) %>%
					ggplot(aes(x = as.numeric(rm_5m), y = Perc_mort)) +
								# geom_line(aes(color = rm_5m, group = Transect_quad_year_unifier), linewidth = 0.1) +
								geom_smooth(aes(color = Time.Since.First.Year, group = First.Year.Time), linewidth = 0.1, method = "glm", method.args = list(family = "binomial"), se = FALSE) +
								theme_Scatter
)
							
(			
RM_glm_rm_plot_densities <- RM_demography_BA %>%
					filter(Time.Since.First.Year != 0) %>%
					group_by(First.Year, SiteID, quadrat) %>%
					mutate(Perc_mort = Count/Starting_counts) %>%
					ggplot(aes(x = as.numeric(rm_5m), y = Perc_mort)) +
								# geom_line(aes(color = rm_5m, group = Transect_quad_year_unifier), linewidth = 0.1) +
								geom_smooth(aes(color = Time.Since.First.Year, group = Time.Since.First.Year), linewidth = 0.1, method = "glm", method.args = list(family = "binomial"), se = FALSE) +
								theme_Scatter
)								
								+ 
								geom_text(x = 46, y = 0.9, label = "R² = 67.4% (F, G, X)", size = 22) +
								geom_text(x = 18, y = 1, label = "a", size = 21) +
								scale_shape(labels = c("AM Focal Tree", "ECM Focal Tree")) +
								labs(title = "Soil", x = "% Relative ECM Basal Area", y = "Primary Saprotroph \nRelative Abundance") +
								scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0, 1)) +
								scale_x_continuous(breaks = seq(20,100,20), limits = c(18, 101)) +
								scale_color_manual(values = c("orange1", "springgreen3")) +
								stat_smooth(data = (Soil_Guild_tot[which(Soil_Guild_tot$Myc_type == "AM"),]), method = "glm", method.args = list(family = "binomial"), color = "orange1", se = FALSE, size = 6) +
								stat_smooth(data = (Soil_Guild_tot[which(Soil_Guild_tot$Myc_type == "ECM"),]), method = "glm", method.args = list(family = "binomial"), color = "springgreen3", se = FALSE, size = 6)

	# Orange = North; Green = South; Circle = AM focal tree; Triangle = ECM focal tree.									
	ggsave("C:/Users/eagar/Dropbox/Andrews Project/R Scripts/Example Figures/Path_abundance_soil.png", plot = Abundance, device = "png", dpi = 800)			
}						
