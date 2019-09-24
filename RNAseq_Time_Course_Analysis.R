################### Limma ####################

### A small number of time points
lev <- c("wt.0hr","wt.6hr","wt.24hr","mu.0hr","mu.6hr","mu.24hr")
f <- factor(targets$Target, levels=lev)
design <- model.matrix(~0+f)
colnames(design) <- lev
fit <- lmFit(eset, design)

# Which genes respond at either the 6 hour or 24 hour times in the wild-type?
cont.wt <- makeContrasts("wt.6hr-wt.0hr", "wt.24hr-wt.6hr", levels=design)

fit2 <- contrasts.fit(fit, cont.wt)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="BH")

#Any two contrasts between the three times would give the same result. The same gene list would be obtained had "wt.24hr-wt.0hr" been used in place of "wt.24hr-wt.6hr"

# Which genes respond differently over time in the mutant relative to the wild-type?
cont.dif <- makeContrasts(Dif6hr=(mu.6hr-mu.0hr)-(wt.6hr-wt.0hr), Dif24hr=(mu.24hr-mu.6hr)-(wt.24hr-wt.6hr), levels=design)

fit2 <- contrasts.fit(fit, cont.dif)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="BH")


### Many time points
# It might be reasonable to represent a time course for a particular gene in a particular condition using a cubic spline curve with a modest number of knots. Choosing effective degrees of freedom to be in range 3~5 is reasonable. Setup a basis for a natural regression spline:
library(splines)
X <- ns(targets$Time, df=5)
# Then fit separate curves for the control and treatment groups:
Group <- factor(targets$Group)
design <- model.matrix(~Group*X)
fit <- lmFit(y, design)
fit <- eBayes(fit)

# This creates a model with 12 parameters, with the last 5 corresponding to interaction, i.e., to differences in the curves between groups. To detect genes with different time trends for treatment vs control:
> topTable(fit, coef=8:12)
# This conducts a moderated F-test for each gene on 5 df, which can detect very general differences between the treatment and control curves.
# Note that for this analysis, it is not necessary to have replicates, nor is it necessary for the two treatment groups to be observed at identical time points.
