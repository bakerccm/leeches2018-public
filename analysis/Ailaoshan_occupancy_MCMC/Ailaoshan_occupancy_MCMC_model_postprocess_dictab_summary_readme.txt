# from the dictab help page:

‘dictab’ creates an object of class ‘dictab’ with the following
components:

Modname: the name of each model of the candidate model set.

 pD: the effective number of estimated parameters for each model.

DIC: the deviance information criterion for each model.

Delta_DIC: the delta DIC of each model, measuring the difference in DIC
     between each model and the top-ranked model.

ModelLik: the relative likelihood of the model given the data
     (exp(-0.5*delta[i])).  This is not to be confused with the
     likelihood of the parameters given the data.  The relative
     likelihood can then be normalized across all models to get
     the model probabilities.

DICWt: the DIC weights, sensu Burnham and Anderson (2002) and
     Anderson (2008). These measures indicate the level of support
     (i.e., weight of evidence) in favor of any given model being
     the most parsimonious among the candidate model set.

Cum.Wt: the cumulative DIC weights.  These are only meaningful if
     results in table are sorted in decreasing order of DIC
     weights (i.e., ‘sort = TRUE’).

Deviance: the deviance of each model.
