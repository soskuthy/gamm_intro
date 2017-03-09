# gamm_intro

This is a practical introduction to GAMMs in the context of linguistic research. This GitHub folder contains the example data sets and a few additional files.

This is an updated version of the introduction. Here is an (incomplete) list of changes. 

*Important changes*

- there was an error in the way ordered factors were set up for difference smooths; this has been fixed
- added discussion of residual autocorrelation and autoregressive models
- the previous version implied that comparison of models with different fixed effects using REML/fREML is possible but a bit dodgy; it's actually very dodgy and shouldn't be done at all, so I've changed the text
- changed words.50 data set so that the "traj" variable is automatically treated as a factor (some people had issues fitting the models in this part of the introduction since mgcv tried to treat "traj" as a numeric variable)

*Smaller changes*

- fixed issues with plots showing basis functions postmultiplied by model coefficients
- more detailed discussion of selection of number of basis functions and estimation methods
- changed words.50 data set so that the "traj" variable is automatically treated as a factor
- removed analysis of words.50 data set using gamm4() (bam() does a better job)
- moved discussion of different GAM(M) functions to appendix
- the function bam() is used consistently throughout the entire introduction in the latest version
