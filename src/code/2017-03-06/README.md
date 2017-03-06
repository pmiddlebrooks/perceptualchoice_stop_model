Scripts to fit Broca and Xena, corrected a starting value issue in a notebook.
========================

Previous experiments, in the final fitting step using at 1-ms resolution, the code used best fits from 5-ms resolution that had fit the STOP trials only (with GO parameters fixed).
This could have caused the odd fitting result (complex models fitting worse).

In notebook Perceptual Choice Stop Model, last section, fixed the code from querying "allFValAndBestX_stopTrials_model" to "allFValAndBestX_allTrials_model".

Also increased number of simulated trials, to make fitting more accurate in general.