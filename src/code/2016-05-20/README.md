Scripts to fit only correct choice GO trials for Broca and Xena
========================

I've tried fitting correct and error trials multiple ways (see 2016-03-25 and 206-03-31), without much success (none of the models fit error trials well).

I want to determine whether this issue is the error trials. I.e., whether the model will fit the correct trials well if I exclude the error trials.

So I will alter job_preprocess_data.m to exclude error trials from observed data.