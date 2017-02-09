Scripts to fit altered data for Broca and Xena
========================

Previous experiments to this point had excluded a bunch of Stop trials because of the way the data was translated.
I tried to fix it, but found out I was still excluding lots of Stop trials. 
I believe I've corrected the bug (in job_preprocess_data.m, I fixed the iAbort variable to be specific for certain types of aborts).
I will run the experiment again - see 2017-01-02.