Scripts to fit GO trials for  Broca and Xena, using correct vs error within one response side rather than treating left vs right as correct vs error (see below).
========================

Up to this point, I've treated each condition separately, fitting left vs right. Thus, in the easiest cyan condition, correct left responses were fit alongside error right responses. 

This set of models fits correct vs error across matching conditions. Thus, in the easiest cyan and magenta conditions, correct left and error left are fit, as well as correct right and error right. 

To do this, observed data must be re-coded. job_preprocess data needs to match correct and error rather than left vs. right