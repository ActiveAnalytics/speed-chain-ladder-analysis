speed-chain-ladder-analysis
===========================

Speed chain ladder analysis with Rcpp


This repository contains the C++ code for a basic chain ladder algorithm. It is designed to complete the run off triangle only and not provide a full posthoc analysis of standard errors etc. It also does not account for missing data in the run-of triangle. It was created for a blog entry at Active Analytics:

http://www.active-analytics.com/blog/speedchainladderanalysiswithrcpp/

The blog entry contains more details of usage and implementation.

# The files

This repository contains two files:

1. chain-ladder.cpp: contains the C++ code for simulating and analysing the chain ladder
2. chain-ladder.r: contains the R code that calls the C++ code

See the R code and the blog for more details.
