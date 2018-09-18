Supplementary Materials for "Alternative Models for Moment Inequalities"

The material here is used to both generate simulated market data described in 3.2 and provide estimates presented in Table 2.  The programs are written in C and Matlab.

Questions can be sent to Robin Lee (rslee@stern.nyu.edu).

Files:

SimulateMarkets.cpp : C program used to generate random markets and output primitives and equilibrium networks/transfers.

SimulateMarkets.pdf : Writeup describing the underlying contracting game and model which is simulated.

Main_ImportAll.m : Imports the data generated from SimulateMarkets.cpp, places them into a Matlab file structure, and estimates Theta0 and the identified set for the entire set of markets.

Main_GenerateDraws.m : Takes imported data and splits it up into sets of 1385 markets (4 observations per market).  

Main_Inequalities_Part1.m and Main_Inequalities_Part2.m :
For each set of 1385 markets, runs the appropriate inequality estimator with different error specifications, and outputs estimates for Table 2 (lines 1-11, 13-14).  Structured to be run as multiple instances on a Sun Grid Engine.

f_*.m :
Supporting functions used in Main_Inequalities_Part1.m and Main_Inequalities_Part2.

