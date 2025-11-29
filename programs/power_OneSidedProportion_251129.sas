* estimate with simulation the probability ('power') that the two-sided 100*(1-alpha_2sided) exact CI;
*  is below a <criterion_at_most> is the true underlying proportion is <p_0>;
* equivalently;
* the prob./power that the proportion is stat.sign. lower than <criterion_at_most>;
*	(exact binomial test with one-sided significance level (1/2)*alpha_2sided); 
**updated version to have nicer tables for output using proc means;

options nodate;

%let n_sim=2000;
%let alpha_2sided=0.10; 
%let confidence_2sided=%sysevalf(100*(1-&alpha_2sided));

* generate &n_sim trials with p_0 as underlying prob and calculate;
* 	the &confidence_2sided exact binomial CI; 
data a; 
call streaminit(37);
do p_0=0.01, 0.02; * the true underlying probability;
do n_subj=40 to 70 by 1; 
do sim=1 to &n_sim;
n_ev=rand('Binomial',p_0,n_subj);
ev=1; n=n_ev;output;
ev=0; n=n_subj - n_ev;output;
end;
end;
end;
run;

proc sort data=a; by p_0 n_subj sim ev ;run;
* calculate a CI for each simulation;
ods exclude all; ods results off;
ods output binomial=proportion;
ods output binomialCLs=cl;
ods output binomialtest=test;
proc freq data=a ;by p_0 n_subj sim; 
weight n/ zeros;* to also count the 0 frequencies;
* p=0.05 is the test null hypothesis, but we 'test' later based on upperbound exact CI;
* level=2 means that sample proportion for the second level of ev (so ev=1) is calculated;
* (via proc sort ev is acending ordered so 0 is the first level and 1 the second level;
table ev /binomial(cl=exact p=0.05 level=2) alpha=&alpha_2sided;* one sided testing;
run;
ods output close;
ods exclude none; ods results on;

ods rtf file="251129_CasperReijnen.doc" style=minimal;
*for several margins calculate wheter CI is below margin or not;
data cl_test; set cl;
do criterion_at_most=0.05,0.10;
if upperCL le criterion_at_most then success=1; else success=0;
output;
end;
run;

title "sample size true proportion p_0, upperbound for the two-sided &confidence_2sided.%-CI at most: ";
proc sort data=cl_test; by p_0 criterion_at_most n_subj;run;
proc means data=cl_test nmiss mean; 
by p_0 criterion_at_most;
class n_subj; 
var success;
run;


ods rtf close;
