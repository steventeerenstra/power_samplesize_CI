%let n_subj=95;
%let n_sim=2000;
%let p_0=0.005;

* here an example for testing one-sided: H0: p>0.05 vs H1: p <= 0.05;

data a; 
do sim=1 to &n_sim;
call streaminit(37);
n_subj=&n_subj;
n_ev=rand('Binomial',&p_0,&n_subj);
ev=1; n=n_ev;output;
ev=0; n=n_subj - n_ev;output;
end;
run;

proc sort data=a; by sim ev ;run;
ods exclude all; ods results off;
ods output binomial=proportion;
ods output binomialCLs=cl;
ods output binomialtest=test;
proc freq data=a ;by sim; 
weight n/ zeros;* to also count the 0 frequencies;
* p=0.05 is the test hypothesis;
table ev /binomial(cl=exact p=0.05 level=2) alpha=0.10 ;* one sided testing;
run;
ods output close;
ods exclude none; ods results on;
data cl_test; set cl;
if  upperCL le 0.05 then success=1; else success=0;
run;

title "sample size &n_subj, proportion &p_0";
proc freq data=cl_test;
table success/missing;run;
