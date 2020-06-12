%let n_sim=1000;
%let n0=22;
%let n1=22;
%let p0=0.3;
%let p1=0.04;
%let twosided_p=0.10;
data a; 
call streaminit(65537);
do sim=1 to &n_sim;
e0=rand("Binomial", &p0, &n0);e1=rand("Binomial", &p1, &n1);
*ctl;
trt=0; suc6=0; count=&n0-e0; output;
trt=0; suc6=1; count=e0; output;
*trt;
trt=1; suc6=0; count=&n1-e1; output;
trt=1; suc6=1; count=e1; output;
end;
run;

ods select none;
ods output ChiSq=Chisq;
ods output FishersExact=FishersExact;
proc freq data=a; by sim;
table trt*suc6/chisq;
weight count;
run;
ods output close;
ods select all;

data ChiSq_rej; set ChiSq; if Statistic="Chi-Square";reject=(prob< &twosided_p);run;
* use table probability for FisherExact;
data FishersExact_rej; set FishersExact; if Name1="P_TABLE"; reject=(nvalue1< &twosided_p);run;
proc means data=ChiSq_rej n mean; var reject;run;

proc means data=FishersExact_rej n mean; var reject;run;
