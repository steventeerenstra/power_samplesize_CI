** precision in HR given a number of events; 
data a; 
alpha=0.10;*two-sided (1-alpha)*100% confidence interval;
do hr=0.9, 0.95, 1.0, 1.05, 1.10;
loghr=log(hr);
do fraction= 0.2, 0.3, 0.4, 0.5, 0.6,0.8, 1.0;
total_samplesize=312+315;
events=fraction*total_samplesize;
logse=	1/sqrt(events*0.5*(1-0.5));
ul_loghr=loghr + probit(1-alpha/2)*logse;
ul_hr= exp(ul_loghr);
ll_loghr=loghr-probit(1-alpha/2)*logse;
ll_hr=exp(ll_loghr);
output;
end; end;
run;
proc print data=a noobs;by hr;var total_samplesize fraction ll_hr ul_hr;run;

*** critical value in terms of HR given a number of events and critical p-value (at interim);
data a; 
alpha=0.0001; *two-sided;
hr=0.7;
theta=0.5; * proportion in experimental group;
events=83;* number of events;
logse=1/sqrt(events*theta*(1-theta)); * standard error on logscale;
loghr=-probit(1-alpha/2)*logse; * to get loghr + probit(1-alpha/2)*logse < log(1);
hr=exp(loghr);
run;
proc print;run;






** expected number of events for different follow-up, so projection of timings of events *****;

data a;
* according to Cox proportional hazards, Fundamentals of clinical trials, Friedman et al., p. 154;
accrual= 26 ; *accrual time;
do fu=1, 10, 26; * followup after finishing recruitment;
median_c= 28.5; * median control group;
median_i=median_c/0.7; * median intervention group;
N_c= 220; * sample size in the control group;
N_i= 220; * sample size intervention group;
* control group;
lambda_c= - log(0.5)/median_c;
lambda_i= -log(0.5)/median_i;
hr=lambda_i/lambda_c;
proportiondeaths_c=(	exp(-lambda_c*fu) -  exp(-lambda_c*(accrual+fu))	) / (lambda_c*accrual);
events_c=N_c*(1- proportiondeaths_c);
* intervention group;
lambda_i= - log(0.5)/median_i;
proportiondeaths_i=(	exp(-lambda_i*fu) -  exp(-lambda_i*(accrual+fu))	) / (lambda_i*accrual);
events_i=N_i*(1- proportiondeaths_i);
* total numbers of events;
events=events_c+events_i;
output;
end;
run;

proc print data=a;run;
***********************************


* valided for OS with 15 mo vs 20 mo, 
* 32 pat/month, and 664 subjects in total, 332 vs 332,
* gives 418 deaths instead of 413 deaths in 36 total study duration, 
* so approx ok;
******************************************************************************;

** PFS of renal trial **;
data a;
* according to Cox proportional hazards, Fundamentals of clinical trials, Friedman et al., p. 154;
accrual= 664/32; *accrual time: within 17 months < 664, i.e. total sample size is expected, so no followup yet;
fu=23-664/32; * followup after finishing recruitment;
median_c= 15; * median control group;
median_i= 20; * median intervention group;
N_c= 332; * sample size in the control group;
N_i=332; * sample size intervention group;
* control group;
lambda_c= - log(0.5)/median_c;
lambda_i= -log(0.5)/median_i;
hr=lambda_i/lambda_c;
proportiondeaths_c=(	exp(-lambda_c*fu) -  exp(-lambda_c*(accrual+fu))	) / (lambda_c*accrual);
events_c=N_c*(1- proportiondeaths_c);
* intervention group;
lambda_i= - log(0.5)/median_i;
proportiondeaths_i=(	exp(-lambda_i*fu) -  exp(-lambda_i*(accrual+fu))	) / (lambda_i*accrual);
events_i=N_i*(1- proportiondeaths_i);
* total numbers of events;
events=events_c+events_i;
run;

proc print data=a;run;


*********************************;
** PFS of liver trial **;
data a;
* according to Cox proportional hazards, Fundamentals of clinical trials, Friedman et al., p. 154;
accrual= 30; *total accrual time, so 576/30=19.2 per month;
fu=0.5; * followup after finishing recruitment;
median_c= 8.2; * median control group;
median_i= 11.39; * median intervention group;
N_c= 192; * sample size in the control group;
N_i=2*192; * sample size intervention group;
* control group;
lambda_c= - log(0.5)/median_c;
lambda_i= -log(0.5)/median_i;
hr=lambda_i/lambda_c;
proportiondeaths_c=(	exp(-lambda_c*fu) -  exp(-lambda_c*(accrual+fu))	) / (lambda_c*accrual);
events_c=N_c*(1- proportiondeaths_c);
* intervention group;
lambda_i= - log(0.5)/median_i;
proportiondeaths_i=(	exp(-lambda_i*fu) -  exp(-lambda_i*(accrual+fu))	) / (lambda_i*accrual);
events_i=N_i*(1- proportiondeaths_i);
* total numbers of events;
events=events_c+events_i;
run;

proc print data=a;run;
 

* intermediate analysis within the accrual time;
data a;
* according to Cox proportional hazards, Fundamentals of clinical trials, Friedman et al., p. 154;
accrual= 23.5 ;
fu=0; * followup after finishing recruitment;
median_c= 8.2; * median control group;
median_i= 11.39; * median intervention group;
N_c= accrual*19.2*(1/3); * sample size in the control group;
N_i=accrual*19.2*(2/3); * sample size intervention group;
* control group;
lambda_c= - log(0.5)/median_c;
lambda_i= -log(0.5)/median_i;
hr=lambda_i/lambda_c;
proportiondeaths_c=(	exp(-lambda_c*fu) -  exp(-lambda_c*(accrual+fu))	) / (lambda_c*accrual);
events_c=N_c*(1- proportiondeaths_c);
* intervention group;
lambda_i= - log(0.5)/median_i;
proportiondeaths_i=(	exp(-lambda_i*fu) -  exp(-lambda_i*(accrual+fu))	) / (lambda_i*accrual);
events_i=N_i*(1- proportiondeaths_i);
* total numbers of events;
events=events_c+events_i;
run;

proc print data=a;run;
 

** power from events, here 850 events gives 84% power for two-sided alpha is 0.025**;
data a;
alpha=0.05; * two-sided alpha;
d=(1-0.17)*850;
theta=0.5; * fraction randomized to experimental, = 0.5 for 1:1 rando;
se=1/sqrt( d*theta*(1-theta) );
hr=0.8;
delta=abs( log(hr) );
power=probnorm(delta/se - probit(1-alpha/2) );*here the division by 2 for one-sided alpha;
run;

proc print;run;

/** formulas
d= number of events;
delta=abs( log(hr) );
se= 1/sqrt(theta*(1-theta)*d)
power=probnorm(delta/se - probit(1-alpha/2) );
d=(z_(1-alpha/2) + z_(1-beta))**2  / (delta**2 * theta * (1-theta)) 
**/


*** number of events needed ***;
data a;
alpha=0.05; beta=0.05;
theta=2/3;
hr=0.5;
delta=abs(log(hr) );
d=( probit(1-alpha/2) + probit(1-beta) )**2  / ( delta**2 * theta * (1-theta)) ;
run;

proc print;run;


***** Point Estimate of hr < 1 with a certain power ****;
data a;
* if the true effect is on log scale log_hr***;
* then it its sampling distribution for a given standard error log_se ****;
* is N(log_hr, log_se), so the point estimate is < 0 with probability power if ****;
* log_hr+ z_{power}*log_se <=0 , because z_{power} has as left (so cumulative) probability "power";  
d= 211; *number of events;
theta=1/3; *fraction in the control arm;
power=0.8;
log_se=1/sqrt( theta*(1-theta)*d  );
log_hr =-probit(power)*log_se;
hr=exp(log_hr);
run;
proc print data=a;run; 


**** power and precision (95-CI) from accrual, follow-up, number of patients etc ***;
** from different scenarios put in a dataset;
data b; input hr median_c N_c N_i accrual fu;
datalines;
0.70  3 100 100 12 6
0.65  3 100 100 12 6
;
run;


data a;set b;
do alpha=0.05, 0.10; * two-sided;
power_aim=0.8;
* according to Cox proportional hazards, Fundamentals of clinical trials, Friedman et al., p. 154;
* provide hr N_c N_i accrual fu (=follow-up after last patient;
theta=N_c/(N_i+N_c);
* events needed;
d_needed=(probit(1-alpha/2) + probit(power_aim))**2 /   ( (log(hr))**2 * theta*(1-theta));
* control group;
lambda_c= - log(0.5)/median_c;
proportiondeaths_c=(	exp(-lambda_c*fu) -  exp(-lambda_c*(accrual+fu))	) / (lambda_c*accrual);
events_c=N_c*(1- proportiondeaths_c);
* intervention group;
lambda_i= hr*lambda_c;
proportiondeaths_i=(	exp(-lambda_i*fu) -  exp(-lambda_i*(accrual+fu))	) / (lambda_i*accrual);
events_i=N_i*(1- proportiondeaths_i);
* total numbers of events;
events=events_c+events_i;
* power;
se=1/sqrt( events*theta*(1-theta) );
delta=abs( log(hr) );
power=probnorm(delta/se - probit(1-alpha/2) );
* precision around the target hr;
hr_95lower=exp( log(hr)- 1.96*se    ); 
hr_95upper=exp( log(hr)+ 1.96*se    ); 
output;
end;
run;

proc print data=a noobs;var alpha power_aim d_needed accrual fu median_c hr N_c N_i theta events power hr_95lower hr hr_95upper; 
run;   


**** power and precision (95-CI) from accrual, follow-up, number of patients etc ***;
** from different scenarios put in a dataset in terms of control proportion that is eventfree at certain milestone;
data b; input hr eventfree_c milestone_c N_c N_i accrual fu;
datalines;
0.6  0.6 60 106 106 36 60
;
run;


data a;set b;
do alpha=0.05, 0.10; * two-sided;
power_aim=0.8;
* according to Cox proportional hazards, Fundamentals of clinical trials, Friedman et al., p. 154;
* provide hr N_c N_i accrual fu (=follow-up after last patient;
theta=N_c/(N_i+N_c);
* events needed;
d_needed=(probit(1-alpha/2) + probit(power_aim))**2 /   ( (log(hr))**2 * theta*(1-theta));
* control group;
lambda_c= - log(eventfree_c)/milestone_c;
proportiondeaths_c=(	exp(-lambda_c*fu) -  exp(-lambda_c*(accrual+fu))	) / (lambda_c*accrual);
events_c=N_c*(1- proportiondeaths_c);
* intervention group;
lambda_i= hr*lambda_c;
proportiondeaths_i=(	exp(-lambda_i*fu) -  exp(-lambda_i*(accrual+fu))	) / (lambda_i*accrual);
events_i=N_i*(1- proportiondeaths_i);
* total numbers of events;
events=events_c+events_i;
* power;
se=1/sqrt( events*theta*(1-theta) );
delta=abs( log(hr) );
power=probnorm(delta/se - probit(1-alpha/2) );
* precision around the target hr;
hr_95lower=exp( log(hr)- 1.96*se    ); 
hr_95upper=exp( log(hr)+ 1.96*se    ); 
output;
end;
run;

proc print data=a noobs;var alpha power_aim d_needed accrual fu eventfree_c milestone_c hr N_c N_i theta events power hr_95lower hr hr_95upper; 
run;
