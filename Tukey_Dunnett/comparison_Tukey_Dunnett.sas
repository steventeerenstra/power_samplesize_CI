data Tukey;
* to get the critical p-value for the Tukey test;
* use prob=0.90 to get the right-tail probability for 0.10 etc;
* nparms is number of groups , say k;
* df is degrees of freedom of the within-group variance from the anova;
* for equal sized groups of n, this is k*n -k;
* check: results agree with https://www.real-statistics.com/statistics-tables/studentized-range-q-table/;
* below example is from ...;
alpha=0.05
q=.; prob=1-alpha; df=44; nparms=4;
q_crit=probmc('range',q,prob,df,nparms);
run;

proc print;run;

*** the HSD smallest difference that is stat.sign. with Tukey test***;
* is q_crit*sqrt( MS_w/n ) ***************;
** where MS_w is the mean square error within, so the Standard error;
** n= number of subjects in each group;

*** example***;
* https://www.real-statistics.com/one-way-analysis-of-variance-anova/unplanned-comparisons/tukey-hsd/ ;
* can be reproduced with;
data a;
q=.; prob=0.95; df=44; nparms=4;
q_crit=probmc('range',q,prob,df,nparms);
run;
proc print;run;


*************** Dunnett's test *************************;
data a;
* to get the critical p-value for the Dunnett's test;
* use prob=0.95 to get the right-tail probability for 0.05 etc;
* nparms is number of groups without the control group, k, so in total k+1 groups;
* df is degrees of freedom of the within-group variance from the anova;
* for equal sized groups of n, this is k*n -k;
* check: results agree with https://www.real-statistics.com/statistics-tables/studentized-range-q-table/;
* below example is from ...;
alpha=0.05;
q=.; prob=1-alpha; df=28; nparms=3;
t_d=probmc('Dunnett2',q,prob,df,nparms);
q_crit=sqrt(2)*t_d;
run;

proc print;run;

*****************************;
** the critical p-value is now by tradition ***;
*** q_crit=sqrt(2)*t_d **************;
** and the minimal statistical significant difference is ***;
** q_crit*sqrt(MS/n)=t_d*sqrt(2*MS/n)******;


*** example ****;
* https://www.real-statistics.com/one-way-analysis-of-variance-anova/unplanned-comparisons/dunnetts-test-2/;
** t_d=2.483 and q_crit=3.51, df=28 and MS_w=SS/df=5025.624/28; 
data b;
q=.; prob=0.95; df=28; nparms=3;
t_d=probmc('Dunnett2',q,prob,df,nparms);
q_crit=sqrt(2)*t_d;
run;



***** comparison of critical value for Tukey test for alpha=0.025 for 4 groups, so 6 comparisons;
***** and critical value for Dunnett's test for alpha=0.5* 0.025=0.0125 for 2 groups vs control;
* fix k=4 groups of n=10, df=nk-n=36;
data Tukey;
alpha=0.025;
q=.; prob=1-alpha; df=36; nparms=4;
q_crit=probmc('range',q,prob,df,nparms);
run;

proc print;run;

*** for Dunnett: 2 groups vs control, each size n=10, so df=3*10-3=27;
*** note that nparms=number of groups without control;
data Dunnett;
alpha=0.0125;
q=.; prob=1-alpha; df=27; nparms=2;
t_d=probmc('Dunnett2',q,prob,df,nparms);
q_crit=sqrt(2)*t_d;
run;
proc print;run;


data compare;
do n=5 to 18;
* critical value Tukey test for 6 pairwise comparisons for full alpha;
* with 4 groups of size n, so df=4*n-4;
alpha=0.025;
q=.; prob=1-alpha; df=4*n-4; nparms=4;
q_crit_Tukey=probmc('range',q,prob,df,nparms);
*critical value for Dunnett's test for 2 x 2 comparisons vs control;
* each at half alpha;
* three groups (one control), so df=3*n-3;
* nparms= number of groups without control for probmc;
alpha=0.0125;
q=.; prob=1-alpha; df=3*n-3; nparms=2;
t_d=probmc('Dunnett2',q,prob,df,nparms);
q_crit_Dunnett=sqrt(2)*t_d;
**;
output;
end;
run;

proc print data=compare;var n q_crit_Tukey q_crit_Dunnett;run;
