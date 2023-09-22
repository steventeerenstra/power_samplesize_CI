

%macro power_its_cont(alpha=0.05,ds_measurements=, ds_out=out);
* for continuous outcome;
/*ds_measurements is a dataset with variables: config, cluster , k, m, delta, sd, n, rho, r, r_dmax, icc
which holds for each cluster (cluster =1,2,..)
        the number of measurements before (k)and after (m) the intervention
        the ICC (rho) of the cluster, cluster size (n),
        the correlation between adjacent measurement (r)
        the correlation between the first and last measurement (r_dmax),
        the distance between the first and last measurement is set to k+m-1
also the configuration (config) can be named, if there are multiple;
*/
* get number of clusters;
data _null_;
  set &ds_measurements NOBS=size;
  call symput("number_clusters",size);
 stop;
 run;

* dataset to calculate the variances of the ITS of each cluster;
data a;
set &ds_measurements;
alpha=&alpha;
sigma2_tot=sd**2;*variance of the cluster;
de=1+(n-1)*rho;*design effect for the cluster;
var_clustermean=(de/n)*sigma2_tot;*variance average in the cluster;
*variance of one ITS (under CS/AR1 correlation);
var_its_CS=var_clustermean*(1-r)*(1/k+1/m);
var_its_AR1=var_clustermean*( (1+r)/(1-r)*(1/k+1/m) - (2*r)/(1-r)**2 *( ((1-r**k)/k)*((1-r**m)/m) + (1-r**k)/k**2 + (1-r**m)/m**2 )  );
 * rr is the correlation in the LEAR model such that ;
 * the correlation is a*rr^k if the distance is k;
 * and if k=1 then a*rr=r and if k=d_max then a*rr^(d_max)=r_dmax;
 * so the correlation decreases from r to r_dmax;
 dmax=k+m-1;*assuming k meas. directly before and m meas. directly after switch;
 rr=(r_dmax/r)**(1/(dmax-1)) ;
 a=r/rr;
 var_its_LEAR=var_clustermean*( (1+(2*a-1)*rr)/(1-rr)*(1/k+1/m) - (2*a*rr)/(1-rr)**2 *( ((1-rr**k)/k)*((1-rr**m)/m) + (1-rr**k)/k**2 + (1-rr**m)/m**2 )  );
%END;
output;
run;

* for each configuration, collect the clusters corresponding to it;
proc sort data=a; by config cluster;run;

* for each configuration, calculate the total variance and average variance, se and power;
data b; set a; by config;
retain ivar_its_CS_sum ivar_its_AR1_sum ivar_its_LEAR_sum
       est_CS_partialsum est_AR1_partialsum est_LEAR_partialsum 0;
* for each configuration;
    *reset the sums and calculate ...;
if first.config then do;
    ivar_its_CS_sum=0; ivar_its_AR1_sum=0;ivar_its_LEAR_sum=0;
    est_CS_partialsum=0; est_AR1_partialsum=0;est_LEAR_partialsum=0;
    end;
    * sum of inverse variances;
ivar_its_CS_sum=ivar_its_CS_sum     +1/var_its_CS;
ivar_its_AR1_sum=ivar_its_AR1_sum   +1/var_its_AR1;
ivar_its_LEAR_sum=ivar_its_LEAR_sum +1/var_its_LEAR;
    * partial sum of weighted estimate (up to the current cluster);
est_CS_partialsum=  est_CS_partialsum   + delta*1/var_its_CS;
est_AR1_partialsum= est_AR1_partialsum  + delta*1/var_its_AR1;
est_LEAR_partialsum=est_LEAR_partialsum +delta*1/var_its_LEAR;
* variance of meta-analysis estimator;
* (i.e. of average over all interrupted time series up to current cluster);
var_meta_its_CS= 1/(ivar_its_CS_sum); *inverse of inverse variances;
var_meta_its_AR1=1/(ivar_its_AR1_sum);
var_meta_its_LEAR=1/(ivar_its_LEAR_sum);
* se of meta-analyis estimator ((for the ITS up to the current cluster);;
se_meta_its_CS=sqrt(var_meta_its_CS);
se_meta_its_AR1=sqrt(var_meta_its_AR1);
se_meta_its_LEAR=sqrt(var_meta_its_LEAR);
* meta-analysis average up to current cluster;
est_meta_its_CS=est_CS_partialsum/(ivar_its_CS_sum);
est_meta_its_AR1=est_AR1_partialsum/(ivar_its_AR1_sum);
est_meta_its_LEAR=est_LEAR_partialsum/(ivar_its_LEAR_sum);
* power according to normal approximation;
power_CS=probnorm(est_meta_its_CS/se_meta_its_CS - probit(1-alpha/2) );
power_AR1=probnorm(est_meta_its_AR1/se_meta_its_AR1 - probit(1-alpha/2) );
power_LEAR=probnorm(est_meta_its_LEAR/se_meta_its_LEAR - probit(1-alpha/2) );
if cluster=1 then do; *the df for the t-distribution is not defined;
 power_CS_df=.;
 power_AR1_df=.;
 power_LEAR_df=.;
 power_inbetween_df=.;
end;
else if cluster > 1 then do;
 power_CS_df=probt(est_meta_its_CS/se_meta_its_CS - tinv(1-alpha/2,cluster-1), cluster -1);
 power_AR1_df=probt(est_meta_its_AR1/se_meta_its_AR1 - tinv(1-alpha/2,cluster-1), cluster -1);
 power_LEAR_df=probt(est_meta_its_LEAR/se_meta_its_LEAR - tinv(1-alpha/2,cluster-1), cluster -1);
end;

run;

title6 "configurations:";
proc print data=&ds_measurements noobs;by config;run;

* only the powers etc from the last cluster give the powers of the whole design;
data &ds_out; set b; by config;
if last.config;
n_cluster=cluster;
run;

title5 "n cluster size, rho ICC in each measurement";
title6 "r=corr. between subsequent measurements; r_dmax= corr. between first and last measurement)";
* se_meta_its_CS power_CS se_meta_its_AR1 power_AR1 se_meta_its_inbetween power_inbetween;
proc print data=&ds_out noobs; by config;
var n_cluster r dmax r_dmax n rho power_LEAR_df power_CS_df power_AR1_df;*  power_CS_df power_AR1_df power_inbetween_df;
; run;
%mend;



/* test case for comparison with old program
title1 "RENEW: 36 clusters, alpha=0.05 (two-sided)";
data init2;
config=1;
do cluster=1 to 36;
k=2; m=2; delta=0.3*0.17; sd=0.5; n=70;rho=0.05; r=0.8; r_dmax=0.5;
output;
end;
run;

title2 "correlation over time (from 0.8 [adjecent periods] to 0.5 [first and last period])";
%power_its_cont(alpha=0.05, ds_measurements=init2);
*/

ods rtf file="230922_RENEW_BenjaminWendt.rtf";
title1 "RENEW: 36 clusters, alpha=0.05 (two-sided), icc=0.05";
title2 "50% of clusters around 35%, 25% of clusters around 10%, 25% of clusters around 50%";
title3 "high correlation between periods";
title4 "relative reduction 17%";
data init2;
config=1;
r=0.8; r_dmax=0.5;k=2;m=2;rho=0.08;reduction=0.17;
* clusters around 35%;
do cluster=1 to 18;
n=70; p_ctl=0.35; delta=p_ctl*reduction; sd=sqrt(p_ctl*(1-p_ctl));
output;
end;
do cluster=19 to 27;
n=80; p_ctl=0.5; delta=p_ctl*reduction; sd=sqrt(p_ctl*(1-p_ctl));
output;
end;
do cluster=28 to 36;
n= 50;p_ctl=0.1; delta=p_ctl*reduction; sd=sqrt(p_ctl*(1-p_ctl));
output;
end;
run;

title2 "correlation over time (from 0.8 [adjacent periods] to 0.5 [between first and last period])";
%power_its_cont(alpha=0.05, ds_measurements=init2);

****  smaller autocorrelation ****; 
title1 "RENEW: 36 clusters, alpha=0.05 (two-sided), icc=0.05";
title2 "50% of clusters around 35%, 25% of clusters around 10%, 25% of clusters around 50%";
title3 "moderate correlation between periods";
title4 "relative reduction 17%";
data init2;
config=1;
r=0.6; r_dmax=0.5;k=2;m=2;rho=0.08;reduction=0.17;
* clusters around 35%;
do cluster=1 to 18;
n=70; p_ctl=0.35; delta=p_ctl*reduction; sd=sqrt(p_ctl*(1-p_ctl));
output;
end;
do cluster=19 to 27;
n=80; p_ctl=0.5; delta=p_ctl*reduction; sd=sqrt(p_ctl*(1-p_ctl));
output;
end;
do cluster=28 to 36;
n= 50;p_ctl=0.1; delta=p_ctl*reduction; sd=sqrt(p_ctl*(1-p_ctl));
output;
end;
run;

title2 "correlation over time (from 0.6 [adjacent periods] to 0.5 [between first and last period])";
%power_its_cont(alpha=0.05, ds_measurements=init2);
ods rtf close;



