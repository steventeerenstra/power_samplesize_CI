/** ods output, title 
ods rtf file="????.rtf"; 
title "???";
**/

* the median, Q1, Q3 are provided on the original scale;
* in case of multiple clusters (so estimation of ICC as well):
*   the dataset contain one record for each cluster;

data a; input raw_med raw_q1 raw_q3;
datalines;
3.0 2.0  6.0
4.0 2.0  7.0
5.0 3.0  9.0
4.0 2.0  6.0
4.0 2.0  7.0
4.0 2.0  8.0
5.0 3.0  9.0
3.0 2.0  5.0
5.0 2.0  9.0
4.0 2.0  7.0
4.0 2.0  7.0
;
run;


* calculating mean and sd on log-scale;
data a1; set a;
* on logscale;
mean=log(raw_med); q1=log(raw_q1); q3=log(raw_q3);
midpoint_q1_q3=(q1+q3)/2; * to check whether the data is lognormally distributed; 
* distance q3 to q1;
dist=probit(0.75)-probit(0.25);
sd=(q3-q1)/(dist);
sd_h=(q3-mean)/(dist/2);* to check whether left and right of the mean the same distance;
sd_l=(mean-q1)/(dist/2); *idem;
sd2=sd**2;
run;
title2 "on natural log scale"; 
proc print data=a1 noobs; var raw_med raw_q1 raw_q3 mean midpoint_q1_q3 q1 q3 sd sd2;
run;
proc means data=a1 mean var; var mean sd2; 
output out=a2 var(mean)=between_sd2 mean(sd2)=within_sd2 ;
run;

* in case only one record, set the between-variance to zero;
data a2; set a2; if between_sd2=. then between_sd2=0;run;

data icc; set a2(drop=_type_ _freq_);
icc=between_sd2/(between_sd2 + within_sd2);
sd2_total=between_sd2 + within_sd2; sd_total=sqrt(sd2_total);
run;
proc odstext data=icc;
p "The between-variance is "||put(between_sd2,10.8)||
	" based on the variance of the medians."/style=[fontsize=12pt];
p "The within-variance is " ||put(within_sd2,10.8)||
	" calculated as the average over all the clusters of their within-variances."/style=[fontsize=12pt]; 
p "The ICC is "||put(icc,10.8)||"."/style=[fontsize=12pt];
run;

** example of cluster RCT, parallel cluster or standard Stepped wedge **;
title3 "SW design, without accounting for transition periods";
title4 " I=number of clusters, T=number of timepoints, N=cluster size";
title5 "power_sw_df=power stepped wedge, power_par=power parallel cluster trial";
data a;set icc;
do I=9,10,11,12; 
do T=I+1;
do N=6,8,10; *number of subjects per cluster;
do	rho = icc;
do theta=log(0.75); *effect= difference in means;
do sigma_tot=sd_total; * total sd; 
	sigma2=(1-rho)*sigma_tot**2/N; * variance of a cluster mean due to sampling N subjects;
	tau2= rho* sigma_tot**2; *between cluster variance;
	U=(I*T)/2;W= (I**2 * T *(2*T-1) ) / (6*(T-1)); V= ( I*T*(2*T-1) ) / 6;
	var_theta= ( I*sigma2*(sigma2 + T*tau2) ) 
                    /
                    ( (I*U-W)*sigma2 + (U**2 +I*T*U - T*W - I*V)*tau2  );

	* cluster stepped wedge;
	power_sw=probnorm(-1.96 + abs(theta) / sqrt(var_theta) );
	*based on I = number of clusters degrees of freedom;
	power_sw_df=probt( tinv(0.025,I) + abs(theta) / sqrt(var_theta), I );
	* parallel cluster;
	es=abs(theta)/sigma_tot; *Cohen's effect size;
	de= (1+ (N-1)*rho);; * the design factor;
	power_par=probnorm(-1.96 + es/sqrt(2) * sqrt( (I*N)/de) ); 
	*based on I=number of clusters as degree of freedom;
	power_par_df=probnorm( tinv(0.025,I-2) + es/sqrt(2) * sqrt( (I*N)/de) ); 
	output;
end;end;end;end;end;end;
run;

* for traffic lighting; 
proc format;
value powerfmt
low -<0.70 = 'lightgray'
0.70-<0.79 = 'green'
0.79-high = 'lightgreen'
;
run;

* rename T and N because already in use as keywords in proc tabulate";
proc tabulate data=a(rename=(I=n_clus T=n_period N=n_clusperiod)); by theta sigma_tot rho;
class n_clus n_period n_clusperiod;
var power_sw_df; *also power for parallel-group cluster RCT: power_sw power_par_df  power_par;
table n_clus*n_period, n_clusperiod*power_sw_df=" "*mean=" "*[style=[background=powerfmt.]]*f=5.4;
 ;run;
ods startpage=no;
proc odstext;
p "Note the above does not account for power loss due to transition periods"/ style=[just=c fontsize=12pt];
run;
ods rtf close;
