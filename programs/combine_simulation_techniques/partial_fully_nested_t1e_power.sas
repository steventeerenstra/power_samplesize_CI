/** estimate power and type I error control for partially / fully nested trial 
    first macros are defined
    - a macro is defined to generate data, 
    - a macro to analyze data (this stores per configuration fixed effects, random effects and convergence as datasets), 
	- an oversight macro that takes in the configurations from a configurations dataset and calls the above two

 ** the results are post-processed to get the estimated type I errors and power
 ** analysis time of this program is faster if batch processed and results are in the .lst file
**/

%macro sim_outcomes(ds_sim=ds_sim,seed=37,
n_sim=, n_clus0=,n_clus1=, n_subj0=, n_subj1=, 
 s2_c=, s2_ct0=, s2_s=, s2_st0=, s2_0ct1=, s2_0st1=, s2_1ct1=,s2_1st1=,
 fe_t0=, fe_t1=, delta=, 
name_config=notspecified
);
*notation: fe=fixed effect (for time=0,1), s2=variances of random effect/components ;
data &ds_sim; length name_config $ 100;
name_config="&name_config"; * unique name to identify the configuration;
k0=&n_clus0; k1=&n_clus1; n0=&n_subj0; n1=&n_subj1;
array n_clus {0:1} k0 k1;
array n_subj {0:1} n0 n1;
* the variance components from the correlations to total SDs; 
s2_c= &s2_c; s2_ct0=&s2_ct0;
s2_s= &s2_s;  s2_st0=&s2_st0;
s2_0ct1= &s2_0ct1; s2_0st1= &s2_0st1;
s2_1ct1= &s2_1ct1;s2_1st1= &s2_1st1;
* set seed to reproduce the simulation;
call streaminit(&seed);
* subjects are (clus,subj)-indexed according to their allocation at follow-up;
* e.g. if no clustering at baseline then s2_c=s2_ct0=0 so subjects are in clusters with ICC=0;
* clusters can be of size 1 if they are individuals;
do sim=1 to &n_sim;
 do tx=0 to 1; * treatment is constant over clusters (Situation 1-3);
   * random intercepts for each cluster, and random cluster x time interaction for t0 and t1;
   do clus=1 to n_clus{tx};rclus=rand("Normal"); rclus_t0=rand("Normal"); rclus_t1=rand("Normal"); 
     * random intercept for each subject, random subject x time interaction for t0 and t1;
     do subj=1 to n_subj{tx};rsubj=rand("Normal"); rsubj_t0=rand("Normal");rsubj_t1=rand("Normal");
	   * generate longitudinal data;
  	   do time=0 to 1; 
	     do hypo=0 to 1;
		 * only a treatment effect for tx=1 at time=1 if hypo=1;
		 fixed_effects=(time=0)*&fe_t0 + (time=1)*&fe_t1 + (time=1)*(hypo=1)*(tx=1)*&delta;
         * a cluster (hence its subjects) has either tx=0 or tx=1, so below: (tx=0)*.. =0 or (tx=1)*.. =0;
		 * therefore we can suffice with one rclus_t1 rsubj_t1 regardless of tx=0,1;
		 random_effects= sqrt(s2_c)*rclus + sqrt(s2_s)*rsubj 
                       + (time=0)*(  sqrt(s2_ct0)*rclus_t0 + sqrt(s2_st0)*rsubj_t0  ) 
                       + (time=1)*(	
									(tx=0)*(  sqrt(s2_0ct1)*rclus_t1 + sqrt(s2_0st1)*rsubj_t1  )
									  +
                       				(tx=1)*(  sqrt(s2_1ct1)*rclus_t1 + sqrt(s2_1st1)*rsubj_t1  )
								   );
		outcome=fixed_effects + random_effects;
		* dummy variables useful for the analyses;
		t1=(time=1); t1_g0=(time=1)*(tx=0); t1_g0_class=t1_g0;t1_g1=(time=1)*(tx=1);t1_g1_class=t1_g1;
		ctl_clus= (tx=1)*clus; * so =0 for all in ctl arm, =j for cluster j in intervention arm;
        output;
end;end;end;end;end;end;
run; 
%mend sim_outcomes;


%macro analyze_outcomes(ds_sim=ds_sim, name_config=name_config, alpha=0.05,
						ds_conv=ds_conv,ds_feparms=ds_feparms, ds_covparms=ds_covparms, ds_cellavg=ds_cellavg);
** analyses must be specified as: by name_config hypo sim , and use the variables as present in &ds_sim;
** below analyses are for independent subjects at baseline and at follow-up: 
	  subjects in the ctl arm and clusters in the intervention arm;  
title6 "post-test analysis";
	ods output ConvergenceStatus= _conv;
	ods output Covparms=_covparms;
	ods output SolutionF=_solutionf;
proc mixed data=&ds_sim(where=(time=1)); * only follow-up measurement;
by name_config hypo sim;
class clus subj t1_g1_class;
model outcome=tx / solution alpha=&alpha ddfm=kr;
random t1_g1 / subject=clus V Vcorr; * random slope for clusters in the intv arm at fu;
repeated / subject=subj(clus) group=t1_g1_class R Rcorr; *different residuals for intv arm at fu;
run;
	ods output close;
	* length statement to create same length for combining these datasets later;
	data _conv_pt; set _conv;length type_analysis $ 100; type_analysis="posttest";run;
	data _cov_pt; length covparm $ 20; length subject $ 20;set _covparms;  length type_analysis $ 100; type_analysis="posttest";run;
	data _fe_pt; length effect $ 20; set _solutionf; length type_analysis $ 100; type_analysis="posttest";if effect="tx";run;


title6 "ancova analysis on individual data";
	ods output ConvergenceStatus= _conv;
	ods output Covparms=_covparms;
	ods output SolutionF=_solutionF;
proc mixed data=&ds_sim; 
by name_config hypo sim;
class clus subj t1_g1_class;
model outcome=t1 t1*tx / solution alpha=&alpha ddfm=kr; 
random t1_g1 / subject=clus V Vcorr; * random slope for clusters in intervention arm at fu;
random intercept /subject=subj(clus) V Vcorr; *random intercept for subjects;
repeated / subject=time(subj*clus) type=vc group=t1_g1_class R Rcorr; *different residuals for intv arm at fu;
run;
	ods output close;
	* length statement to create same length for combining these datasets later;
	data _conv_anc_ind; set _conv;length type_analysis $ 100; type_analysis="ancova indiv";run;
	data _cov_anc_ind; length covparm $ 20;; length subject $ 20;;set _covparms; length type_analysis $ 100; type_analysis="ancova indiv";run;
	data _fe_anc_ind;length effect $ 20; set _solutionf; length type_analysis $ 100; type_analysis="ancova indiv";if effect="t1*tx";run;

title6 "ancova analysis on cluster means";
* calculate cluster-time means;
proc means data=&ds_sim mean var nway noprint; by name_config hypo sim tx; class clus time;
		var outcome t1; 
		output out=_cellavg mean=outcome_clusavg t1 var=outcome_variance t1_variance; run;
* analyze;
	ods output ConvergenceStatus= _conv;
	ods output Covparms=_covparms;
	ods output SolutionF=_solutionf;
proc mixed data=_cellavg;
by name_config hypo sim;
model outcome_clusavg= t1 t1*tx / solution ddfm=kr;
repeated / subject=clus type=csh group=tx; 
* what we loose is that type=cs in ctl arm ;
* and the (residual) variance is the same at baseline (both arms) and ct arm at fu;
* and at t=1: variance s2 (subjects) is related to variance tau2 by tau2=s2/n, with n the number of subjects;
run;
	ods output close;
	* length statement to create same length for combining these datasets later;
	data _conv_anc_clus; set _conv;length type_analysis $ 100; type_analysis="ancova clusavg";run;
	data _cov_anc_clus; length covparm $ 20;; length subject $ 20;set _covparms; length type_analysis $ 100; type_analysis="ancova clusavg";run;
	data _fe_anc_clus; length effect $ 20; set _solutionf; length type_analysis $ 100; type_analysis="ancova clusavg";if effect="t1*tx";run;

title6 " "; * clear title;

*write estimated parameters to output datasets;
data &ds_conv; set _conv_pt _conv_anc_ind _conv_anc_clus;run;
data &ds_feparms; set _fe_pt _fe_anc_ind _fe_anc_clus;run;
data &ds_covparms; set _cov_pt _cov_anc_ind _cov_anc_clus;run;
data &ds_cellavg; set _cellavg; drop t1_variance;run;

*clean-up;
/*proc datasets nolist; 
delete _: ; * delete temporary datasets; 
run; quit;*/

%mend analyze_outcomes;





%macro explore_configs(ds_config=, dir=%str(".";),  n_sim=, seed=37, verbose=0);

** the dataset ds_config must contain at least the variables
** n_sim=, n_clus0=,n_clus1=, n_subj0=, n_subj1=, 
 s2_c=, s2_ct0=, s2_s=, s2_st0=, s2_0ct1=, s2_0st1=, s2_1ct1=,s2_1st1=,
 ftime_t0=, ftime_t1=, delta=;

* possible code to check whether translation of correlations, sd to variances
* went ok;
* stop program if the conditions on correlations are not satisfied;
*if min(s2_0ct1,s2_0st1,s2_1ct1,s2_1st1)< 0 then do;
*	put "ERROR: correlations not in range to get positive variance components, execution stopped";
	* possibly add more details; 
*	stop;
*end;

* set the folder for the output dataset ;
libname dir &dir;
    
*** detemine how many configurations (i.e. how many records in &ds_config) ***;
data _null_;
    dsid= open("&ds_config");
    n_config= attrn(dsid,"nlobs");
    call symput('n_config', trim(left(n_config)));
run;    
%put n_config is &n_config;

*** iterate over the (macro) parameters from the configuration file;
%DO config=1 %TO &n_config;

   	*** get the macro parameters from the dataline with number &config**;
	data _null_; set &ds_config; 
		if _n_=&config;
     	array npar{*} _numeric_; * all numeric variables in an array;
     	do i=1 to dim(npar); * assign variables to macro variables with the same name;
        	call symput(vname(npar{i}), compress(trim(npar{i})));
     	end;
		call symput('name_config',compress(trim(name_config)) );* assign the name of the config; 
		datetime = put(datetime(), datetime19.);call symput('datetime',datetime);* time;
   	run;

  	* keep track of progress of this macro;
  	%put configuration &config of &n_config started at &datetime; 
	%put &name_config;

  	* suppress output, unless otherwise requested;
	%IF &verbose ne 0 %THEN %DO; %END;
	%ELSE %DO; ods exclude all; ods results off; options nonotes;%END;

  	* generate outcomes;
  	%sim_outcomes(ds_sim=ds_sim, seed=&seed, 
		n_sim=&n_sim, n_clus0=&n_clus0,n_clus1=&n_clus1, n_subj0=&n_subj0, n_subj1=&n_subj1, 
 		s2_c=&s2_c, s2_ct0=&s2_ct0, s2_s=&s2_s, s2_st0=&s2_st0, 
		s2_0ct1=&s2_0ct1, s2_0st1=&s2_0st1, s2_1ct1=&s2_1ct1,s2_1st1=&s2_1st1,
 		fe_t0=&fe_t0, fe_t1=&fe_t1, delta=&delta, 
		name_config=&name_config
		);
	
    *title6 "tx by time averages and variances";
	proc sort data=ds_sim; by hypo sim tx clus subj;
	proc means data=ds_sim(where=(hypo=0)) n mean std var; class tx time; var outcome;run;

	%analyze_outcomes(ds_sim=ds_sim, ds_feparms=dir.feparms_&name_config, 
                      ds_covparms=dir.covparms_&name_config,ds_conv=dir.conv_&name_config);

	* output on again;
	ods exclude none; ods results; options notes; 
%END;


%mend explore_configs;


** provide the configurations for which to calculate power, a.o. variance components have to be provided**;
data config;length name_config $ 100;
counter=.;
fe_t0=0; fe_t1=1;sigma_base=2.2; delta=1.3;
do n_clus1=11; 
do n_subj1=5;
n_clus0=n_clus1*n_subj1; 
n_subj0=1;
	do rho_1=0.05;
		do r=0.29;
		* partially nested design;
		* baseline;
		s2_base=sigma_base**2;
		s2_s=r * s2_base; 
		s2_st0=(1-r)* s2_base; 
		* ctl arm @fu;
		s2_0fu=s2_base;* equal residual variance in control arm at follow-up as at baseline;
		s2_0st1=s2_st0;
		s2_c=0; s2_ct0=0; s2_0ct1=0; 
		* intervention arm @fu;
		s2_1fu= s2_base * r / (r-rho_1); *equal test-retest assumption;
		s2_1ct1= rho_1*s2_1fu;
		s2_1st1=s2_1fu - s2_1ct1 - s2_s;
 		* calculation standard error;
		de_1= 1+ (n_subj1-1)*rho_1;
		A_1=de_1*s2_1fu - r**2 *s2_base;
		A_0 = (1-r**2)*s2_base;
		se2= A_1/(n_clus1*n_subj1) + A_0/(n_clus0*n_subj0); 
		se2_posttest=de_1*s2_1fu/(n_clus1*n_subj1) + s2_base/(n_clus0*n_subj0);
		* name configuration;
		counter+1;
		name_config=cat('config',counter);
		output;
end;end;end;end;
run;



* options mprint mlogic symbolgen;* macro debugging options on;
* options nomprint nomlogic nosymbolgen; * macro debugging options off;

** simulate and check power, for one selected config, check estimated random and fixed effects;;
%explore_configs(ds_config=config, n_sim=1000);


** post-processing to calculate rejection rates on H0 and H1;
libname dir ".";


data feparms;
set dir.feparms_: ;*combine feparms_... across all configurations; 
if probt ne . then reject=(probt < 0.05); 
se_emp=stderr;
run;

data conv;
set dir.conv_: ;* combine conv_... across all configurations;
run;

* combine the fixed effects and the non-converged indicator;
proc sort data=feparms; by name_config type_analysis hypo sim;run;
proc sort data=conv;by name_config type_analysis hypo sim; run;
* note the number of records is #name_config * #type_analysis * #hypo (=2) * #sim ; 
data feparms_conv; merge conv feparms;by name_config type_analysis hypo sim;run;

title "calculate rejection rate and estimated effect *ONLY* on converged analyses";
ods exclude all; ods results off;
ods output summary=summary_conv;
proc means data=feparms_conv n mean ; where status=0; * only the converged ones;
class name_config type_analysis hypo ; 
var reject  estimate se_emp; *reject_corr;
run;
ods exclude none; ods results;
ods output close;
* add the simulation parameters to the summary, using name_config (that has to be unique);
proc sort data=config; by name_config;run;
proc sort data=summary_conv; by name_config;run;
data summary_config_conv; merge config summary_conv; by name_config;run;


proc sort data=summary_config_conv; by n_clus1 n_subj1 n_clus0 rho_1 hypo type_analysis  ;
* write output summary;
data dir.summary_config_conv ; set summary_config_conv; se=sqrt(se2);se_posttest=sqrt(se2_posttest); run;
* print summary ;
proc print data=dir.summary_config_conv; by n_clus1 n_subj1 n_clus0 rho_1 ;
var  hypo type_analysis delta reject_n reject_mean  estimate_mean name_config 
     se se_posttest se_emp_mean;* reject_corr_mean;
run;


title "calculate non-convergence rate";
* mean of the non-convergence indicator;
proc sort data=conv; by name_config type_analysis hypo sim;run;
ods exclude all; ods results off;
ods output summary=summary_nonconv; 
proc means data=conv n mean; 
class name_config type_analysis hypo;
var status; *0= converged, 1= not converged;
run;
ods exclude none; ods results;
* add details of configurations ;
proc sort data=config; by name_config; run;
data summary_config_nonconv; merge config summary_nonconv; by name_config;run;

proc sort data=summary_config_nonconv; by n_clus1 n_subj1 n_clus0 rho_1 ;run;
* write output summary;
data dir.summary_config_nonconv; set summary_config_nonconv;run;
* print summary to .lst file;
proc print data=dir.summary_config_nonconv; by n_clus1 n_subj1 n_clus0 rho_1 ;
var  hypo type_analysis delta status_n status_mean;
run;
