*based on 150713 continuous simulation stepped wedge, adapted to run the trial by Mark vd Boogaard & Marieke Zeegers

for 5 clusters do the following: 
 ......

0	0	0	0	0	0	0	0	.	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	0	0	0	0	0	0	0	0	0	.	1	1	1	1	1	1	1	1	1	1	1	1	1
0	0	0	0	0	0	0	0	0	0	0	0	.	1	1	1	1	1	1	1	1	1	1	1
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	.	1	1	1	1	1
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	.	.	1 	1	1

so the sandwich design (groups in rows, time in columns) with implementation period denoted by a "."

note 5 clusters in total, 1 per sequencegroup

first generate a standard sw with two measurements per step:

...
0	0	0	0	0	0	0	0	0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	1	1	1	1	1	1	1	1	1
0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	1	1	1	1	1	1	1
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	1
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1 	1	1
and then make missings

Also:  
* different historical data (till time ...) from prospective data 
* different cluster-period sizes across the clusters
this is done by deleting from a maximum cluster-period size



* data generation depending on parameters;
* rho=ICC, rho_s, rho_c, n=number of subjects in a cluster, sigma_tot=1, n_clus and n_clus' (taken the same): number of clustes;
* fixed time effect tau:
* Note: 5 random effects have to be generated: cluster, cluster x time, subject, subject x time;
* we need expressions to calculate from rho, rho_s, rho_c, n, and sigma2_tot: the variance of these random effects;

* Needed;
* macro that generates, given a set of seeds and given the parameters;
* fixed time effect, fixed treatment effect, random effects for cluster, cluster x time, subject, subject x time, evaluation;
* a dataset;

%macro get_seeds(ds_seeds_source=, seedlength=  , seedstartrow=1, ds_initial_seeds=);

data hulp;
    set &ds_seeds_source;
    if ( (_N_ <= %eval(&seedlength+&seedstartrow-1)) and (_N_ >=&seedstartrow)) ;
run;

proc transpose data=hulp
               out= &ds_initial_seeds
               prefix=seed;
    var prime;
run;

data &ds_initial_seeds; set &ds_initial_seeds; keep seed1-seed&seedlength;run;

proc datasets nolist; delete hulp; run; quit; * clean up ;

%mend get_seeds;; 


%macro append(base=, data=);
proc datasets nolist;
        append base=&base  data=&data  force;;
run;
quit;
%mend append;


***** we need to describe the time effects (eg. via an array)that gets input from a dataset and is set just as with the seeds for simulation?***;
%macro simulate_data(n_clus=, n_subj=, n_rep=, data_seeds_in=, data_seeds_out=, data_time_effects=, data_sim=, sim_start=, sim_end=,
   delta_trt=, rho=, rho_c=, rho_s=,sigma2_tot=1, outcome=y);

/*
The macro "simulate_data" does the simulations with simulation counter from &start_sim to &end_sim;
input:    total number of clusters n_clus (in both groups together), how many of those are on treatment, 
		 number of subjects per cluster n_subj, number of repeated measurements n_rep,  
		number of clusters in a group for a stepped wedge design  defined in the data simulation step;
      		and their time effects (variables time_effect1-time_effect&n_rep, via a dataset)
			and their seeds (variables seed1-seed&seedlength, via a dataset)
          

output:
 1) a dataset of [n_clus + n_rep*n_clus + n_clus*n_subj+n_rep*n_clus*n_subj  ]=(n_rep+1)*n_clus * (1+*n_subj)  variables containing the -standard normal- residuals for 
     n_clus cluster effects, the n_rep*n_clus cluster x time effects, the n_clus*n_subj subject effects, the n_rep*n_clus*n_subj subject x time effects:
      
 2) a dataset of one record with the seeds of the last step (to be used for a next simulation)

Note:  &seedlength has to be [1+ n_clus + n_clus + n_clus*n_subj  ] =1+n_clus*(2+n_subj)
	      (1 seed for generating the cluster residuals, n_clus seeds for generating the subjects for each cluster,  
           n_clus seeds for  generating time effects for each cluster, n_clus*n_subj for generating the time effect for each subject)



/*** generating simulation dataset **************************************************/

* first determine how many seeds we need, which is also the length of the array seed in the dataset:
       (see above);
%local seedlength; %let seedlength=%eval(1+2*&n_clus+&n_subj*&n_clus);  

data &data_sim;
	/*** calculation of variances out of (auto) correlations***/
	sigma2_c_tot=&rho*&sigma2_tot; *total cluster variance;
	sigma2_c = &rho_c*&rho*&sigma2_tot; *cluster variance;
	sigma2_ct = (1-&rho_c)*&rho*&sigma2_tot; *cluster x time variance;

	sigma2_s_tot=(1-&rho)*&sigma2_tot; *total subject variance;
	sigma2_s= &rho_s*(1-&rho)*&sigma2_tot; *subject variance ;
	sigma2_st= (1-&rho_s)*(1-&rho)*&sigma2_tot; *subject x time variance;
			
	/*** initialization of the seeds***/   
  	set &data_seeds_in (keep=seed1-seed&seedlength); * only extract as many seeds as you need;
	/*** obtaining time effects **/ 
	set &data_time_effects (keep=time_effect1-time_effect&n_rep); * only extract as many fixed time effects you need;

	/** store seeds and time effects in an array for reference***/
	array seed{&seedlength} seed1-seed&seedlength; 
    array time_effect{&n_rep} time_effect1-time_effect&n_rep;

	/*** make &n_repcluster x time residuals using an array ***/
	array residuals_cluster_time{1:&n_rep};

		/** simulation ***/
		do sim=&sim_start to &sim_end;
			* generate a standardised cluster residual with seed{1}; 
			do cluster= 1  to &n_clus;
				* generate &n_rep cluster x time residuals for each cluster i 
				* using the seeds with seed{1+n_clus+1}... seed{1+n_clus+ i} ... seed{1+n_clus+n_clus}; 
				call rannor(seed{1}, res_cluster); 
	 			*;
				do time=1 to &n_rep; *time runs from 1 to &n_rep;
					call rannor(seed{1+&n_clus + cluster}, residuals_cluster_time{time}); *cluster x time residual for time=&t; 
	            end;* this do loop closed immediately;

				*generate n_subj subjects for each cluster i with using the seeds 
			    * seed{1+1} ..seed{1+i} ... seed{1+n_clus};
				do subject=1 to &n_subj;  
				call rannor(seed{1+cluster}, res_subject);

					* n_rep timepoints for each subject;
					do time=1 to &n_rep; 
					*generate subject x time residuals for each (cluster, subject) with remaining seeds; 
					call rannor(seed{1+2*&n_clus+ (cluster-1)*&n_subj + subject}, res_subject_time);
					* get the cluster x time residuals for this time point from the array filled above;
					res_cluster_time=residuals_cluster_time{time};  

						/***** build outcome separate records for H0 and H1 ****/
						time_class=time; * class variable to indicate random time effect;

						* this is the random part of the effect to investigate how well the correlations (ICC, auto-subject/cluster correlations) are described;
						random= +sqrt(sigma2_c)*res_cluster + sqrt(sigma2_ct)*res_cluster_time
											 +sqrt(sigma2_s)*res_subject + sqrt(sigma2_st)*res_subject_time;

						* sandwich stepped wedge design;
						* groups of equal size i.e. groupsize=&n_cluster/&n_rep;
								*group=cluster;  
						* treatment effect; 
						pre=8; * number of periods before switching begins;
						trt= (cluster=1)*(time >= pre+2) + (cluster=2)*(time >= pre +4) + (cluster=3)*(time >= pre +6)
							 + (cluster=4)*	(time >= pre +12) + (cluster=5)*(time >= pre +14); 
						* transition periods;
						transition=(cluster=1)*(time=pre+1)+ (cluster=2)*(time = pre+ 3) + (cluster=3)*(time = pre+ 5)
						      + (cluster=4)*(time = pre+ 11) + (cluster=5)*(time=pre+ 12 or time =pre+13);
						* varying clustersizes, use macro n_subj as the maximal cluster size;
						max_clustersize=(cluster=1)*( 40*(time <= pre-1)+8*(time >= pre) ) + 
					 					(cluster=2)*( 25*(time <= pre-1) + 8*(time >= pre) ) +
										(cluster=3)*( 100*(time <= pre-1) + 8*(time >= pre) ) +
										(cluster=4)*( 50*(time <= pre-1) + 20*(time >= pre) ) +
										(cluster=5)*( 15*(time <= pre-1) + 8*(time >= pre) );
						* outcomes;
						if (transition=0 and subject <= max_clustersize) then 
							do;
							hypo=0; &outcome  =  time_effect{time}+ random; output;* H0;
							hypo=1; &outcome  =  time_effect{time}+ random + trt*&delta_trt; output; *H1;
							end;
					end;
				end;
            end;
		end;
run;




/*** finalization ****************************************************************/

* keep the seeds of the last record of the simulation dataset (for use in another dataset);
data &data_seeds_out; 
	set &data_sim point=nobs nobs=nobs;
	keep seed1-seed&seedlength;
	output;
	stop; *else an infinite loop :-) ;
run; 

/* drop the seeds from the residuals dataset;
data &data_resid; set &data_resid; 
	drop seed1-seed&seedlength;
run;
*/
%mend simulate_data;

** analyze by using random effect for cluster (only) i.e. Hussey and Hughes model;
%macro re_analyze_repmeas(ds_in=, ds_out=, outcome=,sim_parameters=%str(;), 
proc_options = %str(method=reml), ddfm=satterth, alpha=0.05);

**  assumes a dataset with subject level data with variables:
        sim hypo cluster time trt &outcome;
** 	&sim_parameters can be used to add simulation parameters to the output dataset
    e.g. number of simulations, clusters, subjects, covariance structure etc and e.g. also predicted power;
** note that we fit with 'noint' option in proc mixed, so time effects are directly 
   estimated (not difference with respect to the baseline=intercept) **;

* sort to get the data corresponding to one hypothesis together;
* and keep only the relevant variables; 
* note subjects are just repeated record within hypo cluster time now with no indexing subject variable; 

proc sort data=_outcomes (keep= sim hypo cluster time trt &outcome) ; 
by sim hypo cluster time ;run;

* analyze data;
ods listing close; * no output to listing;
ods output ConvergenceStatus= _conv0;
ods output SolutionF=_fixed0;
ods output CovParms=_random0;
proc mixed data=_outcomes cl covtest &proc_options;
by sim hypo;
* note: because data sorted on cluster, can leave it out in class statement for efficiency;
class time;
model &outcome= time trt  / solution cl noint ddfm=&ddfm;
random intercept / subject=cluster; 
run;
ods output close; ods listing; quit;

*** checking convergence status **;
data _conv; length Reason $ 200; set _conv0 (drop= pdG pdH);
	call symput('converged', status); * echo convergence to log file;
	format status best12.; * to get decimals;
	rename status=noconv;
	* add simulation parameters if specified;
	&sim_parameters;
run;


* only keep the results for treatment and time effects and the 95%-CI for coverage;
data _fixed(drop= effect time estimate stderr df tvalue probt alpha lower upper); 
set _fixed0; 
by sim hypo;  
array time_est{1:&n_rep}; array time_se{1:&n_rep};
retain time_est1-time_est&n_rep 
       time_se1-time_se&n_rep 
       trt_est 	trt_se 	trt_p	trt_reject trt_se2
	   trt_ll trt_ul;
if compress(effect)='trt' then do; 	* for simulation summary;
	trt_est=estimate; trt_se=stderr; trt_p=probt;
	trt_reject=(trt_p < &alpha); trt_se2=(trt_se)**2; 
	trt_ll=lower; trt_ul=upper;
	end;
else if compress(effect)='time' then do; time_est{time}=estimate; time_se{time}=stderr;end;
if last.hypo then output;
run;

/* to sort out
* keep the variances of random effects and 95%-CI (approximated by Satterthwaite) for coverage;
* note: if we want to have also cluster x time effects and subject random effects, then adapt this;
data _random(drop=covparm subject estimate stderr zvalue probz alpha lower upper); 
set _random0; 
retain s2_c s2_c_ll s2_c_ul
       s2_s s2_s_ll s2_s_ul;
if lowcase(covparm)='intercept' and lowcase(subject)='cluster' then do; 
	s2_c=estimate; s2_c_ll=lower; s2_c_ul=upper;*s2_c_se=stderr;end;
if lowcase(covparm)='residual' then do; 
	s2_s=estimate; s2_s_ll=lower; s2_s_ul=upper;*s2_s_se;
	output;end;* last record so output now;
run;
and add to the merge below later
*/

* write the dataset with results;
* add coverage, note delta_trt has to available from &sim_parameters;
* add type of analysis;
data &ds_out; merge _conv _fixed ; by sim hypo; * note merge by sim and hypo;
	* if no convergence (could be because iteration stops), then no estimates;
	array estimates(*) trt: time: ;
	if noconv=1 then 
	do;
		do i=1 to dim(estimates); estimates(i)=.;end;
		drop i;
	end;
	* add coverage if applicable;
	if noconv=1 then do; cover_trt=.;end;
		else if hypo in (0) then 
			do; cover_trt=( (0 ge trt_ll) and (0 le trt_ul) );
			end;
		else if hypo in (1) then 
			do; cover_trt=( (delta_trt ge trt_ll) and (delta_trt le trt_ul) );
			end;
	analysis="reml";	*random effect analysis;
run;
%mend re_analyze_repmeas;

/* to adapt to incorporate;
* to adjust and implement, if unequal cluster sizes then perhaps use a weighting by inverse variance of cluster means;
**analysis on cluster means **;
%macro clusavg_analyze_repmeas(ds_in=, ds_out=, outcome=,sim_parameters=%str(;), 
proc_options = %str(method=reml), ddfm=kenwardroger, alpha=0.05);

**  assumes a dataset with subject level data with variables:
        sim hypo cluster time trt &outcome;
** 	&sim_parameters can be used to add simulation parameters to the output dataset
    e.g. number of simulations, clusters, subjects, covariance structure etc and e.g. also predicted power;
** note that we fit with 'noint' option in proc mixed, so time effects are directly 
   estimated (not difference with respect to the baseline=intercept) **;

* sort to get the data corresponding to one hypothesis together;
* keep only the relevant variables;
* note subjects are now just repeated record within hypo cluster time  with no indexing subject variable; 
proc sort data=_outcomes (keep= sim hypo cluster time trt &outcome) ; 
by sim hypo cluster time ;run;
* aggregate to cluster level;
proc means data=_outcomes mean nway noprint; class sim hypo cluster time ; 
var &outcome trt; 
output out=_outcomes_clusavg_long mean=clusavg_&outcome trt ; 
run;


* analyze the cluster average data;
ods listing close; * no output to listing;
ods output ConvergenceStatus= _conv0;
ods output SolutionF=_fixed0;
ods output CovParms=_random;
ods output Rcorr=_corr0;
proc mixed data=_outcomes_clusavg_long covtest &proc_options;
	by sim hypo;
	class time;
	model clusavg_&outcome= time trt  / solution cl noint ddfm=&ddfm;
	repeated / subject=cluster type=cs rcorr; * random intercept / subject=cluster; 
	* alternative to repeated time measurement within clusters via random effects;
	* the repeated statement allows however also other correlations such as autoregressive;
run;
ods output close; ods listing; quit;

*** checking convergence status **;
data _conv; length Reason $ 200; set _conv0 (drop= pdG pdH);
	call symput('converged', status); * echo convergence to log file;
	format status best12.; * to get decimals;
	rename status=noconv;
	* add simulation parameters if specified;
	&sim_parameters;
run;

* variances of the random effects;
data _random(drop=covparm subject estimate stderr zvalue probz alpha lower upper); 
set _random0; 
retain s2_c s2_c_ll s2_c_ul
       s2_s s2_s_ll s2_s_ul;
if lowcase(covparm)='intercept' and lowcase(subject)='cluster' then do; 
	s2_c=estimate; s2_c_ll=lower; s2_c_ul=upper;*s2_c_se=stderr;end;
if lowcase(covparm)='residual' then do; 
	s2_s=estimate; s2_s_ll=lower; s2_s_ul=upper;*s2_s_se;
	output;end;* last record so output now;
run;

* only keep the results for treatment and time effects and the 95%-CI for coverage;
data _fixed(drop= effect time estimate stderr df tvalue probt alpha lower upper); 
set _fixed0; 
by sim hypo;  
array time_est{1:&n_rep}; array time_se{1:&n_rep};
retain time_est1-time_est&n_rep 
       time_se1-time_se&n_rep 
       trt_est 	trt_se 	trt_p
	   trt_low trt_up;
if compress(effect)='trt' then do; 	* for simulation summary;
	trt_est=estimate; trt_se=stderr; trt_p=probt;
	trt_reject=(trt_p < &alpha); 
	trt_se2=(trt_se)**2; 
	* to calculate how much the model se is an under/over estimate of the true variance;
	if trt_se >0 then trt_se_inv=1/trt_se; else trt_se_inv=.;		
	trt_ll=lower; trt_ul=upper;
	end;
else if compress(effect)='time' then do; time_est{time}=estimate; time_se{time}=stderr;end;
if last.hypo then output;
run;

* write the dataset with results;
* add coverage, note delta_trt has to available from &sim_parameters;
data &ds_out; merge _conv _fixed _corr; by sim; 
	* if no convergence (could be because iteration stops), then no estimates;
	array estimates(*) trt: time: ;
	if noconv=1 then 
	do;
		do i=1 to dim(estimates); estimates(i)=.;end;
		drop i;
	end;
	* random effects not applicable for cluster average analysis;
	s2_c=.; s2_c_ll=.; s2_c_ul=.;	
	s2_s=.; s2_s_ll=.; s2_s_ul=.;	
	* add coverage, note delta_trt has to available from &sim_parameters;
	if noconv=1 then cover_trt=.;
		else if hypo in (0) then cover_trt=( (0 ge trt_ll) and (0 le trt_ul) );
		else if hypo in (1) then cover_trt=( (delta_trt ge trt_ll) and (delta_trt le trt_ul) );
	analysis="ca__"; *cluster average analysis;
run;
%mend clusavg_analyze_repmeas;
*/

%macro rep_meas(n_clus=,  
n_subj=, n_rep=, 
delta_trt=, data_time_effects=, rho=, rho_c=,rho_s=, outcome=y,sigma2_tot=1, 
n_sim=, sim_block=,seedstartrow=1,
dir='.', store=%str(work), ds_out=);

* store is the directory in which the temporary datasets 
(simulation datasets and the analysis record for each simulation dataset) are kept;
* ds_out is the dataset with all the analyses records of all the simulation datasets;
* ds_out._ is the dataset with power and type I error;


%local sim_parameters; %local n_step; 
%let sim_parameters=%str(n_sim=&n_sim;n_clus=&n_clus;n_subj=&n_subj;
n_rep=&n_rep;delta_trt=&delta_trt;rho=&rho;rho_c=&rho_c;rho_s=&rho_s);

%let n_step=%sysevalf(&n_sim/&sim_block,ceil);

* initialize seeds,  by default the file with seeds is in the current directory;
libname dir &dir;
data randomprimes; set dir.randomprimes;run;
%get_seeds(ds_seeds_source=randomprimes, seedlength=%eval(1+2*&n_clus+&n_subj*&n_clus), seedstartrow=&seedstartrow, ds_initial_seeds=_seeds_in);

%DO step=1 %TO &n_step;  

	* simulate &sim_block simulation datasets;
	%simulate_data(n_clus=&n_clus, n_subj=&n_subj, n_rep=&n_rep, data_seeds_in=_seeds_in, 
	data_seeds_out=_seeds_out, data_time_effects=_time_effects, data_sim=&store.._outcomes, 
	sim_start=(&step-1)*&sim_block + 1, 
	sim_end=(&step)*&sim_block, 
	delta_trt=&delta_trt,rho=&rho, rho_c=&rho_c, rho_s=&rho_s,sigma2_tot=&sigma2_tot, outcome=&outcome);

	* analyze these simulation datasets; 
	* using random effect for cluster;
	%re_analyze_repmeas(ds_in=&store.._outcomes, ds_out=&store.._analyses_re, outcome=&outcome,
	sim_parameters=&sim_parameters);
	/* using cluster averages: note that we have to define a different output dataset then;
	%clusavg_analyze_repmeas(ds_in=&store.._outcomes, ds_out=&store.._analyses_ca, outcome=&outcome,
	sim_parameters=&sim_parameters);
	*/
	* store the results, note two datasets concatenated;
	%IF &step=1 %THEN %DO;data &ds_out; set &store.._analyses_re /*&store.._analyses_ca */;run; %END;
	%ElSE %DO; %append(base=&ds_out, data=&store.._analyses_re); 
			   /*%append(base=&ds_out, data=&store.._analyses_ca);*/%END;

	* update the seeds for the next step;
	data _seeds_in; set _seeds_out;run;

%END; 

* summarize results: power and bias from simulation;
proc means data=&ds_out  noprint;
 class analysis hypo; id n_sim n_clus n_subj n_rep rho rho_c rho_s delta_trt;
 * no var statement so all variables are averaged;
 *var noconv trt_reject trt_est trt_se2 time_est1-time_est&n_rep r_clusavg cover_trt;
	 output out=&ds_out._(where=((hypo ne .) and analysis ne "") drop=_TYPE_ _FREQ_)
	        mean=  var(trt_est)=var_trt_est  n(trt_est)=n ; 
			*for all variables the average and for trt_est also the variance;
			* note that the file name can only be 36 characters long! ;
run;

/* add possibly the theoretical derived power;
data &ds_out._simresult; set &ds_out._simresult;
run;
*/

%mend;

%macro explore_configs(ds_config=, data_time_effects=, seedstartrow=1,dir=, ddfm=satterth, log=no, notify_sim=none);

*dir is the output dataset folder where the results of the simulated configurations are kept;
** the dataset ds_config must contain at least the variables
** n_clus n_subj n_rep delta_trt rho rho_c rho_s n_sim sim_block; 
** additionally an external dataset containing the time effects has to be specified (data_time_effects); 

    
*** detemine how many configurations (i.e. how many records in &ds_config) ***;
data _null_;
    dsid= open("&ds_config");
    n_config= attrn(dsid,"nlobs");
    call symput('n_config', trim(left(n_config)));
run;    
    
%put n_config is &n_config;

  *** iteratively get the running (macro) parameters from the configuration file;
%DO config=1 %TO &n_config;

   *** get the macro parameters from the dataline with number &config**;
   data _null_;
     set &ds_config; if _n_=&config;
     array npar{*} _numeric_; * all numeric variables in an array;
     do i=1 to dim(npar);
        call symput(vname(npar{i}), compress(trim(npar{i})));
     end;        
   run;

   		* keep track of progress of this macro;
   		data _null_;datetime = put(datetime(), datetime19.);call symput('datetime',datetime);run;
		%put configuration &config of &n_config started at &datetime; 


		* set the output disk;
		%IF &dir = %THEN %DO; libname dir "."; %END;%ELSE %DO; libname dir &dir; %END;
		*build output dataset name describing the parameters;
       %local ds_out; 
	   %local d_trt; %let d_trt=%sysevalf(&delta_trt*100, floor);* in terms of percentiles;
	   %local r;%let r=%sysevalf(100*&rho,floor); *correlation in term of percentiles;
	   %local rc;%let rc=%sysevalf(100*&rho_c,floor); *idem;
	   %local rs;%let rs=%sysevalf(100*&rho_s,floor); * idem;
	   * note that the file name can only be 36 characters long! ;
	   %let ds_out= dir.c&n_clus.s&n_subj.e&n_rep.tx&d_trt.r&r.rc&rc.rs&rs;
	   
	   *suppress the log of this simulation using the filename and proc printto statements;
	   filename junk dummy; proc printto log=junk;run;   
	   option nonotes;
		
		%rep_meas(n_clus=&n_clus, n_subj=&n_subj, n_rep=&n_rep, 
                  delta_trt=&delta_trt, 
                  data_time_effects=&data_time_effects,  
				  rho=&rho, rho_c=&rho_c,rho_s=&rho_s,  		
				n_sim=&n_sim, sim_block=&sim_block,ds_out=&ds_out);

 		* put normal log on;
		proc printto;run;
		option notes;		
%END;

* log the ending of this macro;
data _null_;datetime = put(datetime(), datetime19.);call symput('datetime',datetime);run;
%put macro explore_configs completed at &datetime; 
%mend explore_configs;

* debug ;
*options macrogen symbolgen source2;
options nomacrogen nosymbolgen nosource2;

* note that the variance of a single observation is set to 1: sigma2_tot=1, so delta_trt=effect size;
data configs; 
do n_clus=5; * new number of clusters;
do n_subj=100; * actually in the simulation this number is adapted to varying cluster size and missing by transition periods; 
do n_rep=24; 
do delta_trt=0.29; *0.074/0.25512; * new effect size;
do rho=0.035; do rho_c=1; do rho_s=0;  * new subjects at every time;
do n_sim=400; 
do sim_block=25;
	output;
end;end; end; end; end;end;end; end; end;
run;

libname dir ".";
* set fixed time effects;
data _time_effects; array time_effect{26}; 
	do i=1 to 26; time_effect{i}=(i-1)*0.001; drop i;* time effect at time=1 is set to 0;end; 
run; 


%explore_configs(ds_config=configs, data_time_effects=_time_effects);



/*** checks on simulation ***/
* get the last sim number;
data _null_; set configs;
call symput('last_sim',trim(n_sim));
run;


title "check: number of subjects (with non-missing outcome) per cluster-period";
proc tabulate data=_outcomes(where=(sim=&last_sim and hypo=0)); class cluster time; var y;
table cluster, time*y=" "*n;
run;
title "check: treatment assignment per cluster-period (with group)";
proc tabulate data=_outcomes(where=(sim=&last_sim and hypo=0)); class cluster time; var trt;
table cluster, time*trt*mean=" "*f=3.0;
run;
title "check: mean outcome per cluster-period";
proc tabulate data=_outcomes(where=(sim=&last_sim and hypo=0)); class cluster time; var y;
table cluster, time*y*mean;
run;






