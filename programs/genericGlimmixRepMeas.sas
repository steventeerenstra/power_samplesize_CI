
%macro append(base=, data=);
proc datasets nolist;
        append base=&base  data=&data  force;;
run;
quit;
%mend append;


%macro simulate_data(n_clus=, n_rep=, streaminit=, data_time_effects=, data_sim=, sim_start=, sim_end=,
   delta_trt=, sigma2_c=);

/*
The macro "simulate_data" does the simulations with simulation counter from &start_sim to &end_sim;
cluster > subject > time
new subjects in each period, so ***cross-sectionally***;

input:    total number of clusters n_clus (in both groups together), number of measurements n_rep
		  number of clusters in a sequence  and timing of intervention (stepped wedge design)  defined in the data simulation step;
      		and their time effects (variables time_effect1-time_effect&n_rep, via a dataset)
           

output:
a dataset with the simulated outcomes in the variable "outcome";      


/*** generating simulation dataset **************************************************/

data &data_sim;
/*** variance parameters **/
sigma2_c=&sigma2_c; *variance of cluster random effect on the link-transformed scale (here: log);

/*** initialization of the seed***/   
call streaminit(&streaminit);
/*** obtaining time effects **/ 
set &data_time_effects (keep=time_effect1-time_effect&n_rep); * only extract as many fixed time effects you need;
array time_effect{&n_rep} time_effect1-time_effect&n_rep;

do sim=&sim_start to &sim_end;
	do cluster= 1  to &n_clus;
	*cluster level random effect;
	random_c=sqrt(sigma2_c)*rand('normal');
	%clussize_cluster();* code setting the variable "clussize", possibly varying at baseline, so across clusters; 

do time=1 to &n_rep; 		
		time_class=time; * class variable to indicate random time effect;	
		* treatment indicator in stepped wedge design;
		%trt_clusterTime(); * code to set the variable trt depending on cluster, time;
		%clussize_ClusterTimeTrt(); *code for setting the variable clussize depending on cluster, time,trt, possibly if time varying; 
			do subject=1 to clussize;
			/***** build outcome separate records for H0 and H1 ****/
			hypo=0; %data_H0(); output;
			hypo=1; %data_H1(); output;
			end;
		end;
	end;
end;
run;

* incorporate missing data due to design (i.e., implementation periods);
%missing_data();

%mend simulate_data;


** analyze by using random effect for cluster (only) i.e. Hussey and Hughes model;
%macro glimmix_repmeas(ds_in=, ds_out=, sim_parameters=%str(;), 
						  proc_options = %str(method=laplace), alpha=0.05);

** other options would be: method=quad(qpoints=100)
**  assumes a dataset with subject level data with variables: sim hypo cluster time trt outcome;
** 	&sim_parameters can be used to add simulation parameters to the output dataset
    e.g. number of simulations, clusters, subjects, covariance structure etc and e.g. also predicted power;
**  depending on time as categorical, continuous, or absent, if timenote that we fit with 'noint' option in proc mixed, so time effects are directly 
   estimated (not difference with respect to the baseline=intercept) **;

* sort to get the data corresponding to one hypothesis together;
* and keep only the relevant variables; 
* note subjects are just repeated record within hypo cluster time now with no indexing subject variable; 
proc sort data=&ds_in(keep= sim hypo cluster time trt outcome) ; 
by sim hypo cluster time ;run;

* analyze data;
	ods exclude all;ods noresults; * no output to listing;
	ods output ConvergenceStatus= _conv0;
	ods output ParameterEstimates=_fixed0;
	ods output CovParms=_random0;
proc glimmix data=&ds_in &proc_options;
by sim hypo;
* note: because data sorted on cluster, can leave it out in class statement for efficiency;
* time as continous, so no class variable, but then the model with an intercept;
%glimmix_model();  * fixed and random effects or residual; 
run;
	ods output close;  quit;
	ods exclude none; ods results;* output visual;

*** checking convergence status **;
data _conv; length Reason $ 200; set _conv0;
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
retain trt_reject trt_est trt_ll trt_ul trt_se 	trt_p trt_se2 
	   time_est1-time_est&n_rep time_se1-time_se&n_rep ;
array time_est{1:&n_rep}; array time_se{1:&n_rep};
time = time ; * creates time as . if it does not exists (i.e. if time was not a class variable in the analysis;
if compress(effect)='trt' then do; 	* for simulation summary;
	trt_est=estimate; trt_se=stderr; trt_p=probt;
	trt_reject=(trt_p < &alpha); trt_se2=(trt_se)**2; 
	trt_ll=lower; trt_ul=upper;
	end;
else if compress(effect)='time' then do; 
	if time ne . then do;time_est{time}=estimate; time_se{time}=stderr;end; * if time class var. in analysis;
	else do; time_est1=estimate; time_se1=stderr; end; * if time continuous var. in the analysis;
	end; 
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
	/*
	* add coverage if applicable;
	if noconv=1 then do; cover_trt=.;end;
		else if hypo in (0) then 
			do; cover_trt=( (0 ge trt_ll) and (0 le trt_ul) );
			end;
		else if hypo in (1) then 
			do; cover_trt=( (delta_trt ge trt_ll) and (delta_trt le trt_ul) );
			end;
	*/
	analysis="&proc_options";
run;
%mend re_analyze_repmeas;


%macro sim_cluster_repmeas(n_clus=, n_rep=, 
delta_trt=, data_time_effects=, sigma2_c=1, 
n_sim=, sim_block=,seedstartrow=1,
dir='.', store=%str(work), ds_out=, proc_options = %str(method=laplace));

* hybrid simulation: by-processing in blocks;
* store is the directory in which the temporary datasets 
(simulation datasets and the analysis record for each simulation dataset) are kept;
* ds_out is the dataset with all the analyses records of all the simulation datasets;
* ds_out._ is the dataset with power and type I error;

* keep track of progress of this macro;
data _null_;datetime = put(datetime(), datetime19.);call symput('start_datetime',datetime);run;
%put macro sim_study started at &start_datetime; 
* simulatie parameters;
%local sim_parameters; %local n_step; 
%let sim_parameters=%str(n_sim=&n_sim;n_clus=&n_clus;
n_rep=&n_rep;delta_trt=&delta_trt;sigma2_c=&sigma2_c;);
* default output file names;
%IF &ds_out eq %THEN %DO;
	data _null_; delta_trt_100=abs(round(100*&delta_trt));
	call symput('delta_trt_100', strip(put(delta_trt_100,8.)) ); 
	run;
	%let ds_out=dir.c&n_clus._r&n_rep._tx&delta_trt_100;
%END;

%let n_blocks=%sysevalf(&n_sim/&sim_block,ceil);* number of blocks;

%DO step=1 %TO &n_blocks;  

	* simulate &sim_block times for one block, use streaminit=&step to get independent streams;
	%simulate_data(streaminit=&step, n_clus=&n_clus, n_rep=&n_rep, data_time_effects=&data_time_effects, 
    data_sim=&store.._outcomes, 
	sim_start=(&step-1)*&sim_block + 1, 
	sim_end=(&step)*&sim_block, 
	delta_trt=&delta_trt, sigma2_c=&sigma2_c);

	* analyze these simulation datasets; 
	* using random effect for cluster;
	%glimmix_repmeas(ds_in=&store.._outcomes, ds_out=&store.._analyses_re, 
	sim_parameters=&sim_parameters, proc_options=&proc_options);
	/* using cluster averages: note that we have to define a different output dataset then;
	%clusavg_analyze_repmeas(ds_in=&store.._outcomes, ds_out=&store.._analyses_ca, outcome=&outcome,
	sim_parameters=&sim_parameters);
	*/
	* store the results, note two datasets concatenated;
	%IF &step=1 %THEN %DO;data &ds_out; set &store.._analyses_re /*&store.._analyses_ca */;run; %END;
	%ElSE %DO; %append(base=&ds_out, data=&store.._analyses_re); 
			   /*%append(base=&ds_out, data=&store.._analyses_ca);*/%END;
%END; 

* summarize results: power and bias from simulation;
proc means data=&ds_out  noprint;
 class analysis hypo; id n_sim n_clus n_rep ;
 * no var statement so all variables are averaged;
 *var noconv trt_reject trt_est trt_se2 time_est1-time_est&n_rep r_clusavg cover_trt;
	 output out=&ds_out._(where=((hypo ne .) and analysis ne "") drop=_TYPE_ _FREQ_)
	        mean=  var(trt_est)=var_trt_est  n(trt_est)=n ; 
			*for all variables the average and for trt_est also the variance;
			* note that the file name can only be 36 characters long! ;
run;

* log the ending of this macro;
data _null_;datetime = put(datetime(), datetime19.);call symput('datetime',datetime);run;
%put macro rep_meas started at &start_datetime and completed at &datetime; 

%mend;



%macro missing_data();* missing data by design (e.g. implementation periods) depending on "cluster", "time"; %mend; 
%macro clussize_cluster();
/* code to set the variable "clussize" if constant over time but possibly varying over clusters*/
/* %local avg_clussize; %local cv_clussize;%let avg_clussize=10;%let cv_clussize=0;
%if &cv_clussize ne 0 %then
clussize=round( &avg_clussize * rand('gamma', 1/&cv_clussize**2) * &cv_clussize**2) ;
%else clussize=&avg_clussize;
*/
%mend;
%macro clussize_clusterTimeTrt(); 
*code for setting "clussize" depending on "cluster", "time","trt", possibly if time varying; 
clussize=30*(trt=0) + 49*(trt=1);
%mend;
%macro data_H0(); * data generation depending ont time_effect{time} and random effect of cluster: random_c;
p= max(0.00001, min(time_effect{time}+random_c, 0.99999) ); * no  probability equal or above 1 or below or equal 0;
outcome= rand('bernouilli',p);
%mend;
%macro data_H1(); * data generation depending ont time_effect{time} and random effect of cluster: random_c;	
p=max(0.00001, min(time_effect{time}+random_c+trt*&delta_trt, 0.99999) ); * no  probability equal or above 1;
outcome=rand('bernoulli',p);
%mend;		

data data_time_effects; 
array time_effect{10} time_effect1-time_effect10; 
do i=1 to 10; time_effect{i}= 0.23+ (i-1)*( 0.001 )/ 10;
* drift due to time  with +1% over 10 years;
end;run; 

%macro glimmix_model(); 
* glimmix model statement, at least: model outcome = trt <time>/ solution cl;
* time can be specified as class (then specify "noint" in the model) or not or even be left out;
/* binary error, identity link with time as categorical
class time;
model outcome (reference=first)= trt / dist=binary link=identity solution cl noint ;
random intercept / subject=cluster; */
/* binary error, identity link with time as continuous 
model outcome (reference=first)= trt time / dist=binary link=identity solution cl ;
random intercept / subject=cluster; */
/* binary error, identity link without time, i.e. no time trend 
model outcome (reference=first)= trt / dist=binary link=identity solution cl ;
random intercept / subject=cluster; */
%mend;


*for debugging;
options symbolgen mlogic mprint; 
options nosymbolgen nomlogic nomprint;

*********************** 3 measurements SW *******************************************************************;
%macro trt_clusterTime();
*code to define "trt" depending on "cluster","time";
trt=(time > cluster);
%mend;
* with time as continuous;
%macro glimmix_model(); 
* glimmix model statement, at least: model outcome = trt <time>/ solution cl;
* time can be specified as class (then specify "noint" in the model) or not or even be left out;
/* binary error, identity link with time as continuous */
model outcome (reference=first)= trt time / dist=binary link=identity solution cl ;
random intercept / subject=cluster; 
%mend;

%sim_cluster_repmeas(n_clus=2, n_rep=3, delta_trt=0.22, data_time_effects=data_time_effects,
sigma2_c=( 0.01 )**2, proc_options=%str(method=laplace),
n_sim=400, sim_block=200,seedstartrow=1, dir='.', store=%str(work), ds_out=dir.c2_r3_tx22_timecont);


* with time as class;
%macro glimmix_model(); 
* glimmix model statement, at least: model outcome = trt <time>/ solution cl;
* time can be specified as class (then specify "noint" in the model) or not or even be left out;
/* binary error, identity link with time as continuous */
class time;
model outcome (reference=first)= trt time / dist=binary link=identity solution cl noint;
random intercept / subject=cluster; 
%mend;

%sim_cluster_repmeas(n_clus=2, n_rep=3, delta_trt=0.22, data_time_effects=data_time_effects,
sigma2_c=( 0.01 )**2, proc_options=%str(method=laplace),
n_sim=400, sim_block=200,seedstartrow=1, dir='.', store=%str(work), ds_out=dir.c2_r3_tx22_timeclass);

* with time assumed to absent;
%macro glimmix_model(); 
* glimmix model statement, at least: model outcome = trt <time>/ solution cl;
* time can be specified as class (then specify "noint" in the model) or not or even be left out;
/* binary error, identity link with time as continuous */
model outcome (reference=first)= trt / dist=binary link=identity solution cl ;
random intercept / subject=cluster; 
%mend;

%sim_cluster_repmeas(n_clus=2, n_rep=3, delta_trt=0.22, data_time_effects=data_time_effects,
sigma2_c=( 0.01 )**2, proc_options=%str(method=laplace),
n_sim=400, sim_block=200,seedstartrow=1, dir='.', store=%str(work), ds_out=dir.c2_r3_tx22_notimet);



********************** with 5 measurements (2 retrospective) *******************************;
%macro trt_clusterTime();
*code to define "trt" depending on "cluster","time";
trt=(time > cluster+2);
%mend;

* with time as continuous;
%macro glimmix_model(); 
* glimmix model statement, at least: model outcome = trt <time>/ solution cl;
* time can be specified as class (then specify "noint" in the model) or not or even be left out;
/* binary error, identity link with time as continuous */
model outcome (reference=first)= trt time / dist=binary link=identity solution cl ;
random intercept / subject=cluster; 
%mend;

%sim_cluster_repmeas(n_clus=2, n_rep=5, delta_trt=0.22, data_time_effects=data_time_effects,
sigma2_c=( 0.01 )**2, proc_options=%str(method=laplace),
n_sim=400, sim_block=200,seedstartrow=1, dir='.', store=%str(work), ds_out=dir.c2_r5_tx22_timecont);

** with time as class; 
%macro glimmix_model(); 
* glimmix model statement, at least: model outcome = trt <time>/ solution cl;
* time can be specified as class (then specify "noint" in the model) or not or even be left out;
/* binary error, identity link with time as continuous */
class time;
model outcome (reference=first)= trt time / dist=binary link=identity solution cl noint;
random intercept / subject=cluster; 
%mend;

%sim_cluster_repmeas(n_clus=2, n_rep=5, delta_trt=0.22, data_time_effects=data_time_effects,
sigma2_c=( 0.01 )**2, proc_options=%str(method=laplace),
n_sim=400, sim_block=400,seedstartrow=1, dir='.', store=%str(work), ds_out=dir.c2_r5_tx22_timeclass);






* simulation check;
data a; set dir.c2_r5_tx22_timecont;
x=trt_est; y=lag10(x);
run;
proc sgpanel data=a; panelby hypo; scatter x=x y=y;run;

