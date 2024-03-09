* https://www4.stat.ncsu.edu/~davidian/st732/examples/seizure_gee.sas ...
* analyse via log-binomial model or Poisson model fixed effects model (fixed effects for clusters;
* to get Relative Risk;

* try modified Poisson in proc genmod for robust standard errors;


* possibly adapt is there are different treatment effects over the clusters;

** description of analysis of Poisson: https://stats.idre.ucla.edu/r/faq/random-coefficient-poisson-models/ **;

* data generation depending on parameters;


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
%macro simulate_data(n_clus=, data_clussizes=, n_rep=, data_seeds_in=, data_seeds_out=, data_time_effects=, data_sim=, sim_start=, sim_end=,
   data_ctl_values=, data_trt_effects=,data_xotimes=, outcome=y);

/*
The macro "simulate_data" does the simulations with simulation counter from &start_sim to &end_sim;
input:    total number of clusters n_clus (in both groups together),
            trt assignment is set in the data step;
            ! as the simulation uses only fixed effects as clusters, the order of switching is randomized!
            and their cluster sizes (variables clussize1-clussize&n_clus, via a dataset)
            and their time effects (variables time_effect1-time_effect&n_rep, via a dataset)
            and their ctl effects (variables ctl_value1-ctl_value&n_clus, via a dataset): the values (without time trend that a
                cluster would have in the ctl condition;
            and their trt effects (variables trt_effect1-trt_effect&n_clus, via a dataset)
            and their seeds (variables seed1-seed&seedlength, via a dataset)


output:
 1) a dataset with sim hypo cluster time trt &outcome and additional auxiliary parameters;;
 2) a dataset of one record with the seeds of the last step (to be used for a next simulation)

/*** generating simulation dataset **************************************************/

* we need one seed for the random order of the start times of each cluster and a seed for each cluster to generate the subjects;
%local seedlength; %let seedlength=%eval(1+&n_clus);

data &data_sim;
/*** initialization of the seeds***/
  set &data_seeds_in (keep=seed1-seed&seedlength); * only extract as many seeds as you need;
/*** obtaining time effects **/
set &data_time_effects (keep=time_effect1-time_effect&n_rep); * only extract as many fixed time effects you need;
/** obtaining treatment effects **/
set &data_ctl_values (keep=ctl_value1-ctl_value&n_clus); * only extract as many trt effects (for each cluster) you need;
/** obtaining treatment effects **/
set &data_trt_effects (keep=trt_effect1-trt_effect&n_clus); * only extract as many trt effects (for each cluster) you need;
/** obtaining cluster sizes **/
set &data_clussizes (keep=clussize1-clussize&n_clus); *only extract as many cluster sizes as the number of cluster you need;
/** obtaining cross-over times of the clusters **/
set &data_xotimes (keep=xotime1-xotime&n_clus);


/** store seeds and time effects and trt effects and cluster sizes in an array for reference***/
array seed{&seedlength} seed1-seed&seedlength;
array time_effect{&n_rep} time_effect1-time_effect&n_rep;
array trt_effect{&n_clus} trt_effect1-trt_effect&n_clus;
array ctl_value{&n_clus} ctl_value1-ctl_value&n_clus;
array clussize{&n_clus} clussize1-clussize&n_clus;
* initialize a sequence of cross-over time times for the clusters;
array xotime{&n_clus} xotime1-xotime&n_clus; 

do sim=&sim_start to &sim_end;
* first generate a sequence of cross-over times for the clusters;
call ranperm(seed{1}, of xotime1-xotime&n_clus);
    do cluster= 1  to &n_clus;
    do time=1 to &n_rep;
        time_class=time; * class variable to indicate random time effect;
        * treatment indicator in stepped wedge design;
        trt= 1*(time >= xotime{cluster} ); *at and after cross-over time get the intervention;
            do subject=1 to clussize{cluster};
            /***** build outcome separate records for H0 and H1 ****/
            *step 1: the expected value of the outcome at the lowest level given covariates and higher level random effects;
            * is described Link(Expected outcome)= time_effects + trt/no trt + random cluster effect;
            * step 2: generate from a Poisson distribution with that mean using seed 1+cluster (a separate seed for each cluster);
            hypo=0; p= max( min( time_effect{time}+ctl_value{cluster}, 0.99999)  , 0.00001); * 0< probability < 1;
                    call ranbin(seed{1+cluster},1,p,&outcome);
                    output;
            hypo=1; p= max( min( time_effect{time}+ctl_value{cluster}+trt*trt_effect{cluster}, 0.99999)  , 0.00001); * 0< probability < 1;
                    call ranbin(seed{1+cluster},1,p,&outcome);
                    output;
            end;
        end;
end;
end;
run;

* incorporate missing data due to design (i.e., implementation periods) below;


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


** analyze by using fixed effect for cluster;
** as there is no random effect, maximum likelihood estimation is used;
%macro fe_analyze_repmeas(ds_in=, ds_out=, outcome=,sim_parameters=%str(;),
proc_options = %str(), alpha=0.05);


**  assumes a dataset with subject level data with variables:
        sim hypo cluster time trt &outcome;
**  &sim_parameters can be used to add simulation parameters to the output dataset
    e.g. number of simulations, clusters, subjects, covariance structure etc and e.g. also predicted power;
** note that we fit with 'noint' option in proc mixed, so time effects are directly
   estimated (not difference with respect to the baseline=intercept) **;

* sort to get the data corresponding to one hypothesis together;
* and keep only the relevant variables;
* note subjects are just repeated record within hypo cluster time now with no indexing subject variable;

proc sort data=&ds_in(keep= sim hypo cluster time trt &outcome) ;
by sim hypo cluster time ;run;

* analyze data;
    ods exclude all;ods noresults; * no output to listing;
    ods output ConvergenceStatus (nowarn)= _conv0;
    ods output ParameterEstimates=_fixed0;
    ods output CovParms (nowarn)=_random0; * nowarn as we will not estimate random effects;
proc glimmix data=&ds_in &proc_options;
by sim hypo;
* note: cluster as a categorical fixed effect;
* time is linear effect;
class cluster;
model &outcome (reference=first)= time trt cluster / dist=binary link=log solution cl ;
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
data _fixed(drop= effect estimate stderr df tvalue probt alpha lower upper);
set _fixed0;
by sim hypo;
retain time_est time_se
       trt_est  trt_se  trt_p   trt_reject trt_se2
       trt_ll trt_ul;
if compress(effect)='trt' then do;  * for simulation summary;
    trt_est=estimate; trt_se=stderr; trt_p=probt;
    trt_reject=(trt_p < &alpha); trt_se2=(trt_se)**2;
    trt_ll=lower; trt_ul=upper;
    end;
else if compress(effect)='time' then do; time_est=estimate; time_se=stderr;end;
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
    analysis="fixed effect &proc_options";
run;
%mend fe_analyze_repmeas;


%macro rep_meas(n_clus=,
n_rep=,data_clussizes=,data_ctl_values=, data_trt_effects=,data_time_effects=, data_xotimes=, outcome=y,
n_sim=, sim_block=,seedstartrow=1,
dir='.', store=%str(work), ds_out=);

* store is the directory in which the temporary datasets
(simulation datasets and the analysis record for each simulation dataset) are kept;
* ds_out is the dataset with all the analyses records of all the simulation datasets;
* ds_out._ is the dataset with power and type I error;

* keep track of progress of this macro;
data _null_;datetime = put(datetime(), datetime19.);call symput('datetime',datetime);run;
%put macro rep_meas started at &datetime;
* simulatie parameters;
%local sim_parameters; %local n_step;
%let sim_parameters=%str(n_sim=&n_sim;n_clus=&n_clus;n_rep=&n_rep;);
* default output file names;
%IF &ds_out eq %THEN %DO;
    %let ds_out=dir.c&n_clus.r&n_rep;
%END;


%let n_step=%sysevalf(&n_sim/&sim_block,ceil);

* initialize seeds,  by default the file with seeds is in the current directory;
* see macro %sim_data that needs: 1 seed for ordering of the start times, one seed for each cluster;
libname dir &dir;
data randomprimes; set dir.randomprimes;run;
%get_seeds(ds_seeds_source=randomprimes, seedlength=%eval(1+&n_clus), seedstartrow=&seedstartrow, ds_initial_seeds=_seeds_in);

%DO step=1 %TO &n_step;

    * simulate &sim_block simulation datasets;
    %simulate_data(n_clus=&n_clus, n_rep=&n_rep,data_clussizes=&data_clussizes,
     data_time_effects=&data_time_effects, data_ctl_values=&data_ctl_values, data_trt_effects=&data_trt_effects,
	 data_xotimes=&data_xotimes,
     data_seeds_in=_seeds_in, data_seeds_out=_seeds_out,
    sim_start=(&step-1)*&sim_block + 1, sim_end=(&step)*&sim_block, data_sim=&store.._outcomes,
     outcome=&outcome);

    * analyze these simulation datasets;
    * using random effect for cluster;
    %fe_analyze_repmeas(ds_in=&store.._outcomes, ds_out=&store.._analyses_re, outcome=&outcome,
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
%put macro rep_meas completed at &datetime;

%mend;

* clussize is within a year; 
data clus_info;length name $ 30; input name clussize ctl_value;
datalines;
1_Colon         	65  0.15
2_Rectal        	60  0.15
3_Esophageal    	80  0.547
4_Liver         	100 0.148
5_Pancreas      	78  0.33
6_HIPEC         	40  0.30
7_abdom_aneurysm 	35  0.58
8_EVAR				120	0.24
9_Nephrectomy   	60  0.25
10_Cystectomy   	78  0.38
11_craniotomy		55 	0.106
12_hip				50	0.10
13_limb				30	0.30
14_lung				60	0.20
15_thor_aneurysm	150	0.10
16_laryngectomy		65	0.30
17_mandibular		60	0.35
18_DIEP				118	0.22
19_hysterectomy		40	0.10
20_OVHIPEC			80	0.30
;
run;


title "*** total number of patients per year **";
proc means data=clus_info sum; var clussize;run;
/*
data clus_info;length name $ 30; input name clussize ctl_value;
datalines;
1_Esophageal    80  0.547
2_Pancreas      78  0.33
3_Liver         100 0.148
4_HIPEC         40  0.30
5_Rectum        60  0.15
6_Colon         65  0.15
7_Aneurysm      120 0.58
8_Nephrectomy   60  0.25
9_Cystectomy    78  0.38
;
run;
*/


* relative treatment effect;
data clus_info; set clus_info; trt_effect=-0.2*ctl_value;  run;

* ds treatment effects;
proc transpose data=clus_info out= data_trt_effects(drop=_name_) prefix=trt_effect;
var trt_effect; run;
* ds control values;
proc transpose data=clus_info out= data_ctl_values(drop=_name_) prefix=ctl_value;
var ctl_value; run;
* the cluster size per month; 
data clussize_month; set clus_info; clussize=round( 1*clussize/12); run;
proc transpose data=clussize_month out= data_clussizes(drop=_name_) prefix=clussize; var clussize;run;

* adapt the ds time effects (assume a change from -0.03 (-3%) over a 5 year (60 months) period;
* this is because the lowest prevalences are ~10% so the average (over all clusters) time trend
  cannot have to much reduction;
data data_time_effects; array time_effect{200} time_effect1-time_effect200;
do i=1 to 200; time_effect{i}= (i-1)* ( -0.03)/ 60;end;
* drift due to time from 0.25 to 0.20 over a 36 periodes;
* note that time effect on time=1 equals 0 so that (see simulation) time=1 gives
ctl_value{cluster} as baseline;
run;

*options symbolgen mlogic mprint;
options nosymbolgen nomlogic nomprint;


* number of clusters and the numer of months in a period;
data _NULL_;
    if 0 then set clus_info nobs=n;
    call symputx('n_clus',n);
    stop;
run;


**********************************************************************************;
** design: 24 pre, clusters switching in groups of 2 per month, so 10 switche times, 14 months post;  
%let n_rep=%eval(24+ 10+ 14); * 20 clusters with 10 switch times;

* times of crossover to intervention of each cluster;
data data_xotimes; * every two clusters cross-over together;
array xotime{&n_clus} xotime1-xotime&n_clus; do i=1 to &n_clus; xotime{i}= floor((i-1)/2)+1 +24 ; end;
run;

%rep_meas(n_sim=400, sim_block=100, n_clus=&n_clus,n_rep=&n_rep,data_clussizes=data_clussizes, data_ctl_values=data_ctl_values,
data_trt_effects=data_trt_effects,data_time_effects=data_time_effects, data_xotimes=data_xotimes, outcome=y,
seedstartrow=1,
dir='.', store=%str(work), ds_out=dir.pre24_c20g2_post14);




******************************************************************************;
** design: 24 pre, one cluster switching per month (20 months), 4 post; 
%let n_rep=%eval(24+ 20+ 4); 

* design: times of crossover to intervention of each cluster;
data data_xotimes; 
array xotime{&n_clus} xotime1-xotime&n_clus; do i=1 to &n_clus; xotime{i}= i +24 ; end;
run;

%rep_meas(n_sim=400, sim_block=100, n_clus=&n_clus,n_rep=&n_rep,data_clussizes=data_clussizes, data_ctl_values=data_ctl_values,
data_trt_effects=data_trt_effects,data_time_effects=data_time_effects, data_xotimes=data_xotimes, outcome=y,
seedstartrow=1,
dir='.', store=%str(work), ds_out=dir.Pre24_c20g1_post4);



******************************************************************************;
** design: only 12 months pre, 20 months switching, 4 months post; 
%let n_rep=%eval(12+ 20+ 4); 

* design: times of crossover to intervention of each cluster;
data data_xotimes; 
array xotime{&n_clus} xotime1-xotime&n_clus; do i=1 to &n_clus; xotime{i}= i +12 ; end;
run;

%rep_meas(n_sim=100, sim_block=100, n_clus=&n_clus,n_rep=&n_rep,data_clussizes=data_clussizes, data_ctl_values=data_ctl_values,
data_trt_effects=data_trt_effects,data_time_effects=data_time_effects, data_xotimes=data_xotimes, outcome=y,
seedstartrow=1,
dir='.', store=%str(work), ds_out=dir.Pre12_c20g1_post4);




***** with 40% drop-out en flooring ipv rounding, say 48-20-36 design***************************;
** design: 48 pre, one cluster switching per month (20 months), 36 post; 

%let n_rep=%eval(48+ 20+ 36); 

* the cluster size per month, 60% of it, minimum of 2 (for cluster with 30 per year); 
data clussize_month_reduced; set clus_info; clussize=floor( (0.6)*clussize/12); run;
proc transpose data=clussize_month_reduced out= data_clussizes(drop=_name_) prefix=clussize; var clussize;run;

* design: times of crossover to intervention of each cluster;
data data_xotimes;  
array xotime{&n_clus} xotime1-xotime&n_clus; do i=1 to &n_clus; xotime{i}= i +48 ; end;
run;

%rep_meas(n_sim=100, sim_block=100, n_clus=&n_clus,n_rep=&n_rep,data_clussizes=data_clussizes, data_ctl_values=data_ctl_values,
data_trt_effects=data_trt_effects,data_time_effects=data_time_effects, data_xotimes=data_xotimes, outcome=y,
seedstartrow=1,
dir='.', store=%str(work), ds_out=dir.Pre48_c20g1_post36_incl_6_10);


***** with 40% drop-out en flooring ipv rounding, say 48-20-24 design***************************;
** design: 48 pre, one cluster switching per month (20 months), 24 post; 

%let n_rep=%eval(48+ 20+ 24); 

* the cluster size per month, 60% of it, minimum of 2 (for cluster with 30 per year); 
data clussize_month_reduced; set clus_info; clussize=floor( (0.6)*clussize/12); run;
proc transpose data=clussize_month_reduced out= data_clussizes(drop=_name_) prefix=clussize; var clussize;run;

* design: times of crossover to intervention of each cluster;
data data_xotimes;  
array xotime{&n_clus} xotime1-xotime&n_clus; do i=1 to &n_clus; xotime{i}= i +48 ; end;
run;

%rep_meas(n_sim=100, sim_block=100, n_clus=&n_clus,n_rep=&n_rep,data_clussizes=data_clussizes, data_ctl_values=data_ctl_values,
data_trt_effects=data_trt_effects,data_time_effects=data_time_effects, data_xotimes=data_xotimes, outcome=y,
seedstartrow=1,
dir='.', store=%str(work), ds_out=dir.Pre48_c20g1_post24_incl_6_10);


***** with 1/3 drop-out, say 36-20-36 design***************************;
** design: 36 pre, one cluster switching per month (20 months), 36 post; 

%let n_rep=%eval(36+ 20+ 36); 

* the cluster size per month, 66% of it, minimum of 2 (for cluster with 30 per year); 
data clussize_month_reduced; set clus_info; clussize=round( (2/3)*clussize/12); run;
proc transpose data=clussize_month_reduced out= data_clussizes(drop=_name_) prefix=clussize; var clussize;run;

* design: times of crossover to intervention of each cluster;
data data_xotimes;  
array xotime{&n_clus} xotime1-xotime&n_clus; do i=1 to &n_clus; xotime{i}= i +36 ; end;
run;

%rep_meas(n_sim=100, sim_block=100, n_clus=&n_clus,n_rep=&n_rep,data_clussizes=data_clussizes, data_ctl_values=data_ctl_values,
data_trt_effects=data_trt_effects,data_time_effects=data_time_effects, data_xotimes=data_xotimes, outcome=y,
seedstartrow=1,
dir='.', store=%str(work), ds_out=dir.Pre36_c20g1_post36_incl_2_3);



**** with 1/4 dropout, say 36-20-36 design    **************************************;
** clusters switching each separately: 36 months follow-up and 36 pre and 75% of cluster size; 
%let n_rep=%eval(36+ 20+ 36); 

data clussize_month_reduced; set clus_info; clussize=round( (3/4)*clussize/12); run;
proc transpose data=clussize_month_reduced out= data_clussizes(drop=_name_) prefix=clussize; var clussize;run;

* design: times of crossover to intervention of each cluster;
data data_xotimes;  
array xotime{&n_clus} xotime1-xotime&n_clus; do i=1 to &n_clus; xotime{i}= i +36 ; end;
run;

%rep_meas(n_sim=100, sim_block=100, n_clus=&n_clus,n_rep=&n_rep,data_clussizes=data_clussizes, data_ctl_values=data_ctl_values,
data_trt_effects=data_trt_effects,data_time_effects=data_time_effects, data_xotimes=data_xotimes, outcome=y,
seedstartrow=1,
dir='.', store=%str(work), ds_out=dir.Pre36_c20g1_post36_incl_3_4);





****for comparison an ancova design to show that this does not yield sufficient power *************;
** ancova design: 24 pre, 10 clusters switch, after 12 months rest 10 cluster switch, then 12 months follow-up; 
%let n_rep=%eval(24+ 24 ); 

* the cluster size per month, 66% of it, minimum of 2 (for cluster with 30 per year); 
data clussize_month; set clus_info; clussize=round( clussize/12); run;
proc transpose data=clussize_month out= data_clussizes(drop=_name_) prefix=clussize; var clussize;run;

* design: times of crossover to intervention of each cluster;
data data_xotimes;  
array xotime{&n_clus} xotime1-xotime&n_clus; 
do i=1 to &n_clus; 
if i <= 10 then xotime{i}=24+1;
end;
run;

%rep_meas(n_sim=100, sim_block=100, n_clus=&n_clus,n_rep=&n_rep,data_clussizes=data_clussizes, data_ctl_values=data_ctl_values,
data_trt_effects=data_trt_effects,data_time_effects=data_time_effects, data_xotimes=data_xotimes, outcome=y,
seedstartrow=1,
dir='.', store=%str(work), ds_out=dir.Pre24_c20ancova_24months);




