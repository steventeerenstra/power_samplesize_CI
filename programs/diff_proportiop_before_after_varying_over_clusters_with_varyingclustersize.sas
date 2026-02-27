title "";
proc odstext;
p "for Spain: " /style=[fontsize=20pt just=l];
p "          fraction PPI users unknown...0.10 (lower end) taken;" / style=[fontsize=20pt just=l];
p "          fraction inappropriate unknown ... 0.20 (lower end) taken;" / style=[fontsize=20pt just=l];
p "          relative_effect unknown: 0.30 taken (lower end) taken; " / style=[fontsize=20pt just=l];
p "for Hungary: " /style=[fontsize=20pt just=l];
p "          relative effect take small given that elderly more digibeet...0.30 (lower end) taken;" / style=[fontsize=20pt just=l];
p "for other countries: " /style=[fontsize=20pt just=l];
p "          relative effect  ... 0.5 taken;" / style=[fontsize=20pt just=l];


run;

************** sample of clusters ****************************************************;

*** here drawn from the input distribution of averages ***************************;
title "averages per country";
data clus_input0;
length type  $30.;
input type $ av_clus_n av_frac_ppi av_frac_inappr av_relative_effect;
datalines;
Spain 		1538 0.10 	0.20 0.3 
Hungary 	2048 0.20 	0.40 0.3 
Estonia 	1239 0.10 	0.20 0.3 
Netherlands 1694 0.085 	0.60 0.3 
;
run;
proc print;run;

title "example generated";
data clus_input; set clus_input0;
call streaminit(37); 
do n_clus= rand('integer', 4,6); 
	do i=1 to n_clus; 
	clus_n=rand('uniform',0.7,1.3)*av_clus_n;
	frac_ppi=rand('uniform',0.7,1.3)*av_frac_ppi;
	frac_inappr=rand('uniform',0.7,1.3)*av_frac_inappr;
	relative_effect=rand('uniform',0.7,1.3)*av_relative_effect;
    output;
end;end;
run;

proc print;run;


****************** power calculation ***********************************************;

data clus_est; set clus_input;
n_ppi= clus_n*frac_ppi;
n_inappr=frac_inappr*n_ppi;
before_inappr=frac_inappr;
var_before=before_inappr*(1-before_inappr)/n_inappr;
after_inappr=before_inappr*(1-relative_effect);
var_after=after_inappr*(1-after_inappr)/n_inappr;
delta=before_inappr- after_inappr;
var_delta=var_before+ var_after;
run;
proc print;run;


data total_se; set clus_est end=last; 
retain ivar_partialsum 0 est_partialsum 0;
ivar=1/var_delta ; 							 * inverse variance of this cluster which is the weight;
iv_weighted_delta=delta*ivar;
est_partialsum= est_partialsum+  iv_weighted_delta;	*partial sum of weighted country estimates; 
ivar_partialsum= ivar_partialsum+ivar;      * ...........of inverse variances;
if last then do; 
	ivar_sum=ivar_partialsum; 
	est_meta=est_partialsum/ivar_sum; 
	var_meta=(ivar_sum)**(-1); 
	power=probnorm(est_meta/sqrt(var_meta)-1.96); 
end;
run;

proc sort; by type;run; 
proc print; by type av_clus_n av_frac_ppi av_frac_inappr av_relative_effect; 
var type clus_n frac_ppi frac_inappr relative_effect n_ppi n_inappr before_inappr after_inappr delta var_delta est_meta var_meta power;
run;




