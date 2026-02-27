/*** one proportion (success rate) per cluster, varying cluster sizes and successrates **/


ods rtf file="voorna_Dulmen_st260227.doc" style=minimal;

title "";
proc odstext;
p "for Spain: " /style=[fontsize=10pt just=l];
p "          fraction PPI users unknown...0.10 (lower end) taken;" / style=[fontsize=10pt just=l];
p "          fraction inappropriate unknown ... 0.20 (lower end) taken;" / style=[fontsize=10pt just=l];
p "          relative_effect unknown: 0.30 taken (lower end) taken; " / style=[fontsize=10pt just=l];
p "for Hungary: " /style=[fontsize=10pt just=l];
p "          relative effect take small given that elderly more digibeet...0.30 (lower end) taken;" / style=[fontsize=10pt just=l];
p "for other countries: " /style=[fontsize=10pt just=l];
p "          relative effect  ... 0.5 taken;" / style=[fontsize=10pt just=l];


run;

************** sample of clusters ****************************************************;

*** here drawn from the input distribution of averages ***************************;
title "averages per country";
data clus_input0;
length type  $30.;
input type $ av_clus_n av_frac_ppi av_frac_inappr av_frac_suc6;
datalines;
Spain 		1538 0.10 	0.20 0.3 
Hungary 	2048 0.20 	0.40 0.3 
Estonia 	1239 0.10 	0.20 0.3 
Netherlands 1694 0.085 	0.60 0.3 
;
run;
data clus_input0;set clus_input0;
capture=0.8; *only 8 months of the 12 months;
dropout=0.2; *some patients are lost to follow-up;
run;
proc print;run;

title "example generated with 70 to 130% around the country averages";
data clus_input; set clus_input0;
call streaminit(37); 
do totnumber_clus= rand('integer', 4,6); 
	do i=1 to totnumber_clus; 
	clus_n=rand('uniform',0.7,1.3)*av_clus_n;
	frac_ppi=rand('uniform',0.7,1.3)*av_frac_ppi;
	frac_inappr=rand('uniform',0.7,1.3)*av_frac_inappr;
	frac_suc6=rand('uniform',0.7,1.3)*av_frac_suc6;
    output;
end;end;
run;


****************** power calculation ***********************************************;
footnote "n_ppi= clus_n*frac_ppi*capture*(1-dropout)";
data clus_est; set clus_input;
n_ppi= clus_n*frac_ppi*capture*(1-dropout);
n_inappr=frac_inappr*n_ppi;
var_frac_suc6=(frac_suc6)*(1-frac_suc6)/n_inappr;
se_frac_suc6=sqrt(var_frac_suc6);
run;
/* proc sort data=clus_input; by type;
proc print;by type av_clus_n av_frac_ppi av_frac_inappr av_frac_suc6 totnumber_clus; run;
*/
footnote " ";

data total_se; set clus_est end=last; 
retain ivar_partialsum 0 est_partialsum 0;
ivar=1/var_frac_suc6 ; 							 * inverse variance of this cluster which is the weight;
iv_weighted_delta=frac_suc6*ivar;
est_partialsum= est_partialsum+  iv_weighted_delta;	*partial sum of weighted country estimates; 
ivar_partialsum= ivar_partialsum+ivar;      * ...........of inverse variances;
if last then do; 
	ivar_sum=ivar_partialsum; 
	est_meta=est_partialsum/ivar_sum; 
	var_meta=(ivar_sum)**(-1);
	se_meta=sqrt(var_meta);	
end;
run;

proc sort; by type;run; 
proc print; by type av_clus_n av_frac_ppi av_frac_inappr av_frac_suc6; 
var type clus_n frac_ppi frac_inappr n_ppi n_inappr frac_suc6 se_frac_suc6;
run;

data power; set total_se;where est_meta ne .;
do threshold=0, 0.10,0.20,0.25, 0.30,0.35,0.40,0.45;
	if (est_meta-threshold > 0) then power=probnorm((est_meta-threshold)/se_meta-1.96); 
	else power=. ; 
	output;
end; 
run;

proc print;var est_meta se_meta threshold power; run; 

ods rtf close; 


data power; 
* 18 subjects in 20 clusters with estimate at least 30%, so variance; 
do n_clus=4*4, 4*5, 4*6, 4*8; 
est=0.3;
var=0.3*0.7/(18*n_clus);
se=sqrt(var);
do threshold=0.20,0.25;
power=probnorm((est- threshold)/se-1.96); 
output;
end;end;
run;

proc print noobs; var est n_clus se threshold power; run;
