** sample size for log-normal distribution, e.g. PK bio-equivalence trials **;

data power_equivalence;
alpha2=0.10; *two-sided alpha, so 0.10 corresponds to 90% CI;
cv=0.25; * coefficient of variation of the *within-subject* variation sd/mean on original scale; 
r0=1.05; * expected ratio test vs reference on original scale;
n1=26; *number of subjects in sequence1;
rmargin_up=1.25; * upper margin on ratio on the original scale;
rmargin_low=0.8; * lower margin ..;
** assuming log-normal distribution;
sigma_log=sqrt( log(cv**2 + 1) ) ; *cv=sqrt(exp(sigma^2) -1 );
delta_up=log(rmargin_up) - log(r0); * distance between upper and r0;
delta_low=log(r0)-log(rmargin_low);* idem lower and r0;
* variance of difference in sequence 1 (test -> reference) with n1 subjects;
se2= 2*sigma_log**2/n1; *variance under test and under ref period the same;
*** power;
power=probt( delta_up/sqrt(se2) - tinv(1-alpha2/2, n1-1),n1-1 )
      + probt( delta_low/sqrt(se2) - tinv(1-alpha2/2,n1-1),n1-1  )
      -1;
run;

proc print;run;

