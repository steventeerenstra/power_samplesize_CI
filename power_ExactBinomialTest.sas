** based on a macro by Ton de Haan;


********** one-sided testing on the lower limit of the confidence interval *****;
**>>>> sample size needed for a specific power *****************************;
%macro clopperpearson(e=,n=,alpha=0.05);
  /* calculate the clopper-pearson confidence interval */
  e=&e.;n=&n.;alpha=&alpha.;* alpha is two-sided;
  p=e/n;
  if e>0 then LO_CLP=betainv(alpha/2, e, n-e+1);
  else if e=0 then LO_CLP=0;
  if e<n then UP_CLP=betainv(1-alpha/2, e+1, n-e);
  else if e=n then UP_CLP=1;  

%mend clopperpearson;


%let alpha=0.05;*two-sided;
%let power=0.8;
%let prob0=0.25;   
%let lowerlim=0.12; 


data zoek;
  do ss=5 to 100 by 5;  /* first a quick scan of sample sizes */
    do event=0 to ss;
      %clopperpearson(e=event,n=ss,alpha=&alpha); 
      /* if lower limit confidence interval above criteria */
      if LO_CLP >= &lowerlim  then 
        do; *probability to have {event, event+1,...,ss} events given prob0 and sample size ss; 
			power=1-probbnml(&prob0,ss,event-1); 
          if power>&power then do;
            goto bijna  ;
          end;
        end;
      end;
    end;
    bijna:; 
 
  do sss=1 to ss;   /* precise scan of sample sizes */
    do event=0 to sss;
      %clopperpearson(e=event,n=sss,alpha=&alpha);    
      /* if lower limit confidence interval above criteria */
      if LO_CLP >= &lowerlim  then 
        do; *probability to have {event, event+1,...,ss} events given prob0 and sample size ss; 
			power=1-probbnml(&prob0,ss,event-1); 
          if power>&power then do;output;
            goto klaar  ;
          end;
        end;
      end;
    end; 
    klaar:;                  
run;  

proc print;run;

*****>>>>power given sample size ********************************************;

/*
data b;
retain power 0;* set power;
sample_size=72;
prob0=0.44;
lowerlimit=0.25;
do events=0 to sample_size;
	%clopperpearson(e=events,n=sample_size,alpha=0.05);
	prob_binom=pdf("binomial",events,prob0,sample_size);
	if LO_CLP >= lowerlimit then power = power + prob_binom;
	else power=power;
	output;
end;
run;

proc print data=b;run;
*/

data c;
alpha=0.05; * two-sided;
sample_size=80;
prob0=0.34;
lowerlimit=0.20;
do events=1 to sample_size; *note we start at as 0 will never give a CI above the lower limit;
	%clopperpearson(e=events,n=sample_size,alpha=alpha);
	power=1-probbnml(prob0,sample_size,events-1); 
	if LO_CLP >= lowerlimit then goto klaar;
end;
klaar: output;
run;

proc print data=c;run;

***************************************************************************************;
