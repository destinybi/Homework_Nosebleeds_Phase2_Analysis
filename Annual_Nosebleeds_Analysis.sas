****************************************************************************;
****************************************************************************;
* author: Bi Bozhu*;
* date: 2019-09-30*;
* purpose: Bayer nosebleed data anlysis;
* version: V1;
****************************************************************************;

****************************************************************************;

**********************  Clear previous SAS log and list **********************;
%macro ClearLogList();
DM 'log; "clear";';
DM 'Output; "clear";';
run;
%mend;
%ClearLogList();

%let folder=C:\Projects\Nosebleed_case\Scenario;

libname data "&folder.\data";

%include "&folder.\Vuong_test.sas";

*******************************************************************************;
*import all the data files and integrated them into one dataset;
*******************************************************************************;

%macro input_data(table);

proc import datafile="&folder.\data\&table..csv" out=&table. dbms=csv replace;
run;

proc sort data=&table.;
	by subject;
run;

%mend;

%input_data(efficacy);
%input_data(subject);
%input_data(randomization);

data baseline_efficacy_info;
	merge subject randomization efficacy;
	by subject;
	log_mucus_viscosity=log(mucus_viscosity+1);   *log transformation of mucus viscosity to be close to normal distribution ;
	if nosebleeds<previous_year then response=1;  *define the response ;
		else response=0;              
	if nosebleeds=0 then duration_year=0;   
		else duration_year=log(nosebleeds/duration*365); *to estimate the nosebleeds for patients with <365 days duration;
run;

*******************************************************************************;
*Exploratory data analysis to profile the data and take actions if needed;
*******************************************************************************;


*profile the continuous variables by arms ;
ods graphics on;

proc means data = baseline_efficacy_info n mean var min max; 
	class arm;
	var mucus_viscosity ;
run;

proc sort data=baseline_efficacy_info;
	by arm;
run;

proc boxplot data=baseline_efficacy_info;  *create box plot ;
	plot mucus_viscosity*arm;
run;

*profile the categorical variables by arms, and test the difference of arms at baseline ;

proc freq data=baseline_efficacy_info;
	table arm*country arm*eye_colour arm*tissue_use arm*previous_year/plots=all chisq norow ; *output the plots and chisq value;
run;

*profile the quantitive/continous variables by arms;

proc ttest data=baseline_efficacy_info;
	class arm;
	var mucus_viscosity log_mucus_viscosity nosebleeds;
run;

proc ttest data=baseline_efficacy_info;
	where nosebleeds=0;
	class arm;
	var mucus_viscosity log_mucus_viscosity;
run;


*******************************************************************************;
*Poisson model to quantify the impact of predictors; 
*******************************************************************************;

%macro fit_model(dist);

proc genmod data=baseline_efficacy_info;
	class arm tissue_use country;  *define categorical variables;
	model nosebleeds=arm tissue_use mucus_viscosity  previous_year country /dist=&dist. offset=duration_year type1 ; *poission or negative binomial regression ;
	zeromodel arm tissue_use mucus_viscosity  previous_year country; *this is for zero inflated model to fit the logistic model, estimating the excess zeros incidence ;
	output out=out_&dist. pred=pred_&dist. pzero=p0;
run;

%mend;

%fit_model(poisson);
%fit_model(zip); *zip=zero inflated possion;
%fit_model(negbin); *zip=negative binomial regression;
%fit_model(zinb);*zip=zero inflated nagative binomial regression ;

*test the difference at baseline for subgroup and poisson regression for subgroup with mucus_viscosity>1;
*firstly, test if the baseline info are aligned for subgroups, and subsquently, fit the poisson regression;

proc freq data=baseline_efficacy_info;
	where mucus_viscosity>=1; *median is 1 ;
	table arm*country arm*eye_colour arm*tissue_use arm*previous_year /plots=all chisq norow ; *output the plots and chisq value;
run;

proc ttest data=baseline_efficacy_info;
	where mucus_viscosity>=1;
	class arm;
	var mucus_viscosity log_mucus_viscosity nosebleeds;
run;

proc genmod data=baseline_efficacy_info;
	where mucus_viscosity>=1;
	class arm tissue_use country;
	model nosebleeds=arm tissue_use  mucus_viscosity previous_year country /dist=poisson offset=duration_year type1; *offset is to adjust the non-365 days patients ;
	output out=out_poi pred=pred_poi;
run;

*test the difference at baseline for subgroup and poisson regression for subgroup with mucus_viscosity>1 and tissue_use=HIGH;

proc freq data=baseline_efficacy_info;
	where mucus_viscosity>=1 and tissue_use="HIGH"; *median is 1 ;
	table arm*country arm*eye_colour arm*tissue_use arm*previous_year /plots=all chisq norow ; *output the plots and chisq value;
run;

proc ttest data=baseline_efficacy_info;
	where mucus_viscosity>=1 and tissue_use="HIGH"; 
	class arm;
	var mucus_viscosity log_mucus_viscosity nosebleeds;
run;

proc genmod data=baseline_efficacy_info;
	where mucus_viscosity>=1 and tissue_use="HIGH";
	class arm tissue_use country;
	model nosebleeds=arm tissue_use  mucus_viscosity previous_year country /dist=poisson offset=duration_year type1 ; 
	output out=out_poi pred=pred_poi;
run;

ods graphics off;

/* this is to for Vuong test to compare the models  */

/*proc genmod data=baseline_efficacy_info;*/
/*	class arm tissue_use country;*/
/*	model nosebleeds=arm tissue_use mucus_viscosity  previous_year country /dist=poisson offset=duration_year type1 SCALE=DEVIANCE; *zip=zero inflated possion;*/
/*	output out=out_poi pred=pred_poi;*/
/*run;*/
/**/
/*proc genmod data=out_poi;*/
/*	class arm tissue_use country;*/
/*	model nosebleeds=arm tissue_use mucus_viscosity  previous_year country /dist=zip offset=duration_year type1 SCALE=DEVIANCE; *zip=zero inflated possion;*/
/*	zeromodel  ;*/
/*	output out=out pred=pred_zip pzero=p0;*/
/*run;*/
/**/
/**/
/*%vuong(data=out, response=nosebleeds,*/
/*       model1=poi, p1=pred_poi, dist1=poi, scale1=1.00, pzero2=p0, */
/*       model2=zip, p2=pred_zip, dist2=zip, scale2=1.00,*/
/*       nparm1=5,   nparm2=5);*/

/*****calculate the pvalue;*/
/*data pvalue;*/
/*  df = 429; chisq = 321.87 ;*/
/*  pvalue = 1 - probchi(chisq, df);*/
/*run;*/
/*proc print data = pvalue noobs;*/
/*run;*/

****************************** Auto-Save the Log and List ************************************;
%macro AutoSave();
DM 'LOG; log; FILE "&folder.\SAS Logs and Lists\nosebleed_analysis.log" replace;';  
DM 'Output; Output; FILE "&folder.\SAS Logs and Lists\nosebleed_analysis.lst" replace;';  
run;
%mend;
%AutoSave();
