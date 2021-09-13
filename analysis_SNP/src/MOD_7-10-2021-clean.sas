
options mprint nodate nocenter formdlim=' ' NOSYNTAXCHECK; 

%let dir = C:\Users\lxw391\TBL Dropbox\Lily Wang\pathwayPCA_multiOmics\LWtest\GWAS_pathways; 

libname h "&dir\DATA"; libname r "&dir\RESULT"; 

%include "&dir\spatial_model.sas"; 

%macro sort (dat, by); proc sort data=&dat; by &by; run; %mend; 

%let ngs=2833; 

%macro q (results_dat);

data &results_dat; run;
                                     
%do gs=1 %to &ngs;  

PROC IMPORT OUT=geneset&gs DATAFILE= "C:\Users\lxw391\TBL Dropbox\Lily Wang\pathwayPCA_multiOmics\GWAS_pathway_method\gene_set_tables\c2cp\c2_cp_geneset_&gs..xlsx" REPLACE; GETNAMES=yes; RUN;

data temp; set geneset&gs; if _n_ ne 1; if chisq = 0 then chisq = 1.221234e-8; run; 

data temp2; set temp; if _n_ = 1; run;

data _null_; set temp2; 
call symput ("ngenes",  compress(ngenes)); 
call symput ("nsnps", compress(nsnps) ); 
call symput ("pathway", compress(pathway) ); run;

	title "Pathway = &pathway"; 

	%modc_sp_summary (dat=&pathway, input=temp, out_genes=one_mod2d, outstatus=one_status_mod2d);

	data one_mod; format datset $50.; format Reason $51.; merge one_mod2d one_status_mod2d; by datset; geneset = &gs; run;

	data &results_dat; format datset $50.; format Reason $51.; set &results_dat one_mod; run; 

proc datasets library=work kill; run; quit; 

%end; 
%mend; 

%q (r.results); 


