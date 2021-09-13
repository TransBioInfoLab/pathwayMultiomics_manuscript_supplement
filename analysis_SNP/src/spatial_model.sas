/****************
 spatial model
****************/
%macro modc_sp_summary (dat, input, out_genes, outstatus);

data gene_snp; set &input;  
keep chromosome genename snpid location; run; %sort (gene_snp nodup, genename snpid);

* make indicator vars for each gene; 
data c; set &input; keep snpid genename ; run;
%sort (c nodup, snpid genename); 

data d; set c; large_gene=1; run; %sort(d, snpid); 
proc transpose data=d out=e(drop=_name_) ; by snpid; id genename; var large_gene; run;

proc contents data=e out=out noprint; run;

%sort (out, name); 
data _null_; set out(where=(name ne "snpid")) end=eof;
call symput ("genename"||compress(_n_), compress(name));
if eof then call symput ('ngenes',compress(_n_)); run;

* take out duplicates-snps on multiple genes; 
data noduplicates; set &input; drop genename gene_nsnps; run;
proc sort data=noduplicates nodup; by snpid; run; 

* calculate gene indicator vars; 
%sort (e, snpid); %sort(noduplicates, snpid); 
data f; merge e noduplicates; by snpid; 
%do i=1 %to &ngenes; if &&genename&i = . then &&genename&i=0; %end;
locs = location/1e7;
run;

%sort (f, chromosome snpid); 

* compute design matrix for each gene; 
data one; set f; keep snpid chromosome; run;
%sort (one nodup, snpid); %sort (gene_snp, snpid); 

data two; merge one gene_snp; by snpid; ind=1; keep genename chromosome ind; run;
%sort (two nodup, genename); 

proc glmmod data=two outparm=parm outdesign=design; 
class chromosome; model ind = chromosome/noint; run;

data _null_; set parm nobs=nk; call symput ('nchrom',nk); run; 

data design2; set design; desg = catx(" ", of col1-col%eval(1*&nchrom)); run;

data _null_; set design2 end=eof; call symput ("gene"||compress(_n_), desg); run;

* fit mixed model; 
proc glimmix data=f maxopt=100; 
class chromosome; 
model chisq = /s dist=gamma link=log; 
random %do g=1 %to &ngenes; &&genename&g %end; 
/type=sp(exp) (locs) gcoord=mean subject=chromosome;

*overall; 
estimate "all genes " int 1/upper; 

* first gene; 
estimate "&genename1" |&genename1 1  /subject &gene1 ; 

%do k=2 %to %eval(&ngenes); 
estimate "&&genename&k" | &&genename&k 1/ subject &&gene&k;
%end; 

nloptions tech=nrridg; 
ods output ConvergenceStatus=status_mod2d; 
ods output estimates=est_mod2d; 
run;

data _null_; set status_mod2d; call symput ("status", compress(status)); run;

data &outstatus; format datset $50.; format Reason $51.;  set status_mod2d;
datset = "&dat"; mod="spatial"; run; 

%if &status = 0 %then %do; *converged;
	data &out_genes; format datset $50.; set est_mod2d;  mod="spatial";  
	if Label = "all genes"; 
	datset="&dat"; run; 
%end;

%else %do; *did not converge; 
	data &out_genes; format datset $50.; datset = "&dat"; run; 
%end; 

title;  
%mend; 
