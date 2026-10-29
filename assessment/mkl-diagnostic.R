# Diagnostic only: unchanged Brown rt_modality_full fixture and tolerances.
# Usage: Rscript mkl-diagnostic.R [repository-root] [output-directory]
# Optional: MIXEFF_DIAGNOSTIC_LIBRARY=/path/to/library
#           MIXEFF_DIAGNOSTIC_TIGHT_ENGINE=true (one additional bounded fit)
args <- commandArgs(TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath('.')
outdir <- if (length(args) > 1L) args[[2]] else file.path(tempdir(), 'brown-mkl-diagnostic')
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
lib <- Sys.getenv('MIXEFF_DIAGNOSTIC_LIBRARY')
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
library(mixeff)
stopifnot(requireNamespace('lme4', quietly = TRUE))
options(digits = 16)
cat('UTC:', format(Sys.time(), tz='UTC'), '\nROOT:',root,'\nOUTPUT:',outdir,'\n')
print(sessionInfo()); print(extSoftVersion())
print(Sys.getenv(c('MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','RAYON_NUM_THREADS')))
for (p in c('mixeff','lme4','Matrix','minqa','nloptr','testthat')) {
 if (requireNamespace(p,quietly=TRUE)) cat(p,as.character(packageVersion(p)),find.package(p),'\n')
}
lock <- file.path(root,'src/rust/upstream/mixeff-rs.lock')
if (file.exists(lock)) cat(readLines(lock),sep='\n')
fixture <- file.path(root,'tests/fixtures/brown_rt_dummy_data.csv.gz')
cat('\nFixture MD5:',tools::md5sum(fixture),'\n')
d <- read.csv(fixture,stringsAsFactors=FALSE)
d$PID<-factor(d$PID); d$stim<-factor(d$stim)
d$modality<-ifelse(d$modality=='Audio-only',0,1)
f <- RT ~ 1 + modality + (1 + modality | PID) + (1 + modality | stim)
stopifnot(nrow(d)==21679L,nlevels(d$PID)==53L,nlevels(d$stim)==543L,
          identical(sort(unique(d$modality)),c(0,1)))
cat('FORMULA:',deparse(f),'REML=TRUE\n')
cat('PID/modality cell counts:\n'); print(summary(as.numeric(table(d$PID,d$modality))))
cat('stim/modality cell counts:\n'); print(summary(as.numeric(table(d$stim,d$modality))))
fits <- list(); records <- list()
run <- function(label, expr) {
 cat('\nBEGIN',label,'\n'); flush.console()
 warnings <- character(); started<-proc.time()[['elapsed']]
 value <- tryCatch(withCallingHandlers(force(expr),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 elapsed<-proc.time()[['elapsed']]-started
 cat('END',label,'elapsed_seconds=',elapsed,'\n')
 if(length(warnings)) print(warnings)
 if(inherits(value,'condition')) {cat('ERROR:',conditionMessage(value),'\n');return(NULL)}
 vc<-as.data.frame(if(inherits(value,'mm_fit')) mixeff::VarCorr(value) else lme4::VarCorr(value))
 beta<-if(inherits(value,'mm_fit')) mixeff::fixef(value) else lme4::fixef(value)
 cat('logLik:',as.numeric(logLik(value)),' REMLcriterion:',-2*as.numeric(logLik(value)),' sigma:',sigma(value),'\n')
 print(beta); print(vc)
 if(inherits(value,'mm_fit')) {print(value$fit$optimizer);print(mixeff::optimizer_certificate(value));print(value$artifact$optimizer_certificate)} else print(value@optinfo)
 fits[[label]] <<- value
 records[[label]] <<- list(elapsed=elapsed,warnings=warnings,logLik=as.numeric(logLik(value)),sigma=sigma(value),fixef=beta,varcorr=vc)
 saveRDS(value,file.path(outdir,paste0(label,'.rds')))
 value
}
mm <- run('mixeff_default',mixeff::lmm(f,d,REML=TRUE,control=mixeff::mm_control(verbose=-1)))
ref <- run('lme4_default',lme4::lmer(f,d,REML=TRUE,control=lme4::lmerControl(optimizer='bobyqa')))
stopifnot(!is.null(mm),!is.null(ref))
# Group/name joins reproduce the tested estimand (SD, not variance), with
# uniqueness checks that the original helper does not currently enforce.
vm<-as.data.frame(mixeff::VarCorr(mm)); vr<-as.data.frame(lme4::VarCorr(ref))
vm<-vm[is.na(vm$var2),c('grp','var1','sdcor')];vr<-vr[is.na(vr$var2),c('grp','var1','sdcor')]
key<-function(v)paste(v$grp,v$var1,sep='/')
stopifnot(!anyDuplicated(key(vm)),!anyDuplicated(key(vr)),setequal(key(vm),key(vr)))
vm<-vm[match(key(vr),key(vm)),]
cat('\nEXACT DIAGONAL SD ALIGNMENT (helper tolerance 0.3):\n')
print(data.frame(component=key(vr),mixeff_SD=vm$sdcor,lme4_SD=vr$sdcor,
                 absolute_gap=abs(vm$sdcor-vr$sdcor),relative_gap=abs(vm$sdcor-vr$sdcor)/abs(vr$sdcor)))
cat('Fitted max absolute gap:',max(abs(fitted(mm)-fitted(ref))),
    'response scaled:',max(abs(fitted(mm)-fitted(ref)))/max(abs(fitted(ref))),'\n')
# Reconstruct lme4 theta from ABSOLUTE covariance summaries and residual SD.
# Validate the converter on the reference before evaluating engine theta.
cnms<-lme4::getME(ref,'cnms')
theta_from_vc<-function(object) {
 v<-as.data.frame(if(inherits(object,'mm_fit'))mixeff::VarCorr(object) else lme4::VarCorr(object))
 result<-numeric()
 for(g in names(cnms)) {
  nm<-cnms[[g]]; stopifnot(length(nm)==2L)
  S<-matrix(0,2,2,dimnames=list(nm,nm))
  for(j in 1:2) {r<-which(v$grp==g & is.na(v$var2) & v$var1==nm[j]);stopifnot(length(r)==1L);S[j,j]<-v$vcov[r]}
  r<-which(v$grp==g & !is.na(v$var2) & ((v$var1==nm[1]&v$var2==nm[2])|(v$var1==nm[2]&v$var2==nm[1])))
  stopifnot(length(r)==1L);S[1,2]<-S[2,1]<-v$vcov[r]
  S<-S/sigma(object)^2
  a<-sqrt(S[1,1]);b<-if(a>0)S[2,1]/a else 0
  if(a==0)stopifnot(abs(S[2,1])<1e-12)
  remaining<-S[2,2]-b*b;stopifnot(remaining>= -1e-10*max(1,max(abs(S))))
  result<-c(result,a,b,sqrt(max(0,remaining)))
 }
 names(result)<-names(lme4::getME(ref,'theta'));result
}
theta_ref<-theta_from_vc(ref);theta_mm<-theta_from_vc(mm)
stopifnot(max(abs(theta_ref-lme4::getME(ref,'theta')))<1e-8)
cat('\nTHETA IN VERIFIED LME4 GROUP/BASIS ORDER:\n');print(cbind(lme4=theta_ref,mixeff=theta_mm))
devfun<-lme4::lmer(f,d,REML=TRUE,devFunOnly=TRUE,control=lme4::lmerControl(optimizer='bobyqa'))
cat('COMMON lme4 REML OBJECTIVE: reference=',devfun(theta_ref),'engine=',devfun(theta_mm),'\n')
cat('ENGINE objective discrepancy vs common evaluator=',devfun(theta_mm)+2*as.numeric(logLik(mm)),'\n')
# BOBYQA stopping controls from lme4's public interface, with explicit cap.
tight<-lme4::lmerControl(optimizer='bobyqa',optCtrl=list(rhobeg=0.01,rhoend=1e-9,maxfun=20000))
r1<-run('lme4_tight_from_reference',lme4::lmer(f,d,REML=TRUE,start=list(theta=theta_ref),control=tight))
r2<-run('lme4_tight_from_mixeff',lme4::lmer(f,d,REML=TRUE,start=list(theta=theta_mm),control=tight))
if(identical(Sys.getenv('MIXEFF_DIAGNOSTIC_TIGHT_ENGINE'),'true')) {
 run('mixeff_tight',mixeff::lmm(f,d,REML=TRUE,control=mixeff::mm_control(verbose=-1,ftol_rel=1e-12,ftol_abs=1e-8,max_feval=5000)))
}
if(!is.null(r1)&&!is.null(r2)) {
 cat('\nCOMMON OBJECTIVES AT ALL RECORDED COVARIANCES:\n')
 objectives<-vapply(fits,function(z)devfun(theta_from_vc(z)),numeric(1));print(objectives)
 cat('Objective excess above best:',objectives-min(objectives),'\n')
 best<-fits[[which.min(objectives)]];tb<-theta_from_vc(best)
 grid<-seq(0,1,length.out=11)
 cat('LINE BETWEEN DEFAULT MIXEFF AND BEST OBSERVED THETA (not a profile likelihood):\n')
 print(data.frame(fraction=grid,objective=vapply(grid,function(a)devfun((1-a)*theta_mm+a*tb),numeric(1))))
 cat('COMMON OBJECTIVE CENTRAL DIFFERENCE GRADIENTS (h=1e-5*max(1,abs(theta_j))):\n')
 gradient<-function(theta) vapply(seq_along(theta),function(j){h<-1e-5*max(1,abs(theta[j]));a<-b<-theta;a[j]<-a[j]+h;b[j]<-b[j]-h;if(b[j]<0&&lme4::getME(ref,'lower')[j]==0)return(NA_real_);(devfun(a)-devfun(b))/(2*h)},numeric(1))
 print(cbind(mixeff=gradient(theta_mm),reference=gradient(theta_ref),best=gradient(tb)))
}
saveRDS(records,file.path(outdir,'records.rds'))
cat('\nDONE: diagnostics only; no test tolerance changed.\n')
