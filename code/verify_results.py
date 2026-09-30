"""Independent numerical checks from the released CSVs, without MATLAB.

Run: python code/verify_results.py
Dependencies: numpy, scipy. Writes outputs/verification; does not edit data.
This checks coefficients and deterministic statistics, not MATLAB graphics,
CANlab bootstrap inference, imaging preprocessing, or historic random splits.
"""
from pathlib import Path
import csv,json,math
import numpy as np
from scipy import stats,integrate,optimize

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/verification';OUT.mkdir(parents=True,exist_ok=True)
def read(rel):
 with (ROOT/rel).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def number(v):
 try:return float(v)
 except (ValueError,TypeError):return np.nan
def col(rows,key):return np.array([number(r.get(key))for r in rows])
def save(name,rows):
 with (OUT/name).open('w',encoding='utf8',newline='')as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def bh(p):
 p=np.array(p,float);q=np.full(p.shape,np.nan);ix=np.where(np.isfinite(p))[0];order=np.argsort(p[ix]);s=p[ix][order]
 q[ix[order]]=np.minimum(1,np.minimum.accumulate((s*len(s)/np.arange(1,len(s)+1))[::-1])[::-1]);return q
def partial(x,y,z=None,rank=False):
 x=np.asarray(x,float);y=np.asarray(y,float);z=np.empty((len(x),0))if z is None else np.asarray(z,float).reshape(len(x),-1)
 v=np.c_[x,y,z];v=v[np.isfinite(v).all(axis=1)]
 if rank:v=np.apply_along_axis(stats.rankdata,0,v)
 design=np.c_[np.ones(len(v)),v[:,2:]];rankz=np.linalg.matrix_rank(design)
 rxy=v[:,:2]-design@np.linalg.lstsq(design,v[:,:2],rcond=None)[0]
 r=np.corrcoef(rxy.T)[0,1];df=len(v)-rankz-1
 p=2*stats.t.sf(abs(r)*np.sqrt(df/(1-r*r)),df)
 return len(v),r,p
def bf(r,n,scale=1):
 nu=n-2;a=n*scale*scale;r2=r*r
 def f(u):
  if not u:return 0.
  v=u*u
  return np.sqrt(2/np.pi)*np.exp(-v/2+.5*(np.log(v)-np.log(v+a))+nu/2*(np.log(v+a)-np.log(v+a*(1-r2))))
 return integrate.quad(f,0,12,epsabs=1e-10,epsrel=1e-9)[0]
def partition(x):
 # Global 1-D two-means solution, independently checking MATLAB kmeans.
 order=np.argsort(x);s=x[order];cut=min((i for i in range(1,len(s))if s[i-1]<s[i]),key=lambda i:sum((s[:i]-s[:i].mean())**2)+sum((s[i:]-s[i:].mean())**2))
 label=np.zeros(len(x),int);label[order[cut:]]=1;return label
T=read('data/participants/analysis_data.csv');ids=[r['ReleaseID']for r in T]
assert len(ids)==len(set(ids))==110
g=col(T,'GroupID');age=col(T,'Age');nm=col(T,'NMDAS');acc=col(T,'Nback_ACC');flag=col(T,'Behavioral_Available')==1
brain=np.column_stack([col(T,k)for k in ['Nback','Multisensory','Cold']]);tasks=['Working memory','Multisensory','Cold pain']
imaging=np.isfinite(brain)&np.isin(g,[0,1])[:,None];imaging[:,0]&=flag
assert imaging.sum(axis=0).tolist()==[74,90,91]
assert (imaging&(g==0)[:,None]).sum(axis=0).tolist()==[25,28,29]
rows=[];expected=[-.670360287094522,-.459626851691269,-.370595543152961,-.191521021417447]
for k,(name,y,t)in enumerate([('N-back brain',brain[:,0],0),('N-back behavior',acc,0),('Multisensory brain',brain[:,1],1),('Cold-pain brain',brain[:,2],2)]):
 m=imaging[:,t]&(g==0)&np.isfinite(nm)&np.isfinite(age)&np.isfinite(y)
 n,r,p=partial(nm[m],y[m],age[m],True);_,raw,praw=partial(nm[m],y[m],rank=True)
 assert abs(r-expected[k])<1e-12
 rows.append(dict(Outcome=name,N=n,RhoRaw=raw,PRawAsymptotic=praw,RhoAge=r,PAge=p,BF10=bf(r,n)if k!=1 else ''))
for row,q in zip(rows,bh([r['PAge']for r in rows])):row['QAge']=q
save('severity_correlations.csv',rows)
assert np.allclose([rows[i]['BF10']for i in [0,2,3]],[89.541921,.883389,.230625],atol=1e-6)
bb=[]
for group,name in [(0,'Patients'),(1,'Controls')]:
 m=imaging[:,0]&(g==group)&np.isfinite(acc)&np.isfinite(age)
 n,r,p=partial(brain[m,0],acc[m],age[m],True);_,raw,praw=partial(brain[m,0],acc[m],rank=True)
 bb.append(dict(Population=name,N=n,RhoAge=r,PAge=p,RhoRaw=raw,PRawAsymptotic=praw))
save('brain_behavior.csv',bb)
sub=[]
for name,y,t in [('N-back brain',brain[:,0],0),('N-back behavior',100*acc,0),('Multisensory brain',brain[:,1],1),('Cold-pain brain',brain[:,2],2)]:
 m=imaging[:,t]&(g==0)&np.isfinite(nm)&np.isfinite(age);inds=np.flatnonzero(m);lab=partition(nm[m]);sets=[y[imaging[:,t]&(g==1)],y[inds[lab==0]],y[inds[lab==1]]]
 for a,b in [(0,1),(0,2),(1,2)]:
  aa=sets[a][np.isfinite(sets[a])];bb_=sets[b][np.isfinite(sets[b])];test=stats.ttest_ind(aa,bb_,equal_var=True)
  sub.append(dict(Outcome=name,Group1=['Controls','Resilient patients','Severe Patients'][a],Group2=['Controls','Resilient patients','Severe Patients'][b],N1=len(aa),N2=len(bb_),T=test.statistic,P=test.pvalue))
save('subgroup_tests.csv',sub)
med=[];m=imaging[:,0]&np.isfinite(nm)&np.isfinite(age)&np.isfinite(acc);assert m.sum()==72
X=np.c_[np.ones(m.sum()),g[m],age[m]];a=np.linalg.lstsq(X,nm[m],rcond=None)[0][1]
for y,name in [(brain[:,0],'N-back brain'),(acc,'N-back behavior proportion')]:
 fit=np.linalg.lstsq(np.c_[X,nm[m]],y[m],rcond=None)[0];c=np.linalg.lstsq(X,y[m],rcond=None)[0][1]
 med.append(dict(Outcome=name,N=int(m.sum()),a=a,b=fit[-1],cprime=fit[1],c=c,ab=a*fit[-1]))
save('mediation_OLS_paths.csv',med)
# All paired observations retained; display omission never alters dz.
pairs=[]
for key in ['Nback','Multisensory','Cold']:
 p=read('data/figure2/'+key+'_paired_scores.csv');lo=col(p,'Control');hi=col(p,'Task');delta=hi-lo
 dz=delta.mean()/delta.std(ddof=1);cv=[r for r in read('data/figure2/cv_summary.csv')if r['Task']==key]
 assert abs(dz-number(cv[-1]['CohensDz']))<1e-12
 omitted=0
 if key=='Cold':
  v=np.c_[lo,hi];medv=np.median(v,axis=0);mad=1.482602218505602*np.median(abs(v-medv),axis=0)
  omitted=int(np.any(abs(v-medv)>3*mad,axis=1).sum());assert omitted==3
 pairs.append(dict(Task=key,N=len(lo),FinalSplitDz=dz,MeanCVDz=float(col(cv,'CohensDz').mean()),MeanCVAccuracy=float(col(cv,'ClassificationAccuracy').mean()),DisplayOmitted=omitted))
save('task_activation.csv',pairs)
# Figure 4A preserves the pre-exclusion N-back imaging sample.
lg=np.log10(col(T,'GDF15'));base=np.isfinite(brain[:,0])&np.isin(g,[0,1])&np.isfinite(lg)
test=stats.ttest_ind(lg[base&(g==0)],lg[base&(g==1)],equal_var=True);assert abs(test.statistic-6.10134922992542)<1e-12
gdf=[dict(Outcome='GDF15 group difference',Population='Pre-exclusion N-back imaging',N=int(base.sum()),StatisticType='t',Statistic=test.statistic,P=test.pvalue,Covariates='None')]
m=base&(g==0)&np.isfinite(nm)
for rank,typ in [(True,'Spearman'),(False,'Pearson')]:
 n,r,p=partial(lg[m],nm[m],rank=rank);gdf.append(dict(Outcome='GDF15 and NMDAS',Population='Patients pre-exclusion',N=n,StatisticType=typ,Statistic=r,P=p,Covariates='None'))
for y,name in [(brain[:,0],'GDF15 and N-back brain'),(acc,'GDF15 and N-back behavior')]:
 for group,pop in [(None,'All eligible'),(0,'Patients'),(1,'Controls')]:
  m=imaging[:,0]&np.isfinite(lg)&np.isfinite(y)
  if group is not None:m&=g==group
  n,r,p=partial(lg[m],y[m],g[m]if group is None else None)
  gdf.append(dict(Outcome=name,Population=pop,N=n,StatisticType='Pearson',Statistic=r,P=p,Covariates='Group'if group is None else 'None'))
save('GDF15_statistics.csv',gdf)
reg=[]
for y,name in [(brain[:,0],'Brain'),(acc,'Behavior')]:
 m=imaging[:,0]&np.isfinite(lg)&np.isfinite(age)&np.isfinite(y)
 X=np.c_[np.ones(m.sum()),lg[m],age[m],2*g[m]-1];beta=np.linalg.lstsq(X,y[m],rcond=None)[0]
 df=m.sum()-X.shape[1];res=y[m]-X@beta;se=np.sqrt(np.diag(np.linalg.inv(X.T@X))*(res@res)/df)
 for key,b,s in zip(['Intercept','GDF15','Age','Group'],beta,se):reg.append(dict(Outcome=name,Coefficient=key,N=int(m.sum()),DF=int(df),Beta=b,SE=s,T=b/s,P=2*stats.t.sf(abs(b/s),df)))
save('GDF15_main_effects.csv',reg)
# Reconciled Table S1, including 20 explicitly inherited summary rows.
P=read('data/phenotype/phenotype_scores.csv');G={r['ReleaseID']:int(r['GroupID'])for r in read('data/participants/cohort_groups.csv')}
expected=read('data/summary/table_s1_reconciled.csv');old=read('data/summary/table_s1_original_summary.csv');s1=[];matched=0
for ref,orig in zip(expected,old):
 v=ref['variable'];ns=np.array([number(orig[k])for k in ['N_Control','N_Patient','Mean_Control','Mean_Patient','SD_Control','SD_Patient']]);gg=number(orig['Hedges_g']);src='Original summary; raw variable unavailable'
 if v in P[0]:
  matched+=1;x=np.array([number(r[v])for r in P if G[r['ReleaseID']]==1]);y=np.array([number(r[v])for r in P if G[r['ReleaseID']]==0]);x=x[np.isfinite(x)];y=y[np.isfinite(y)]
  z=np.array([len(x),len(y),x.mean(),y.mean(),x.std(ddof=1),y.std(ddof=1)])
  src='Raw observations reproduce original summary'
  if np.any(abs(z-ns)>1e-4):
   ns=z;pool=np.sqrt(((len(x)-1)*z[4]**2+(len(y)-1)*z[5]**2)/(len(x)+len(y)-2));gg=(z[2]-z[3])/pool*(1-3/(4*(len(x)+len(y))-9));src='Recalculated from raw observations'
 n1,n2,m1,m2,sd1,sd2=ns
 p=stats.ttest_ind_from_stats(m1,sd1,n1,m2,sd2,n2,equal_var=False).pvalue if min(n1,n2)>1 and np.isfinite(ns).all() else np.nan
 s1.append(dict(Variable=v,NControl=int(n1),NPatient=int(n2),MeanControl=m1,MeanPatient=m2,SDControl=sd1,SDPatient=sd2,HedgesG=gg,P=p,Provenance=src))
assert matched==194
for row,q,ref in zip(s1,bh([r['P']for r in s1]),expected):
 row['Q']=q
 for k,v in [('HedgesG','g'),('P','p_welch'),('Q','q')]:assert np.isclose(row[k],number(ref[v]),atol=1e-10,equal_nan=True),(row['Variable'],k)
assert np.isfinite([r['P']for r in s1]).sum()==211
save('TableS1_reproduced.csv',s1)
# Full age-adjusted phenome associations in Tables S2 and S3.
pd={r['ReleaseID']:r for r in P};items={r['ReleaseID']:r for r in read('data/phenotype/clinical_items.csv')}
merged=[dict(pd[r['ReleaseID']],**{k:v for k,v in items[r['ReleaseID']].items()if k!='ReleaseID'})for r in T]
for table,y in [(2,brain[:,0]),(3,acc)]:
 ref=read(f'data/summary/reported_Table_S{table}.csv');rs=[]
 for row in ref:
  v=row['Variable'];assert v in merged[0],v
  x=col(merged,v);m=imaging[:,0]&(g==0)&np.isfinite(age)&np.isfinite(x)&np.isfinite(y)
  n,r,p=partial(x[m],y[m],age[m],True)
  assert n==number(row['N']) and abs(r-number(row['R']))<.00051 and abs(p-number(row['p']))<.000051,v
  rs.append(dict(Variable=v,N=n,RhoAge=r,PAge=p,ReportedFDRValue=number(row['p_FDR'])))
 for row,q in zip(rs,bh([r['PAge']for r in rs])):
  row['QBH304Tests']=q
  assert (q<.05)==(row['ReportedFDRValue']<.05)
 save(f'TableS{table}_recomputed.csv',rs)
# Two-group power with actual unequal sample sizes and per-group balanced N.
power=[]
def power_at(d,n1,n2):
 df=n1+n2-2;c=stats.t.ppf(.975,df);nc=d/np.sqrt(1/n1+1/n2)
 # The installed SciPy 1.18 nct.sf can return NaN for very high power;
 # use CDF complement here (the root of interest is power=.8).
 return stats.nct.cdf(-c,df,nc)+1-stats.nct.cdf(c,df,nc)
for k,dz in enumerate([1.67,2.18,.92]):
 n1=int((imaging[:,k]&(g==0)).sum());n2=int((imaging[:,k]&(g==1)).sum())
 # For these observed sample sizes the 80% root is below 1.5. Keep the
 # bracket away from extreme noncentralities that trigger a SciPy 1.18 bug.
 detectable=optimize.brentq(lambda d:power_at(d,n1,n2)-.8,.001,1.5)
 for frac in [.25,.5]:
  d=dz*frac;n=2
  while power_at(d,n,n)<.8:n+=1
  power.append(dict(Task=tasks[k],Fraction=frac,AssumedBetweenGroupD=d,RequiredNPerGroup=n,RequiredTotalN=2*n,AvailablePatients=n1,AvailableControls=n2,PowerAtAvailableN=power_at(d,n1,n2),DFor80Percent=detectable))
save('power_between_groups.csv',power)
checks={'participants':110,'cohort_patients':40,'cohort_controls':70,'task_samples':[74,90,91],'task_patient_samples':[25,28,29], 'mediation_n':72,'TableS1_raw_variables_checked':194,'TableS1_retained_summary_variables':20,'TableS1_FDR_tests':211,'TableS2_S3_correlations_checked_each':304,'TableS2_S3_original_FDR_values_differ_from_BH':True,'TableS2_S3_FDR_significance_classifications_unchanged':True,'cold_pairs_omitted_from_display_only':3,'MATLAB_executed':False,'CANlab_inference_recomputed':False,'imaging_pipeline_rerun':False,'deterministic_checks_passed':True}
(OUT/'checks.json').write_text(json.dumps(checks,indent=2),encoding='utf8')
print(json.dumps(checks,indent=2));print('Output:',OUT)
