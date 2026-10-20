#!/usr/bin/env python3
"""MASC-style binomial random-intercept analysis of per-mouse cluster counts.
Plots order clusters from lowest to highest log2 odds ratio in each comparison.

Keeps original comparisons: KO at each age vs ALL WT mice, and pooled KO vs ALL WT.
Clusters 6, 14, 15, 16 are excluded. BH-FDR is applied across clusters
separately for each comparison. Requires pandas, numpy, scipy, matplotlib.

Model for each cluster: y_mouse ~ Binomial(n_mouse, p_mouse),
logit(p_mouse) = intercept + beta*KO + u_mouse, u_mouse ~ N(0, sigma^2).
Marginal likelihood integrates u_mouse using Gauss-Hermite quadrature.
Wald z tests and 95% CIs use the inverse observed Hessian of marginal log likelihood.
"""
from pathlib import Path
import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import expit, logsumexp, gammaln
# Set the full path to your CSV here (Windows example):
INPUT_CSV = r"C:\Users\ioaman\OneDrive - Karolinska Institutet\Skrivbordet\MASC\masc_counts_template.csv"
# Results will be saved in a folder next to the CSV.
OUTPUT_DIR = None

EXCLUDED = {6, 14, 15, 16}
N_QUAD = 30
NODES, WEIGHTS = np.polynomial.hermite.hermgauss(N_QUAD)
LOG_WEIGHTS = np.log(WEIGHTS) - 0.5 * np.log(np.pi)

def fit_cluster(sub):
    y = sub['n_cells_in_cluster'].to_numpy(dtype=float)
    n = sub['total_cells_in_scope'].to_numpy(dtype=float)
    x = (sub['Group'] == 'KO').to_numpy(dtype=float)
    # One aggregated binomial observation per mouse; animal-specific random intercept.
    def nll(theta):
        a, b, log_sd = theta
        sd = np.exp(log_sd)
        eta = a + b*x
        logits = eta[:, None] + np.sqrt(2.0)*sd*NODES[None, :]
        # Binomial combinatorial terms do not depend on parameters and are omitted.
        loglik = y[:, None]*(-np.logaddexp(0, -logits)) + (n-y)[:, None]*(-np.logaddexp(0, logits))
        return -np.sum(logsumexp(loglik + LOG_WEIGHTS[None, :], axis=1))
    p0 = np.clip((y+0.5)/(n+1), 1e-6, 1-1e-6)
    starts = [(float(np.log(p0[x==0].mean()/(1-p0[x==0].mean()))),
               float(np.log(p0[x==1].mean()/(1-p0[x==1].mean())) - np.log(p0[x==0].mean()/(1-p0[x==0].mean()))), np.log(s))
              for s in (0.05, 0.2, 0.6)]
    fits = [minimize(nll, s, method='L-BFGS-B', bounds=[(-15,15),(-15,15),(-9,3)],
                    options={'maxiter':2000,'ftol':1e-11}) for s in starts]
    fit = min((f for f in fits if np.isfinite(f.fun)), key=lambda f:f.fun)
    from statsmodels.tools.numdiff import approx_hess
    H = approx_hess(fit.x, nll, epsilon=1e-4)
    # Profile out nuisance intercept and random-effect variance via inverse Hessian.
    try:
        cov = np.linalg.inv(H)
        se = np.sqrt(cov[1,1]) if cov[1,1] > 0 else np.nan
    except np.linalg.LinAlgError:
        se = np.nan
    beta = float(fit.x[1]); z = beta/se if np.isfinite(se) and se>0 else np.nan
    from scipy.stats import norm
    p = float(2*norm.sf(abs(z))) if np.isfinite(z) else np.nan
    boundary = bool(fit.x[2] < -8.9 or fit.x[2] > 2.9)
    reliable = bool(fit.success and np.isfinite(se) and not boundary)
    return dict(log_odds=beta, se_log_odds=se, z_stat=z, p_value=p,
                random_intercept_sd=float(np.exp(fit.x[2])), converged=bool(fit.success),
                boundary_variance=boundary, reliable=reliable,
                n_animals_KO=int(x.sum()), n_animals_WT=int((1-x).sum()),
                optimizer_message=str(fit.message))

def analyze(data, label):
    rows=[]
    for cl, sub in data.groupby('Cluster'):
        if sub.Group.nunique()!=2: continue
        out=fit_cluster(sub)
        rows.append(dict(comparison=label, cluster=int(cl), **out))
    results=pd.DataFrame(rows).sort_values('cluster').reset_index(drop=True)
    from statsmodels.stats.multitest import multipletests
    # Unreliable fits must not be interpreted as statistically significant.
    mask=results.reliable & results.p_value.notna()
    results['q_value']=np.nan
    if mask.any():
        results.loc[mask,'q_value']=multipletests(results.loc[mask,'p_value'],method='fdr_bh')[1]
    scale=np.log(2)
    results['log2_or']=results.log_odds/scale
    results['lcl_log2_or']=(results.log_odds-1.96*results.se_log_odds)/scale
    results['ucl_log2_or']=(results.log_odds+1.96*results.se_log_odds)/scale
    results['odds_ratio']=np.exp(results.log_odds)
    return results

def plot(df, output, title):
    # Sort separately for each comparison, from lowest to highest log2 odds ratio.
    # This changes plot order only; statistical results remain unchanged.
    df = df.sort_values("log2_or", ascending=True, na_position="last").reset_index(drop=True)
    fig,ax=plt.subplots(figsize=(10,5))
    x=np.arange(len(df)); y=df.log2_or.to_numpy()
    lower=df.lcl_log2_or.to_numpy(); upper=df.ucl_log2_or.to_numpy()
    ax.errorbar(x,y,yerr=np.vstack([y-lower,upper-y]),fmt='o',capsize=3)
    ax.axhline(0,linestyle='--',color='gray')
    ax.set_xticks(x,df.cluster.astype(str)); ax.set_xlabel('Cluster')
    ax.set_ylabel('log₂ odds ratio (KO vs WT)');ax.set_title(title)
    fig.tight_layout();fig.savefig(output);plt.close(fig)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input', default=INPUT_CSV)
    parser.add_argument('--output', default=OUTPUT_DIR)
    args=parser.parse_args()
    data=pd.read_csv(args.input)
    required=['animal_id','Group','age','Cluster','n_cells_in_cluster','total_cells_in_scope']
    if any(c not in data.columns for c in required):
        raise ValueError('Missing columns: '+str(set(required)-set(data.columns)))
    data=data.loc[~data.Cluster.isin(EXCLUDED)].copy()
    if data[required].isna().any().any(): raise ValueError('Missing data in required columns')
    if data.duplicated(['animal_id','Cluster']).any(): raise ValueError('Duplicate mouse-cluster rows')
    if (data.n_cells_in_cluster>data.total_cells_in_scope).any(): raise ValueError('Cluster count exceeds total')
    if (data.total_cells_in_scope<=0).any(): raise ValueError('Nonpositive total')
    if not set(data.Group.unique())=={'WT','KO'}: raise ValueError('Group labels must be WT and KO')
    outdir=Path(args.output) if args.output else Path(args.input).expanduser().resolve().parent / 'MASC_mixed_results'
    outdir.mkdir(parents=True,exist_ok=True)
    frames=[]
    for age in sorted(data.loc[data.Group=='KO','age'].unique()):
        comp=data.loc[((data.Group=='KO')&(data.age==age))|(data.Group=='WT')]
        frames.append(analyze(comp,f'{age}_vs_all_WT'))
    per_age=pd.concat(frames,ignore_index=True)
    pooled=analyze(data,'AllAges_KO_vs_all_WT')
    per_age.to_csv(outdir/'per_age_masc_mixed.csv',index=False)
    pooled.to_csv(outdir/'pooled_masc_mixed.csv',index=False)
    for age,df in per_age.groupby('comparison'):
        plot(df,outdir/f'log2OR_CI_{age}.pdf',age)
    plot(pooled,outdir/'log2OR_CI_overall.pdf','Pooled KO vs all WT')
    for name,df in [('per_age',per_age),('pooled',pooled)]:
        print(name, 'fits:',len(df),'reliable:',int(df.reliable.sum()),'flagged:',int((~df.reliable).sum()))
        if (~df.reliable).any():
            print('CHECK flagged fits:',df.loc[~df.reliable,['comparison','cluster','random_intercept_sd','optimizer_message']].to_string(index=False))
    print('Results saved to',outdir.resolve())
    print('CAUTION: Small animal counts; Wald inference may be unstable. Review diagnostics before publication.')

if __name__=='__main__': main()
