import numpy as np
import pandas as pd

def parse_positions(m):
    return [int(t[1:-1]) for t in str(m).split(":") if len(t)>=3 and t[1:-1].isdigit()]
def parse_aa(m):
    return {int(t[1:-1]):t[-1] for t in str(m).split(":") if len(t)>=3 and t[1:-1].isdigit()}

def load(path, enc, max_n=1500, seed=0):
    df=pd.read_csv(path).dropna(subset=["phenotype"])
    if len(df)>max_n: df=df.sample(max_n,random_state=seed)
    y=df["phenotype"].values.astype(float)
    if enc=="mutation_list":
        aas=df["mutant"].apply(parse_aa)
        pos=sorted({p for d in aas for p in d})
        wt={p:"_" for p in pos}
        cols=[]; names=[]; groups=[]
        for p in pos:
            letters=sorted({d.get(p,wt[p]) for d in aas})
            if len(letters)<2: continue
            for L in letters[1:]:                       # drop one ref level per position
                cols.append([1.0 if d.get(p,wt[p])==L else 0.0 for d in aas])
                names.append(f"p{p}_{L}"); groups.append(f"p{p}")
        X=np.array(cols,float).T
    else:
        feats=[c for c in df.columns if c!="phenotype"]
        cols=[]; names=[]; groups=[]
        for c in feats:
            s=df[c]
            if pd.api.types.is_numeric_dtype(s):
                cols.append(s.values.astype(float)); names.append(c); groups.append(c)
            else:
                levels=sorted(s.astype(str).unique())
                for L in levels[1:]:
                    cols.append((s.astype(str)==L).values.astype(float))
                    names.append(f"{c}={L}"); groups.append(c)
        X=np.column_stack(cols)
    # drop zero-variance columns
    keep=[i for i in range(X.shape[1]) if np.std(X[:,i])>1e-12]
    X=X[:,keep]; names=[names[i] for i in keep]; groups=[groups[i] for i in keep]
    nominal_d=len(set(groups))
    return X,y,names,groups,nominal_d
