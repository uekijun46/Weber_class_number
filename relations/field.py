# Common definitions for the fields k_{p^n} (p odd): Gaussian periods, small elements, norms (fixed point, NOT rigorous).
# Elements are written alpha = c + sum_i a_i eta_i, where eta_i = sum_{h in H} zeta^(g^i h), i mod N = p^n,
# zeta = exp(2 pi i / p^(n+1)), H the subgroup of order p-1 of (Z/p^(n+1))^*, g a generator.  Note sum_i eta_i = 0.
import random
class Field:
    def __init__(self,p,n):
        self.p,self.n=p,n; self.N=N=p**n; self.M=M=p**(n+1); phi=(p-1)*N
        prs=set([p]+[d for d in range(2,p) if (p-1)%d==0 and all(d%e for e in range(2,int(d**.5)+1))])
        self.g=next(x for x in range(2,M) if x%p and all(pow(x,phi//r,M)!=1 for r in prs))
        self.H=[pow(self.g,N*t,M) for t in range(p-1)]
        self.exps=[[pow(self.g,i,M)*h%M for h in self.H] for i in range(N)]
        self._E=None
    def split(self,q): return q!=self.p and pow(q,self.p-1,self.M)==1
    def fixed(self,S=176):
        if self._E is None:
            from mpmath import mp, mpf, cos, pi
            mp.prec=S+80
            self._E=[int(sum(cos(2*pi*e/self.M) for e in ex)*mpf(2)**S) for ex in self.exps]; self._S=S
        return self._E,self._S
    def norm_fixed(self,c,coef):                      # coef: list of (i, a_i); returns |N(alpha)| by rounding (not rigorous)
        E,S=self.fixed(); N=self.N; prod=1
        for j in range(N):
            v=c<<S
            for i,a in coef: v+=a*E[(i+j)%N]
            prod*=v
        sh=S*N; prod=abs(prod); return (prod+(1<<(sh-1)))>>sh
    def periods_mod(self,q,rnd=random):               # eta_i modulo one prime above the split prime q (search only)
        M,p=self.M,self.p; f=1
        while pow(q,f,M)!=1: f+=1
        if f==1:
            while True:
                z=pow(rnd.randrange(2,q),(q-1)//M,q)
                if pow(z,M//p,q)!=1: break
            return [sum(pow(z,e,q) for e in ex)%q for ex in self.exps]
        import flint
        ctx=flint.fq_default_ctx(q,f)
        while True:
            z=ctx.random_element()**((q**f-1)//M)
            if z**(M//p)!=ctx.one(): break
        zp=[ctx.one()]
        for _ in range(M-1): zp.append(zp[-1]*z)
        R=[]
        for ex in self.exps:
            s=ctx.zero()
            for e in ex: s+=zp[e]
            R.append(int(s.to_list()[0]) if hasattr(s,'to_list') else int(str(s)))
        return R
