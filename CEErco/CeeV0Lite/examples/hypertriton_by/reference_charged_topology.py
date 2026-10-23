"""Uniform By charged-mother extension for a v0.6.2 toy, not the stock Finder.

Measured prior: x=(r_D,p_D,E,r_P), full 10x10 covariance.
Exact constraint: r_D-r_P=A(S)p_D. S=arc length/|p| is a free nuisance.
At each iteration eliminate S analytically in the linearized constraint,
update against the SAME prior, then form a full 8D covariance at decay.
No mass or lifetime constraint. The production-end transport includes d/dS.
"""
import numpy as np

KAPPA=0.00299792458  # cm, GeV/c, Tesla

def helix_matrices(s, by, charge=1):
    s=np.atleast_1d(s); k=KAPPA*by*charge; a=k*s
    sn=np.sin(a); co=np.cos(a)
    si=s*np.sinc(a/np.pi)
    om=s*(a/2)*np.sinc(a/(2*np.pi))**2
    A=np.zeros((len(s),3,3)); Ap=A.copy(); App=A.copy()
    A[:,0,0]=A[:,2,2]=si; A[:,1,1]=s
    A[:,0,2]=om; A[:,2,0]=-om
    Ap[:,0,0]=Ap[:,2,2]=co;Ap[:,1,1]=1
    Ap[:,0,2]=sn;Ap[:,2,0]=-sn
    App[:,0,0]=App[:,2,2]=-k*sn
    App[:,0,2]=k*co;App[:,2,0]=-k*co
    return A,Ap,App

def fit(q,c,pv,pvc,by=.5,maxiter=50):
    q=np.asarray(q); n=len(q); x0=np.column_stack((q,pv)); x=x0.copy()
    C=np.zeros((n,10,10));C[:,:7,:7]=c;C[:,7:,7:]=pvc
    s=np.sum((q[:,:3]-pv)*q[:,3:6],axis=1)/np.sum(q[:,3:6]**2,axis=1)
    converged=np.zeros(n,bool);bad=np.zeros(n,bool); nit=np.zeros(n,int)
    eye=np.eye(10)[None,:,:]
    for it in range(maxiter):
        A,Ap,_=helix_matrices(s,by)
        v=np.einsum('nij,nj->ni',Ap,x[:,3:6])
        g=x[:,:3]-x[:,7:]-np.einsum('nij,nj->ni',A,x[:,3:6])
        H=np.zeros((n,3,10));H[:,:,:3]=np.eye(3);H[:,:,3:6]=-A;H[:,:,7:]=-np.eye(3)
        CH=C@H.transpose(0,2,1);W=H@CH
        # Eigenvalues are checked, never replace a singular input by a fake covariance.
        ev=np.linalg.eigvalsh(W); singular=(ev[:,0]<=1e-14*np.maximum(ev[:,-1],1e-30))|~np.isfinite(ev).all(axis=1)
        bad|=singular; W[singular]=np.eye(3)
        V=np.linalg.inv(W);Vv=np.einsum('nij,nj->ni',V,v)
        den=np.sum(v*Vv,axis=1);bad|=(den<=0)|~np.isfinite(den);den=np.where(den>0,den,1.)
        Q=V-Vv[:,:,None]*Vv[:,None,:]/den[:,None,None]
        b=g+np.einsum('nij,nj->ni',H,x0-x)
        ds=np.sum(Vv*b,axis=1)/den
        xn=x0-np.einsum('nij,nj->ni',CH@Q,b);sn=s+ds
        change=np.max(np.abs(xn-x)/(1+np.abs(x)),axis=1)
        As,_,_=helix_matrices(sn,by)
        err=np.linalg.norm(xn[:,:3]-xn[:,7:]-np.einsum('nij,nj->ni',As,xn[:,3:6]),axis=1)
        newly=(change<1e-8)&(np.abs(ds)/(1+np.abs(s))<1e-8)&(err<1e-7)&~bad
        nit[(nit==0)&newly]=it+1
        converged|=newly
        x=xn;s=sn
        if np.all(converged|bad):break
        # Prevent runaway background candidates from contaminating other lanes.
        runaway=(np.abs(s)>1e5)|~np.isfinite(x).all(axis=1)|~np.isfinite(s)
        bad|=runaway;x[runaway]=x0[runaway];s[runaway]=0
    # Re-evaluate exactly at the final linearization point.
    A,Ap,_=helix_matrices(s,by);v=np.einsum('nij,nj->ni',Ap,x[:,3:6])
    H=np.zeros((n,3,10));H[:,:,:3]=np.eye(3);H[:,:,3:6]=-A;H[:,:,7:]=-np.eye(3)
    CH=C@H.transpose(0,2,1);W=H@CH;W[bad]=np.eye(3)
    V=np.linalg.inv(W);Vv=np.einsum('nij,nj->ni',V,v);den=np.sum(v*Vv,axis=1);den=np.maximum(den,1e-30)
    Q=V-Vv[:,:,None]*Vv[:,None,:]/den[:,None,None]
    Jx=eye-CH@Q@H; Js=np.einsum('ni,nij->nj',Vv,H)/den[:,None]
    J=np.concatenate((Jx[:,:7,:],Js[:,None,:]),axis=1)
    c8=J@C@J.transpose(0,2,1);q8=np.column_stack((x[:,:7],s))
    g=x[:,:3]-x[:,7:]-np.einsum('nij,nj->ni',A,x[:,3:6])
    b=g+np.einsum('nij,nj->ni',H,x0-x)
    chi=np.einsum('ni,nij,nj->n',b,Q,b)
    valid=converged&~bad&(np.linalg.norm(g,axis=1)<1e-7)&(x[:,6]**2>np.sum(x[:,3:6]**2,axis=1))
    return q8,c8,chi,valid,nit

def observables(q,c):
    p=np.linalg.norm(q[:,3:6],axis=1);E=q[:,6];s=q[:,7]
    m=np.sqrt(np.maximum(E*E-p*p,1e-30));L=s*p;ct=s*m
    gm=np.zeros_like(q);gm[:,3:6]=-q[:,3:6]/m[:,None];gm[:,6]=E/m
    gl=np.zeros_like(q);gl[:,3:6]=s[:,None]*q[:,3:6]/p[:,None];gl[:,7]=p
    gt=s[:,None]*gm;gt[:,7]=m
    sig=lambda g: np.sqrt(np.maximum(np.einsum('ni,nij,nj->n',g,c,g),0))
    return m,sig(gm),L,sig(gl),ct,sig(gt)

def endpoint(q,c,by,production=True):
    """Decay -> production (or inverse); S is the total signed flight parameter."""
    s=q[:,7];sgn=1 if production else -1
    A,R,Rp=helix_matrices(sgn*s,by)
    # rP=rD-A(S)pD; pP=Ap(S)pD, inverse uses -S.
    out=q.copy();out[:,:3]-=np.einsum('nij,nj->ni',A,q[:,3:6]);out[:,3:6]=np.einsum('nij,nj->ni',R,q[:,3:6])
    J=np.broadcast_to(np.eye(8),(len(q),8,8)).copy();J[:,:3,3:6]=-A;J[:,3:6,3:6]=R
    J[:,:3,7]=-sgn*np.einsum('nij,nj->ni',R,q[:,3:6]);J[:,3:6,7]=sgn*np.einsum('nij,nj->ni',Rp,q[:,3:6])
    return out,J@c@J.transpose(0,2,1)
