"""Independent seven-dimensional output fit.

The fit starts from ten measurements (mother q7, PV xyz). S has no measurement
or prior: eliminate its local variation by a weighted rank-one projection.
Only q7 and its MARGINAL covariance are returned. No C++ fitting routine,
8D output, or stored S covariance is used by this implementation.
"""
import numpy as np
KAPPA=0.00299792458

def geometry(s,p,by,charge):
    k=KAPPA*by*charge; angle=k*s
    a=s*np.sinc(angle/np.pi)
    b=.5*s*angle*np.sinc(angle/(2*np.pi))**2
    A=np.zeros((len(s),3,3));A[:,0,0]=A[:,2,2]=a;A[:,1,1]=s
    A[:,0,2]=b;A[:,2,0]=-b
    v=p.copy();v[:,0]=np.cos(angle)*p[:,0]+np.sin(angle)*p[:,2]
    v[:,2]=-np.sin(angle)*p[:,0]+np.cos(angle)*p[:,2]
    return A,v

def fit7(q,c,pv,pvc,by=.5,charge=1,cross=None,maxiter=30):
    n=len(q); y=np.column_stack((q,pv)); x=y.copy()
    C=np.zeros((n,10,10));C[:,:7,:7]=c;C[:,7:,7:]=pvc
    if cross is not None:C[:,:7,7:]=cross;C[:,7:,:7]=cross.transpose(0,2,1)
    s=np.einsum('ni,ni->n',q[:,:3]-pv,q[:,3:6])/np.sum(q[:,3:6]**2,axis=1)
    done=np.zeros(n,bool);bad=np.zeros(n,bool);nit=np.zeros(n,int)
    def linearize(xx,ss,cc,yy):
        A,v=geometry(ss,xx[:,3:6],by,charge)
        g=xx[:,:3]-xx[:,7:]-np.einsum('nij,nj->ni',A,xx[:,3:6])
        H=np.zeros((len(xx),3,10));H[:,:,:3]=np.eye(3);H[:,:,3:6]=-A;H[:,:,7:]=-np.eye(3)
        CH=cc@H.transpose(0,2,1); W=H@CH
        ev=np.linalg.eigvalsh(W);bad0=(ev[:,0]<=1e-14*np.maximum(ev[:,-1],1e-30))|~np.isfinite(ev).all(axis=1)
        W[bad0]=np.eye(3);V=np.linalg.inv(W)
        Vv=np.einsum('nij,nj->ni',V,v);den=np.einsum('ni,ni->n',v,Vv)
        bad0|=(den<=0)|~np.isfinite(den);den=np.where(den>0,den,1.)
        Q=V-Vv[:,:,None]*Vv[:,None,:]/den[:,None,None]
        b=g+np.einsum('nij,nj->ni',H,yy-xx)
        return H,CH,Q,b,v,Vv,den,bad0,g
    for it in range(maxiter):
        ii=np.flatnonzero(~done&~bad)
        if not len(ii):break
        H,CH,Q,b,v,Vv,den,bad0,g=linearize(x[ii],s[ii],C[ii],y[ii])
        ds=np.einsum('ni,ni->n',Vv,b)/den
        xn=y[ii]-np.einsum('nij,nj->ni',CH@Q,b);sn=s[ii]+ds
        p=np.linalg.norm(xn[:,3:6],axis=1)
        bad0|=(np.abs(sn)*p>1000)|(np.abs(KAPPA*by*charge*sn)>=np.pi)|~np.isfinite(xn).all(axis=1)
        A,_=geometry(sn,xn[:,3:6],by,charge)
        err=np.linalg.norm(xn[:,:3]-xn[:,7:]-np.einsum('nij,nj->ni',A,xn[:,3:6]),axis=1)
        done0=(np.max(abs(xn-x[ii])/(1+abs(x[ii])),axis=1)<1e-8)&(abs(ds)/(1+abs(s[ii]))<1e-8)&(err<1e-7)&~bad0
        x[ii]=xn;s[ii]=sn;bad[ii]|=bad0;done[ii]=done0;nit[ii]=it+1
    ok=done&~bad&(x[:,6]**2>np.sum(x[:,3:6]**2,axis=1))
    out=np.full((n,7),np.nan);cov=np.full((n,7,7),np.nan);chi=np.full(n,np.nan)
    ii=np.flatnonzero(ok)
    if len(ii):
        H,CH,Q,b,*_=linearize(x[ii],s[ii],C[ii],y[ii])
        J=np.eye(10)[None,:,:]-CH@Q@H
        cov[ii]=J[:,:7,:]@C[ii]@J[:,:7,:].transpose(0,2,1)
        out[ii]=x[ii,:7];chi[ii]=np.einsum('ni,nij,nj->n',b,Q,b)
    return out,cov,chi,ok,nit
