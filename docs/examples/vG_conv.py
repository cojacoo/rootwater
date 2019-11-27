#Van Genuchten Conversions
#(cc) jackisch@kit.edu

import numpy as np
import pandas as pd

# standard parameters after Carsel & Parrish 1988
carsel=pd.DataFrame(
[[  'C', 30.,  15.,  55.,   0.068,   0.38,   0.008*100.,   1.09,    0.200/360000.],
[  'CL', 37.,  30.,  33.,   0.095,   0.41,   0.019*100.,   1.31,    0.258/360000.],
[   'L', 40.,  40.,  20.,   0.078,   0.43,   0.036*100.,   1.56,    1.042/360000.],
[  'LS', 13.,  81.,   6.,   0.057,   0.43,   0.124*100.,   2.28,   14.592/360000.],
[   'S',  4.,  93.,   3.,   0.045,   0.43,   0.145*100.,   2.68,   29.700/360000.],
[  'SC', 11.,  48.,  41.,   0.100,   0.38,   0.027*100.,   1.23,    0.121/360000.],
[ 'SCL', 19.,  54.,  27.,   0.100,   0.39,   0.059*100.,   1.48,    1.308/360000.],
[  'SI', 85.,   6.,   9.,   0.034,   0.46,   0.016*100.,   1.37,    0.250/360000.],
[ 'SIC', 48.,   6.,  46.,   0.070,   0.36,   0.005*100.,   1.09,    0.021/360000.],
['SICL', 59.,   8.,  33.,   0.089,   0.43,   0.010*100.,   1.23,    0.071/360000.],
[ 'SIL', 65.,  17.,  18.,   0.067,   0.45,   0.020*100.,   1.41,    0.450/360000.],
[  'SL', 26.,  63.,  11.,   0.065,   0.41,   0.075*100.,   1.89,    4.421/360000.]],
columns=['Typ','Silt','Sand','Clay','thr','ths','alpha','n','ks'],index=np.arange(12).astype(int)+1)

#DEBUG: in the calculations type errors occur.
#This should be fixed some time. However, I override warnings for now:
np.seterr(all='ignore')

# conversions
def ku_psi(psi, ks, alpha, n, m=None, l=0.5):
    #Calculate unsaturated hydraulic conductivity (ku) from matrix head (psi)
    if m is None:
        m=1.-1./n
    v = 1. + (alpha*np.abs(psi))**n
    ku = ks* v**(-1.*m*l) * (1. - (1. - 1/v)**m)**2
    return ku

def ku_thst(thst, ks, alpha, n, m=None, l=0.5):
    #Calculate unsaturated hydraulic conductivity (ku) relative saturation (theta*)
    if m is None:
        m=1.-1./n
    ku = ks*thst**l * (1 - (1-thst**(1/m))**m)**2#
    return ku

def ku_theta(theta, ths, thr, ks, alpha, n, m=None):
    #Calculate unsaturated hydraulic conductivity (ku) from matrix head (psi)
    if m is None:
        m=1.-1./n
    th_star=thst_theta(theta,ths,thr)
    ku = ku_thst(th_star,ks, alpha, n, m)
    return ku

def thst_theta(theta,ths,thr):
    #Calculate relative saturation (theta*) from soil moisture (theta)
    th_star=(theta-thr)/(ths-thr) #
    return th_star

def theta_thst(th_star,ths,thr):
    #Calculate soil moisture (theta) from relative saturation (theta*)
    theta=th_star*(ths-thr)+thr
    return theta

def theta_psi(psi,ths,thr,alpha,n,m=None):
    #Calculate soil moisture (theta) from matrix head (psi)
    if m is None:
        m=1.-1./n
    theta=theta_thst(thst_psi(psi,alpha,n,m),ths,thr)
    return theta

def psi_thst(th_star,alpha,n,m=None):
    #Calculate matrix head (psi) from relative saturation (theta*)
    if m is None:
        m=1.-1./n
    psi = -1./alpha * ( (1-th_star**(1/m))/(th_star**(1/m)) )**(1./n)
    if (np.iterable(psi)):
        if any(np.isinf(psi)):
            if type(alpha)==float:
                psi[np.isinf(psi)]=(-1./alpha * ( (1-0.98**(1/m))/(0.98**(1/m)) )**(1./n))
            else:
                psi[np.isinf(psi)]=(-1./alpha * ( (1-0.98**(1/m))/(0.98**(1/m)) )**(1./n))[np.isinf(psi)]
    return psi

def psi_theta(theta,ths,thr,alpha,n,m=None):
    #Calculate matrix head (psi) from soil moisture (theta)
    if m is None:
        m=1.-1./n
    th_star=thst_theta(theta,ths,thr)
    psi= -1. * ( (1 - th_star**(1./m)) / (th_star**(1./m)) )**(1./n) / alpha
    if (np.iterable(psi)):
        if any(np.isinf(psi)):
            if type(alpha)==float:
                psi[np.isinf(psi)]=(-1./alpha * ( (1-0.98**(1/m))/(0.98**(1/m)) )**(1./n))
            else:
                psi[np.isinf(psi)]=(-1./alpha * ( (1-0.98**(1/m))/(0.98**(1/m)) )**(1./n))[np.isinf(psi)]
    return psi

def thst_psi(psi,alpha,n,m=None):
    #Calculate relative saturation (theta*) from matrix head (psi)
    if m is None:
        m=1.-1./n
    th_star = (1./(1.+(np.abs(psi)*alpha)**n))**m#
    return th_star

def c_psi(psi,ths,thr,alpha,n,m=None):
    #Calculate water capacity (c) from matrix head (psi)
    if m is None:
        m=1.-1./n
    c=-1.*(ths-thr)*n*m* alpha**n * np.abs(psi)**(n-1.) * (1+(alpha*np.abs(psi))**n)**(-1.*m - 1.)
    #y=m*(1./(1+np.abs(psi*alpha)**n))**(m+1.) * n *(np.abs(psi)*alpha)**(n-1.) *alpha
    #c=(ths-thr)*y
    return c

def dpsidtheta_thst(th_star,ths,thr,alpha,n,m=None):
    #Calculate matrix head (psi) from relative saturation (theta*)
    if m is None:
        m=1.-1./n
    if (type(th_star)==float) | (type(th_star)==np.float64):
        if th_star>0.9899:
            th_star=0.9899
        if th_star<0.01:
            th_star=0.01
    else:
        th_star[th_star>0.9899]=0.9899
        th_star[th_star<0.01]=0.01
    
    th_star1=th_star-0.01
    th_star+=0.01
    psi = -1./alpha * ( (1-th_star**(1/m))/(th_star**(1/m)) )**(1./n)
    psist = -1./alpha * ( (1-th_star1**(1/m))/(th_star1**(1/m)) )**(1./n)

    if np.iterable(psi):
        psi[np.isinf(psi)]=-1./alpha * ( (1-0.98**(1/m))/(0.98**(1/m)) )**(1./n)
    if np.iterable(psist):
        psist[np.isinf(psist)]=-1./alpha * ( (1-0.98**(1/m))/(0.98**(1/m)) )**(1./n)
    
    theta=th_star*(ths-thr)+thr
    thetast=th_star1*(ths-thr)+thr

    dpsidtheta=(psist-psi)/(thetast-theta)
    return dpsidtheta

def D_psi(psi,ks,ths,thr,alpha,n,m=None):
    #Calculate diffusivity (D) from matrix head (psi)
    if m is None:
        m=1.-1./n
    psix=np.array([psi-0.05,psi,psi+0.05])
    if isinstance(ks,np.float):
        kus=ku_psi(psix, ks, alpha, n, m)
        dth=np.diff(theta_thst(thst_psi(psix,alpha,n,m),ths,thr),axis=0)[0]
    else:
        kus=ku_psi(psix, ks.repeat(3).reshape(np.shape(psix)), alpha.repeat(3).reshape(np.shape(psix)), n.repeat(3).reshape(np.shape(psix)), m.repeat(3).reshape(np.shape(psix)))
        dth=np.diff(theta_thst(thst_psi(psix,alpha.repeat(3).reshape(np.shape(psix)),n.repeat(3).reshape(np.shape(psix)),m.repeat(3).reshape(np.shape(psix))),ths.repeat(3).reshape(np.shape(psix)),thr.repeat(3).reshape(np.shape(psix))),axis=0)[0]
    
    if len(np.shape(kus))==1:
        D=kus[1]*0.1/dth
    else:
        D=kus[1,:]*0.1/dth
    return D

def dDdtheta_thst(th_star,ths,thr,ks,alpha,n,m=None):
    #Calculate matrix head (psi) from relative saturation (theta*)
    if m is None:
        m=1.-1./n
    if (type(th_star)==float) | (type(th_star)==np.float64):
        if th_star>0.9899:
            th_star=0.9899
        if th_star<0.01:
            th_star=0.01
    else:
        th_star[th_star>0.9899]=0.9899
        th_star[th_star<0.01]=0.01
    
    th_star1=th_star-0.01
    th_star+=0.01
    D=D_thst(th_star,ths,thr,ks,alpha,n,m)
    Dst=D_thst(th_star1,ths,thr,ks,alpha,n,m)
        
    theta=th_star*(ths-thr)+thr
    thetast=th_star1*(ths-thr)+thr

    dDdtheta=(Dst-D)/(thetast-theta)
    return dDdtheta



def D_theta(theta,ths,thr,ks,alpha,n,m=None):
    #Calculate diffusivity (D) from soil moisture (theta)
    if m is None:
        m=1.-1./n

    the=(theta-thr)/(ths-thr)
    Dd=(ks*(1-m)*(the**(0.5-(1/m)))) / (alpha*m*(ths-thr)) *( (1-the**(1/m))**(-1*m) + (1-the**(1/m))**m -2 )
    return Dd

def D_thst(thst,ths,thr,ks,alpha,n,m=None):
    #Calculate diffusivity (D) from soil moisture (theta)
    if m is None:
        m=1.-1./n

    Dd=(ks*(1.-m)*(thst**(0.5-(1./m)))) / (alpha*m*(ths-thr)) *( (1.-thst**(1./m))**(-1.*m) + (1.-thst**(1./m))**m -2. )
    return Dd

def dcst_thst(thst,ths,thr,ks,alpha,n,m=None):
    #Calculate diffusivity (D as ku*dpsi/dthst) from soil moisture (theta)
    if m is None:
        m=1.-1./n
    c=c_psi(psi_thst(thst,alpha,n,m),ths,thr,alpha,n,m)
    ku=ku_thst(thst, ks, alpha,n,m)
    D=-ku/(c*theta_thst(thst,ths,thr))
    return D


# wrapper
def th_psi_f(psi,sample,mc):
    if (isinstance(psi,pd.DataFrame)) or (isinstance(psi,pd.Series)):
        theta=theta_psi(psi.values, mc.soilmatrix.ts[sample].values, mc.soilmatrix.tr[sample].values, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    else:
        theta=theta_psi(psi, mc.soilmatrix.ts[sample].values, mc.soilmatrix.tr[sample].values, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    return theta

def psi_th_f(theta,sample,mc):
    if (isinstance(theta,pd.DataFrame)) or (isinstance(theta,pd.Series)):
        psi=psi_theta(theta.values, mc.soilmatrix.ts[sample].values, mc.soilmatrix.tr[sample].values, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    else:
        psi=psi_theta(theta, mc.soilmatrix.ts[sample].values, mc.soilmatrix.tr[sample].values, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    return psi

def psi_ths_f(thst,sample,mc):
    if (isinstance(thst,pd.DataFrame)) or (isinstance(thst,pd.Series)):
        psi=psi_thst(thst.values, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    else:
        psi=psi_thst(thst, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    return psi

def D_psi_f(psi,sample,mc):
    if (isinstance(psi,pd.DataFrame)) or (isinstance(psi,pd.Series)):
        D=D_psi(psi.values, mc.soilmatrix.ks[sample].values, mc.soilmatrix.ts[sample].values, mc.soilmatrix.tr[sample].values, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    else:
        D=D_psi(psi, mc.soilmatrix.ks[sample].values, mc.soilmatrix.ts[sample].values, mc.soilmatrix.tr[sample].values, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    return D

def D_thst_f(thst,sample,mc):
    if (isinstance(thst,pd.DataFrame)) or (isinstance(thst,pd.Series)):
        D=D_thst(thst.values, mc.soilmatrix.ts[sample].values, mc.soilmatrix.tr[sample].values, mc.soilmatrix.ks[sample].values, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    else:
        D=D_thst(thst, mc.soilmatrix.ts[sample].values, mc.soilmatrix.tr[sample].values, mc.soilmatrix.ks[sample].values, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    return D

def ku_psi_f(psi,sample,mc):
    if (isinstance(psi,pd.DataFrame)) or (isinstance(psi,pd.Series)):
        ku=ku_psi(psi.values, mc.soilmatrix.ks[sample].values, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    else:
        ku=ku_psi(psi, mc.soilmatrix.ks[sample].values, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    return ku

def Dku_thst_f(thst,sample,mc):
    if (isinstance(thst,pd.DataFrame)) or (isinstance(thst,pd.Series)):
        psi=psi_thst(thst.values,mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    else:
        psi=psi_thst(thst,mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)

    ku=ku_psi(psi, mc.soilmatrix.ks[sample].values, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    D=D_psi(psi, mc.soilmatrix.ks[sample].values, mc.soilmatrix.ts[sample].values, mc.soilmatrix.tr[sample].values, mc.soilmatrix.alpha[sample].values, mc.soilmatrix.n[sample].values)
    theta=theta_thst(thst,mc.soilmatrix.ts[sample].values, mc.soilmatrix.tr[sample].values)
    return (D, ku, theta)


# create look-up tables
def create_lookup(ths,thr,ks,alpha,n,m=None):
    n_soils=len(ths)
    if m is None:
        m=1.-1./n

    psi=np.empty((100,n_soils))
    theta=np.empty((100,n_soils))
    ku=np.empty((100,n_soils))
    D=np.empty((100,n_soils))
    for i in np.arange(n_soils):
        psi[:,i]=psi_thst(np.arange(0.01,1.01,0.01),alpha[i],n[i],m[i])
        theta[:,i]=theta_thst(np.arange(0.01,1.01,0.01),ths[i],thr[i])
        ku[:,i]=ku_thst(np.arange(0.01,1.01,0.01),ks[i],alpha[i],n[i],m[i])
        D[:,i]=D_thst(np.arange(0.01,1.01,0.01),ths[i],thr[i],ks[i],alpha[i],n[i],m[i])

    D[-1,:]=D[-2,:]
    return [psi,theta,ku,D]




