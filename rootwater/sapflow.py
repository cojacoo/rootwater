"""
The Sap flow toolbox
====================

Sap velocity and sap flow is a very interesting means to monitor xylem water 
dynamics in higher plants (i.e. trees in our case). Since a mere flow conversion
is likely to result in erroneous assumptions (Čermák et al., 2004), we have 
compiled a couple of sapwood-related functions to estimate active sapwood area 
for sap velocity conversion to sap flow.

.. note::
    This toolbox is by no means complete nor exhaustive. Please regard it as 
    helper functions which require throughout testing and deserve substantial
    extension to further application cases.

.. note::
    To get started and for direct application (tested for beech trees) use rootwater.sapflow.sap_calc
    and provide a pandas.DataFrame with measured sap velocity from East30 sensors 
    (in cm/h).

References
----------
Čermák, J., J. Kučera, and N. Nadezhdina (2004), Sap flow measurements with some 
thermodynamic methods, flow integration within trees and scaling up from sample 
trees to entire forest stands, Trees, 18(5), 529–546, doi:10.1007/s00468-004-0339-6.

"""
import json
import sys
import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

sys.path.append(os.path.dirname(__file__))
from gebauer_params import gp
#gp : dictionary for all valid tree names. Each name 
#    has to be key to a nested dict that defines the 
#    four Weibull parameters a,b,c,d


def roessler(r, tree='beech'):
    r"""Estimate bark thickness

    Returns the estimated bark thickness as published by Rössler (2008)

    Parameters
    ----------
    r : float
        tree radius at breast height (in cm)
    tree : str
        Tree name, for which to calculate bark thickness.
        Can be one of ['beech', 'oak']

    Returns
    -------
    db : float
        bark thickness (in mm)
    
    Raises
    ------
    NotImplementedError : if tree is not in ('oak', 'beech')

    References
    ----------
    Rössler, G.: Rindenabzug richtig bemessen, Forstzeitung, 4, p. 21, 2008.
    """
    if tree=='beech':
        db = 2.61029 + 0.28522 * 2 * r
    elif tree=='oak':
        db = 9.88855 + 0.56734 * 2 * r
    else:
        raise NotImplementedError("Only 'beech' and 'oak' supported.")
    
    return db / 10


def gebauer(r, tree='beech'):
    r"""Sap-wood thickness

    Calculates sap-wood thickness as published by Gebauer et al. (2008)

    Parameters
    ----------
    r : float
        tree radius at breast height (in cm)
    tree : str
        Tree name, for which to calculate bark and sapwood thickness.
        Can be one of ['beech', 'oak']

    Returns
    -------
    th : float
        sap-wood thickness (in mm)

    Raises
    ------
    NotImplementedError : if tree is not in ('oak', 'beech')

    References
    ----------
    Gebauer, T., Horna, V., and Leuschner, C.: Variability in radial sap flux
    density patterns and sapwood area among seven co-occurring temperate 
    broad-leaved tree species, Tree Physiol., 28, 1821–1830, 2008.

    """
    r = r - roessler(r, tree=tree) / 2.
    
    if tree=='beech':
        As = 0.778 * (2*r)**1.917
    elif tree=='oak':
        As = 0.065 * (2*r)**2.264
    else:
        raise NotImplementedError("Only 'beech' and 'oak' supported.")
    
    return -1.*(np.sqrt((np.pi*r**2 - As)/np.pi)-r)


def gebauer_weibull(x,a,b,c,d):
    r"""4-parameter Weibull function after Gebauer

    Calls Weibull distribution function as published by 
    Gebauer et al. (2008).

    Parameters
    ----------
    x : float or numpy.ndarray
        realtive sampling points for distribution function
    a, b, c, d : float
        Weibull function parameters, which are tree-specific in this case

    References
    ----------
    Gebauer, T., Horna, V., and Leuschner, C.: Variability in radial sap flux
    density patterns and sapwood area among seven co-occurring temperate 
    broad-leaved tree species, Tree Physiol., 28, 1821–1830, 2008.

    """
    #Weibull function with 4 tree-specific parameters abcd
    return (c-1)/c + (a*((c-1)/c)**((c-1)/c))*np.exp(-1.*((x-d)/b + (c-1/c)**(1/c))**c) * ((x-d)/b + (c-1/c)**(1/c))**(c-1)


def get_default_gp():
    r"""read default gp

    Loads default Weibull distribution parameters as published by 
    Gebauer et al. (2008). Can be used by 
    :func:`gebauer_weibull <rootwater.sf.gebauer_weibull>` to 
    calculate sap velocity distribution in sapwood

    References
    ----------
    Gebauer, T., Horna, V., and Leuschner, C.: Variability in radial sap flux
    density patterns and sapwood area among seven co-occurring temperate 
    broad-leaved tree species, Tree Physiol., 28, 1821–1830, 2008.

    """
    with open(GP_PATH, 'r') as fs:
        params = json.load(fs)
    
    # if needed, you could validate the params here
    return params


def gebauer_rel(r, tree='beech', n_points=50):
    r"""relative flux density

    Calculates relative flux density as a function 
    of depth on sapwood for n_points. 

    Parameters
    ----------
    r : float
        tree radius at breast height (in cm)
    tree : str
        Tree name, for which to calculate Weibull function.
        Tree name has to be in gp.keys()
    n_points : int
        Number of points for solving Weibull. 
        This is the resolution over depth.

    Returns
    -------
    sv : numpy.ndarray
        relative flux density at n_points

    References
    ----------
    Gebauer, T., Horna, V., and Leuschner, C.: Variability in radial sap flux
    density patterns and sapwood area among seven co-occurring temperate 
    broad-leaved tree species, Tree Physiol., 28, 1821–1830, 2008.

    """
    
    x=np.arange(n_points)/ n_points *gebauer(r)
    
    p = gp.get(tree)

    if p is None:
        raise ValueError('Tree %s is unknown' % tree)
    return gebauer_weibull(x, p['a'], p['b'], p['c'], p['d'])


def recko(r,hydra=False):
    r"""sap-wood thickness after Račko et al. (2018)

    Calculates sap-wood thickness as a function of the tree radius 
    at breast height.

    Parameters
    ----------
    r : float
        tree radius at breast height (in cm)
    hydra : bool
        selects if only the hydrated area is returned (when True)

    Returns
    -------
    sapwood : float
        sapwood thickness (in cm)
    
    References
    ----------
    Račko, V., O. Mišíková, P. Hlaváč, and V. Deáková (2018), 
    Can bark stripping cause red heartwood formation in beech stems? 
    iForest - Biogeosciences and Forestry, 11(2), 251–258, doi:10.3832/ifor2147-011.
    """
    
    if hydra:
        #only hydrated area
        return 0.34*r - 2.714
    else:
        return 0.3748*r


def galvac(r):
    r"""sap-wood thickness after Galvac et al. (1990)

    Calculates sap-wood thickness as a function of the tree radius 
    at breast height.

    Parameters
    ----------
    r : float
        tree radius at breast height (in cm)
    
    Returns
    -------
    sapwood : float
        sapwood thickness (in cm)
    
    References
    ----------
    Glavac, V., Koenies, H. & Ebben, U. Holz als Roh- und Werkstoff (1990) 48: 437. https://doi.org/10.1007/BF02627628
    """
    
    #returns sap-wood thickness
    dbh = r*2
    As = 0.6546*dbh**2 + 0.5736*dbh - 40.069
    return -1.*(np.sqrt((np.pi*r**2 - As)/np.pi)-r)


def gebauer_act(r,perc=0.95,tree='beech'):
    r"""Active sapwood area based on percentile of Weibull distribution

    Calculates the "zero" sap velocity limit as given percentile of relative 
    flux velocity distribution as a Weibull function after Gebauer et al. (2008)

    Parameters
    ----------
    r : float or numpy.ndarray
        tree radius at breast height (in cm)
    perc : float
        percentile to define the "zero" sap velocity limit
    tree : str
        Tree name, for which to calculate Weibull function.
        Tree name has to be in gp.keys()
    
    Returns
    -------
    act_sap : float or numpy.ndarray
        depth of "zero" sap velocity limit as active sapwood in tree (in cm)

    References
    ----------
    Gebauer, T., Horna, V., and Leuschner, C.: Variability in radial sap flux
    density patterns and sapwood area among seven co-occurring temperate 
    broad-leaved tree species, Tree Physiol., 28, 1821–1830, 2008.

    """
    if type(r)==float:
        return np.where(np.cumsum(gebauer_rel(r,tree))/sum(gebauer_rel(r,tree))>perc)[0][0]/50.*gebauer(r,tree)
    else:
        act_sap = r*np.nan
        for i in np.arange(len(r))[1:]:
            act_sap[i] = np.where(np.cumsum(gebauer_rel(r[i],tree))/sum(gebauer_rel(r[i],tree))>perc)[0][0]/50.*gebauer(r[i],tree)
        return act_sap


def sap_volume(r,s1,s2,vout=False,perc=0.95,tree='beech'):
    r"""Estimate sap flow from sap velocity in inner sapwood measured with East30 sensors

    Calculates the sap flow after Gebauer et al. (2008) based on sap velocity measurements 
    by fitting of Gebauer-Weibull function to measured sap velocity at mid and inner point
    through a scaling factor (but not changing the empirical, tree-specific parameters).

    Parameters
    ----------
    r : float
        tree radius at breast height (in cm)
    s1 : float or pandas.Series (datetime index is preferable)
        sap velocity measurement at mid point of East30 sensor (cm/time)
    s2 : float or pandas.Series (datetime index is preferable)
        sap velocity measurement at inner point of East30 sensor (cm/time)
    vout : bool
        True returns aggregated volume flux (cm3/time), False returns velocity distribution (cm/time)
    perc : float
        percentile to define the "zero" sap velocity limit
    tree : str
        Tree name, for which to calculate bark thickness and Weibull function.
        Tree name has to be in gp.keys()

    Returns
    -------
    return : float or pandas.Series
        aggregated volume flux (cm3/time) (if vout is True), or
        returns velocity distribution (cm/time) (if vout is False)

    References
    ----------
    Gebauer, T., Horna, V., and Leuschner, C.: Variability in radial sap flux
    density patterns and sapwood area among seven co-occurring temperate 
    broad-leaved tree species, Tree Physiol., 28, 1821–1830, 2008.

    """

    xi = np.arange(50)/50.*gebauer(r,tree)
    
    def aply_geb(s1x):
        dummy = s1x*gebauer_rel(r,tree)
        er1 = (dummy[np.arange(50)/50.*gebauer(r,tree) >= 1.8][0]-s1)**2
        er2 = (dummy[np.arange(50)/50.*gebauer(r,tree) >= 3.][0]-s2)**2
        return np.sqrt(0.2*er1+er2)
    
    res = minimize_scalar(aply_geb)
    dummy = res.x*gebauer_rel(r,tree)
    v3 = dummy[(np.arange(50)/50.*gebauer(r,tree) > 2.4) & (np.arange(50)/50.*gebauer(r,tree) <= gebauer_act(r,perc,tree))]
    
    rx = (np.arange(50)/50.*gebauer(r,tree))[(np.arange(50)/50.*gebauer(r,tree) > 2.4) & (np.arange(50)/50.*gebauer(r,tree) <= gebauer_act(r,perc,tree))]
    Ax = rx*np.nan
    for i in np.arange(len(Ax)):
        Ax = A_circ(r,[rx[i]-0.01*gebauer(r,tree),rx[i]+0.01*gebauer(r,tree)])
    
    if vout:
        return dummy
    else:
        return np.sum(Ax * v3)


def A_circ(r,sens=[0.,1.1],tree='beech'):
    r"""Calculate area of circular ring

    Simple geometrical calculation of a circular ring area as reference cross-section
    for sap flow calculation.

    Parameters
    ----------
    r : float
        tree radius at breast height (in cm)
    sens : list of floats
        outer and inner point of ring (in cm)
    tree : str
        Tree name, for which to calculate and subtract bark thickness.
        Can be one of ['beech', 'oak'].

    Returns
    -------
    return : float
        area of ring as 

    """
    
    #East30 thermocouplers location at needles 5, 18 and 30 mm
    r = r - roessler(r, tree)/2. #correct for bark
    return np.pi*((r-sens[0])**2) - np.pi*((r-sens[1])**2)


def sap_calc(SV,r,perc=0.95,tree='beech'):
    r"""Wrapper for sap flow calculation with rootwater.sapflow.sap_volume

    Calculates the sap flow after Gebauer et al. (2008) based on measured sap velocity 
    with East30 sensors.
    
    Parameters
    ----------
    SV : pandas.DataFrame
        sap velocity (in cm/h) in three columns ordered inner, mid, outer point
    r : float
        tree radius at breast height (in cm)
    perc : float
        percentile to define the "zero" sap velocity limit
    tree : str
        Tree name, for which to calculate bark thickness and Weibull function.
        Tree name has to be in gp.keys()

    Returns
    -------
    return : pandas.DataFrame
        sap volume flux (cm3/h)

    """
    
    Sap = SV.copy()*np.nan
    colx = SV.columns[:3]
    for i in Sap.index:
        Sap.loc[i,colx[0]] = sap_volume(r,SV.loc[i,colx[1]],SV.loc[i,colx[0]],False,perc,tree)
        Sap.loc[i,colx[1]] = SV.loc[i,colx[1]]*A_circ(r,[1.1,2.4],tree)
        Sap.loc[i,colx[2]] = SV.loc[i,colx[2]]*A_circ(r,[0.,1.1],tree)

    return Sap




def stackplot(A):
    import matplotlib.pyplot as plt
    r"""plot stacked time series (of first three columns of the provided dataframe)

    
    Parameters
    ----------
    A : pandas.DataFrame (preferrably with datetime index)
        DataFrame with three columns to be stacked
    
    Returns
    -------
    plot

    """
    
    #plot stacked time series of first three columns of the dataframe A
    # These are the "Tableau 20" colors as RGB.  
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.  
    for i in range(len(tableau20)):  
        r, g, b = tableau20[i]  
        tableau20[i] = (r / 255., g / 255., b / 255.)  
    
    tableau10=tableau20[0::2]

    plt.fill_between(A.index,A.iloc[:,2],facecolor=tableau10[0],alpha=0.7,color='b',lw=0,label=A.columns[2])
    plt.fill_between(A.index,A.iloc[:,[1,2]].sum(axis=1),A.iloc[:,2],facecolor=tableau10[2],alpha=0.7,color='g',lw=0,label=A.columns[1])
    plt.fill_between(A.index,A.iloc[:,:3].sum(axis=1),A.iloc[:,[1,2]].sum(axis=1),facecolor=tableau10[3],alpha=0.7,color='y',lw=0,label=A.columns[0])

