import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

#sapwood functions to estimate active sapwood area for sap velocity conversion

def roessler(r, tree='beech'):
    #estimate bark thickness
    #Rössler, G.: Rindenabzug richtig bemessen, Forstzeitung, 4, p. 21, 2008.
    if tree=='beech':
        db = 2.61029 + 0.28522 * 2 *r
    elif tree=='oak':
        db = 9.88855 + 0.56734 * 2 * r
    
    return db/10


def gebauer(r, tree='beech'):
    #estimate sap-wood thickness
    #Gebauer, T., Horna, V., and Leuschner, C.: Variability in radial sap flux density patterns and sapwood area among seven co-occurring temperate broad-leaved tree species, Tree Physiol., 28, 1821–1830, 2008.
    #r :: radius of tree at breast height
    r = r - roessler(r, tree='beech')/2.
    if tree=='beech':
        As = 0.778 * (2*r)**1.917
    elif tree=='oak':
        As = 0.065 * (2*r)**2.264
    return -1.*(np.sqrt((np.pi*r**2 - As)/np.pi)-r)


def gebauer_weibull(x,a,b,c,d):
    #Weibull function with 4 tree-specific parameters abcd
    return (c-1)/c + (a*((c-1)/c)**((c-1)/c))*np.exp(-1.*((x-d)/b + (c-1/c)**(1/c))**c) * ((x-d)/b + (c-1/c)**(1/c))**(c-1)

#Gebauer, T., Horna, V., and Leuschner, C.: Variability in radial sap flux density patterns and sapwood area among seven co-occurring temperate broad-leaved tree species, Tree Physiol., 28, 1821–1830, 2008.
#Parameters for 4-parameter Weibull distribution of sap velocity distribution in sapwood
gp = {
      'beech': {'name':'Fagus sylvatica', 'a': 2.69, 'b': 3.42, 'c': 1.00, 'd': 2.44},
      'hornbeam': {'name':'Carpinus betulus', 'a': 1.37, 'b': 5.88, 'c': 2.43, 'd': 2.79},
      'limeA': {'name':'Tilia sp. (A)', 'a': 1.62, 'b': 6.35, 'c': 2.71, 'd': 3.28},
      'limeB': {'name':'Tilia sp. (B)', 'a': 1.11, 'b': 4.52, 'c': 1.67, 'd': 1.88},
      'sycamoremaple': {'name':'Acer pseudoplatanus', 'a': 1.44, 'b': 8.98, 'c': 3.47, 'd': 3.42},
      'maple': {'name':'Acer campestre', 'a': 1.74, 'b': 4.86, 'c': 1.94, 'd': 2.50},
      'ash': {'name':'Fraxinus excelsior', 'a': 1.00, 'b': 1.44, 'c': 1.54, 'd': 0.42}
     }

def gebauer_rel(r, tree='beech'):
    #relative flux density as function of depth on sapwood
    #r :: radius of tree at breast height
    #returns flux distribution at 50 points along tree radius
    x=np.arange(50)/50.*gebauer(r)
    return gebauer_weibull(x,gp[tree]['a'],gp[tree]['b'],gp[tree]['c'],gp[tree]['d'])


def recko(r,hydra=False):
    #estimate sap-wood thickness
    #Račko, V., O. Mišíková, P. Hlaváč, and V. Deáková (2018), Can bark stripping cause red heartwood formation in beech stems? iForest - Biogeosciences and Forestry, 11(2), 251–258, doi:10.3832/ifor2147-011.
    #r :: radius of tree at breast height
    #returns sap-wood thickness
    if hydra:
        #only hydrated area
        return 0.34*r - 2.714
    else:
        return 0.3748*r


def galvac(r):
    #estimate sap-wood thickness
    #after Galvac et al. 1989
    #r :: radius of tree at breast height
    #returns sap-wood thickness
    dbh = r*2
    As = 0.6546*dbh**2 + 0.5736*dbh - 40.069
    return -1.*(np.sqrt((np.pi*r**2 - As)/np.pi)-r)


def gebauer_act(r,perc=0.95,tree='beech'):
    '''
    Estimate of active sapwood area based on percentile of Weibull distribution
    
    r :: radius of the tree
    perc :: percentile to localise transition to inactive sapwood 
    returns depth of active sapwood
    '''
    if type(r)==float:
        return np.where(np.cumsum(gebauer_rel(r,tree))/sum(gebauer_rel(r,tree))>perc)[0][0]/50.*gebauer(r,tree)
    else:
        dummy = r*np.nan
        for i in np.arange(len(r))[1:]:
            dummy[i] = np.where(np.cumsum(gebauer_rel(r[i],tree))/sum(gebauer_rel(r[i],tree))>perc)[0][0]/50.*gebauer(r[i],tree)
        return dummy


def sap_volume(r,s1,s2,vout=False,perc=0.95,tree='beech'):
    '''
    Fitting of Gebauer-Weibull function to measured sap velocity in order to return an estimate of the sap flux for the inner sap velocity measurement

    r :: radius of the tree
    s1 :: value or time series of sap velocity measurements at mid point of East30 sensor
    s2 :: value or time series of sap velocity measurements at inner point of East30 sensor
    tree :: trivial tree species name for Gebauer parameters
    vout :: True returns aggregated volume, False returns velocity distribution 
    '''
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
    '''
    Calculate area of circular ring
    r :: radius of tree at BH
    sens :: list of two locations marking the ring (measured from outer perimeter)
    tree :: can be beech or oak and refers to Roessler parameters
    
    returns area of circular reference sapwood
    '''
    #East30 thermocouplers location at needles 5, 18 and 30 mm
    r = r - roessler(r, tree)/2. #correct for bark
    return np.pi*((r-sens[0])**2) - np.pi*((r-sens[1])**2)


def sap_calc(SV,r,perc=0.95,tree='beech'):
    ''' 
    Wrapper for sap flow calculation
    SV :: dataframe of sap velocity (cm/h) of East30 sensors with three columns ordered 'inner', 'mid', 'outer'
    r  :: radius as breast height (cm)
    perc :: percentile for passive sapwood
    tree :: trivial tree species name for Gebauer and Roessler parameters

    returns sap flow data frame
    '''

    Sap = SV.copy()*np.nan
    colx = SV.columns[:3]
    for i in Sap.index:
        Sap.loc[i,colx[0]] = sap_volume(r,SV.loc[i,colx[1]],SV.loc[i,colx[0]],False,perc,tree)
        Sap.loc[i,colx[1]] = SV.loc[i,colx[1]]*A_circ(r,[1.1,2.4],tree)
        Sap.loc[i,colx[2]] = SV.loc[i,colx[2]]*A_circ(r,[0.,1.1],tree)

    return Sap




def stackplot(A):
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

    fill_between(A.index,A.iloc[:,2],facecolor=tableau10[0],alpha=0.7,color='b',lw=0,label=A.columns[2])
    fill_between(A.index,A.iloc[:,[1,2]].sum(axis=1),A.iloc[:,2],facecolor=tableau10[2],alpha=0.7,color='g',lw=0,label=A.columns[1])
    fill_between(A.index,A.iloc[:,:3].sum(axis=1),A.iloc[:,[1,2]].sum(axis=1),facecolor=tableau10[3],alpha=0.7,color='y',lw=0,label=A.columns[0])

