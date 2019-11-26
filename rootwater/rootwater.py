# coding=utf-8
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf 
import scipy.ndimage.filters as spf
import datetime
from astral import Location
import hydroeval as he

# function to calculate change in soil moisture as root water uptake

def fRWU(ts,lat=49.70764, lon=5.897638, diffx=3, slope_diff=3, maxdiffs=0.25, mintime=3.5):
    '''
    Calulate a daily root water uptake estimate from a soil moisture time series
    (cc) c.jackisch@tu-braunschweig.de
    
    ts :: time series of one soil moisture sensor (assumes vol.%) as a pandas series (requires high temporal resolution and a time zone aware index)
    lat/lon :: latitude/longitude of location (degree)
    
    further parameters:
    diffx :: number of time steps to evaluate change in moisture to
    slope_diff :: min difference factor of slope between night and day lin. regession (in case of night decrease)
    maxdiffs :: max of soil moistue difference to assume no significant other water transport
    mintime :: minmimal time of a day or night period (h)
    
    returns a data frame with time series of daily RWU estimates and daily references:
    rwu :: root water uptake with extrapolated night changes
    rwu_nonight :: neglecting nocturnal changes
    lm_night :: slope of linear model during night
    lm_day :: slope of linear model during day
    step_control :: control values (1111 means all criteria met)
    evalx :: control values for time references
    eval_nse :: control values for diurnal step shape as nash-sutcliffe efficiency 
    tin :: start of previous night
    tout :: start of day 
    tix :: start of next night
    '''
    
    # use astral to get sunrise/sunset time references as a function of the date    
    l = Location()
    l.latitude = lat
    l.longitude = lon
    l.timezone = str(ts.index.tz)
    
    #sunrise sunset
    def sunr(dd):
        # give date and return time of sunrise
        sunrise = pd.to_datetime(l.sun(date=dd)['sunrise'])
        return sunrise
        
    def suns(dd):
        # give date and return time of sunset
        sunset = pd.to_datetime(l.sun(date=dd)['sunset'])
        return sunset

    # get unique days in time series
    ddx = ts.resample('1d').mean().index.date
    
    # get frequencies of ts
    freqx = (pd.Series(ts.index[1:]) - pd.Series(ts.index[:-1])).value_counts()
        
    # get change in soil moisture as smoothed diff
    dif_ts = pd.Series(spf.gaussian_filter1d(ts.diff(diffx),1))
    dif_ts.index = ts.index
    
    # create empty dataframe for RWU calculation and evaluation
    RWU = pd.DataFrame(np.zeros((len(ddx),10))*np.nan)
    RWU.index = pd.to_datetime(ddx)
    RWU.columns = ['rwu','rwu_nonight','lm_night','lm_day','step_control','evalx','eval_nse','tin','tout','tix']
    
    def startstopRWU(dd):
        # give soilmoisture ts and date, return time of end of RWU
        try:
            # find first steps of not declining soil moisture after sunset
            tsx = dif_ts.loc[suns(dd-datetime.timedelta(hours=24))-datetime.timedelta(hours=5):suns(dd)+datetime.timedelta(hours=2)]
            stopRWU = tsx.index[np.where(tsx.values >=0)[0][0]]
            startRWU = tsx.index[np.where(tsx.values <=-0.02)[0][np.where(tsx.values <=-0.02)[0]>np.where(tsx.values >=0)[0][0]][0]-1]
            stop2RWU = tsx.index[np.where(tsx.values > 0)[0][np.where(tsx.values > 0)[0] > np.where(tsx.values <=-0.02)[0][np.where(tsx.values <=-0.02)[0]>np.where(tsx.values >=0)[0][0]][0]][0]]
            return [stopRWU,startRWU,stop2RWU, 1]
        except:
            # in case soil moisture keeps falling without stepping assume 1 hour after sunset/sunrise but return warning flag
            stopRWU = ts.index[ts.index.get_loc(suns(dd-datetime.timedelta(hours=24))+datetime.timedelta(hours=1), method='nearest')]
            startRWU = ts.index[ts.index.get_loc(sunr(dd)+datetime.timedelta(hours=1), method='nearest')]
            stop2RWU = ts.index[ts.index.get_loc(suns(dd)+datetime.timedelta(hours=1), method='nearest')]
            return [stopRWU,startRWU,stop2RWU, 0]
    
    def idstep_startstop(dd):
        # return astro reference times for idealised step
        tin = ts.index[ts.index.get_loc(suns(dd-datetime.timedelta(hours=24))-datetime.timedelta(hours=1.5), method='nearest')]
        tout = ts.index[ts.index.get_loc(sunr(dd)+datetime.timedelta(hours=2), method='nearest')]
        tix = ts.index[ts.index.get_loc(suns(dd)-datetime.timedelta(hours=0.5), method='nearest')]
        return [tin,tout,tix]
    
    def dayRWU(dd):
        # get reference times
        [tin,tout,tix,evalx] = startstopRWU(dd)
        
        # check for soil moisture differences and min time spans
        if ((tout-tin).seconds<mintime*3600.) | ((tix-tout).seconds<mintime*3600.):
            return [np.nan, np.nan, np.nan, np.nan, 2, evalx, tin, tout, tix]
        if any(dif_ts.loc[tin:tix]>maxdiffs):
            return [np.nan, np.nan, np.nan, np.nan, 3, evalx, tin, tout, tix]
        
        # build linear extrapolation model of night time change
        df=pd.DataFrame([np.arange(len(ts.loc[tin:tout-datetime.timedelta(hours=1)])),ts.loc[tin:tout--datetime.timedelta(hours=1)].values]).T 
        df.columns=['x','y']
        try:
            mod = smf.ols(formula='y ~ x', data=df) 
            res = mod.fit()
        except:
            return [np.nan, np.nan, np.nan, np.nan, 0, evalx, tin, tout, tix]
                   
        
        # build linear model of day time change
        df2=pd.DataFrame([np.arange(len(ts.loc[tout:tix])),ts.loc[tout:tix].values]).T 
        df2.columns=['x','y'] 
        try:
            mod2 = smf.ols(formula='y ~ x', data=df2) 
            res2 = mod2.fit()
        except:
            return [np.nan, np.nan, res.params.x, np.nan, 0, evalx, tin, tout, tix]
        
        # create dummy time series for night time extrapolation
        dummy = pd.date_range(tin,tix, freq=freqx.index[0])
        fuse = pd.Series(data = res.params.Intercept+res.params.x*np.arange(len(dummy)), index = dummy)
        
        # control of assumptions of a step
        step_control = 0
        if res.params.x / ((6.*3600.)/freqx.index[0].seconds) > -0.5/6.: #night slope shall be more than minus 0.5 vol.% per 6h
            step_control += 10
        if res.params.x / ((6.*3600.)/freqx.index[0].seconds) < 1/6.: #night slope shall be less than plus 1 vol.% per 6h
            step_control += 100
        if res2.params.x < 0: #day slope must be negative 
            step_control += 1000
        if res2.params.x < slope_diff*res.params.x: #day slope must be at least 3 times more steep than night (if night was negative)
            step_control += 1
        
        # only if step_control == 1111 fully valid results:
        rwu = fuse.loc[tix]-ts.loc[tix]
        rwu_nonight = ts.loc[tout]-ts.loc[tix]
        return [rwu, rwu_nonight, res.params.x, res2.params.x, step_control, evalx, tin,tout,tix]
        
    def dayRWU2(dd,crit_nse=0.5):
        # perform comparison to idealised step before evaluation
        # get reference times
        [tin,tout,tix,evalx] = startstopRWU(dd)
        [dtin,dtout,dtix] = idstep_startstop(dd)

        # construct idealised step reference
        idx = pd.date_range(dtin, dtix, freq='30min')
        dummy = pd.Series(np.zeros(len(idx))*np.nan,index = idx)
        
        dummy[dtin] = ts.loc[dtin]
        dummy[dtout+datetime.timedelta(hours=1)] = ts.loc[dtin]+0.01
        dummy[dtix-datetime.timedelta(hours=2.5)] = ts.loc[dtix]
        
        dummy = dummy.interpolate()
        dummyx = pd.concat([dummy,ts.loc[tin-datetime.timedelta(hours=0.5):tix+datetime.timedelta(hours=0.5)]],axis=1)
        dummyx.columns = ['ideal','obs']
        
        dummyx = dummyx.dropna()
        
        # compare observed soil moisture dynamics with idealised step
        evaly = he.nse_c2m(dummyx.obs.values,dummyx.ideal.values)

        #if evaly >= crit_nse:
        [rwu, rwu_nonight, resparamsx, res2paramsx, step_control, evalx, tin2,tout2,tix2] = dayRWU(dd)

        return [rwu, rwu_nonight, resparamsx, res2paramsx, step_control, evalx, evaly, tin,tout,tix]
        
    for dd in ddx[:-1]:
        RWU.loc[dd] = dayRWU2(dd)
    
    return RWU

def dfRWUc(dummyd,tz='Etc/GMT-1',safeRWU=True):
    '''
    Wrapper to quickly apply fRWU to a dataframe with soil moisture values.
    dummyd :: input data frame with time stamp index and columns of soil moisture
    tz :: tiime zone which is required for the astral solar reference
    safeRWU :: flag if quality controls are applied when True
    
    returns three dataframes with RWU, RWU_without nocturnal correction, step shape NSE
    '''
    dummyd = dummyd.tz_localize(tz)
    dummyc = dummyd.columns
    
    #first column
    dummz = fRWU(dummyd[dummyc[0]])
    if safeRWU:
        dummz.loc[dummz.step_control<1100,'rwu'] = np.nan #refuse values based on too much night increase and no day decrease
        dummz.loc[dummz.rwu<0.,'rwu'] = np.nan #refuse values less than zero
    dummx = dummz.rwu
    
    if safeRWU:
        dummz.loc[dummz.step_control<1100,'rwu_nonight'] = np.nan #refuse values based on too much night increase and no day decrease
        dummz.loc[dummz.rwu_nonight<0.,'rwu_nonight'] = np.nan #refuse values less than zero
    dummy = dummz.rwu_nonight
    dummc = dummz.eval_nse
    
    #process further columns
    for i in dummyc[1:]:
        dummz = fRWU(dummyd[i])
        if safeRWU:
            dummz.loc[dummz.step_control<1100,'rwu'] = np.nan #refuse values based on too much night increase and no day decrease
            dummz.loc[dummz.rwu<0.,'rwu'] = np.nan #refuse values less than zero
        dummx = pd.concat([dummx,dummz.rwu],axis=1)
        
        if safeRWU:
            dummz.loc[dummz.step_control<1100,'rwu_nonight'] = np.nan #refuse values based on too much night increase and no day decrease
            dummz.loc[dummz.rwu_nonight<0.,'rwu_nonight'] = np.nan #refuse values less than zero
        dummy = pd.concat([dummy,dummz.rwu_nonight], axis=1)
    
        dummc = pd.concat([dummc,dummz.eval_nse],axis=1)

    dummx.columns = dummyd.columns
    dummy.columns = dummyd.columns
    dummc.columns = dummyd.columns
    return [dummx, dummy, dummc]