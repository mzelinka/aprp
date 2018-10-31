#!/usr/bin/env cdat
"""
This script demonstrates how use the APRP method to compute TOA SW radiation anomalies due to individual components
for a short (2-year) period of MPI-ESM-LR using the difference between sstClimAerosol and sstClim runs.

One should difference longer periods for more robust results -- these are just for demonstrative purposes

Variables used in this script (standard CMIP parlance):
    clt,rsdt,rsut,rsutcs,rsds,rsus,rsdscs,rsuscs

This script written by Mark Zelinka (zelinka1@llnl.gov) on 31 October 2018

Reference:
Taylor, K. E. et al. (2007), Estimating shortwave radiative forcing and response in 
    climate models, J. Clim., 20(11), 2530-2543, doi:10.1175/JCLI4143.1.  
"""
 
#IMPORT STUFF:
#=====================
import cdms2 as cdms
import cdutil
import MV2 as MV
import numpy as np
import pylab as pl
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap

 
###########################################################################
# HELPFUL FUNCTIONS FOLLOW
###########################################################################
###########################################################################
def albedo(c,a_clr,a_oc,mu_clr,mu_cld,ga_clr,ga_cld):

    """
    This function is called by APRP
    Equation numbers refer to Taylor et al. (2007)
    """

    mu_oc=mu_clr*mu_cld # Eq. 14
    ga_oc=1-(1-ga_clr)*(1-ga_cld) # Eq. 13
    A_clr=(mu_clr*ga_clr) + ((mu_clr*a_clr*(1-ga_clr)**2)/(1-(a_clr*ga_clr))) # Eq. 7
    A_oc= (mu_oc*ga_oc)   + ((mu_oc*a_oc*(1-ga_oc)**2)/(1-(a_oc*ga_oc))) # Eq. 7
    A=(1-c)*A_clr + c*A_oc # Eq. 15

    return A 
    
###########################################################################
def parameters(SWupsfccs,SWdnsfccs,SWdn,SWupcs,SWupsfcoc,SWdnsfcoc,SWupoc):

    """
    This function is called by APRP
    Equation numbers refer to Taylor et al. (2007)
    """
    
    #calculating clear sky para
    a_clr=SWupsfccs/SWdnsfccs # albedo
    Q=SWdnsfccs/SWdn # ratio of incident sfc flux to TOA insolation
    mu_clr=SWupcs/SWdn + Q*(1-a_clr) # Eq. 9
    ga_clr=(mu_clr-Q)/(mu_clr-a_clr*Q) # Eq. 10

    # calculating overcast para
    a_oc=SWupsfcoc/SWdnsfcoc # albedo
    Q=SWdnsfcoc/SWdn # ratio of incident sfc flux to TOA insolation
    mu_oc=SWupoc/SWdn + Q*(1-a_oc) # Eq. 9
    ga_oc=(mu_oc-Q)/(mu_oc-a_oc*Q) # Eq. 10

    # calculating cloud parameters through clear sky and overcast parameters
    mu_cld=mu_oc/mu_clr  # Eq. 14 sometimes this is greater than 1??
    ga_cld=(ga_oc-1)/(1-ga_clr)+1  # Eq. 13

    return (a_clr,mu_clr,ga_clr,a_oc,mu_cld,ga_cld)    

###########################################################################
def APRP(ctl_clt,ctl_rsdt,ctl_rsut,ctl_rsutcs,ctl_rsds,ctl_rsus,ctl_rsdscs,ctl_rsuscs,\
         pert_clt,pert_rsdt,pert_rsut,pert_rsutcs,pert_rsds,pert_rsus,pert_rsdscs,pert_rsuscs):
    """
    Perform APRP calculations 
    
    Reference: Taylor, K. E. et al. (2007), Estimating shortwave radiative forcing and response in 
    climate models, J. Clim., 20(11), 2530-2543, doi:10.1175/JCLI4143.1.
    
    input: total cloud cover and SW radiative fluxes at TOA and SFC for clear- and all-sky conditions
           -standard CMIP nomenclature, with prefixes that refer to the control (ctl_*) and perturbed climate (pert_*) 
    
    output: TOA SW anomalies due to changes in:
            -surface albedo (for all-, clear-, and overcast-sky conditions)
            -clouds (total change and contributions from changing cloud cover, scattering, and absorption)
            -non-cloud atmosphere (e.g., from changes in water vapor, aerosols, ozone)
            
    functions called: parameters() and albedo()
    """
    
    # Make sure the cld fractions are expressed as fraction and not percent
    if MV.max(ctl_clt)>1.:
        ctl_clt=ctl_clt/100.
    if MV.max(pert_clt)>1.:
        pert_clt=pert_clt/100.

    # Derive overcast conditions
    ctl_rsut_cld=(1/ctl_clt)*(ctl_rsut-(1-ctl_clt)*ctl_rsutcs)
    ctl_rsds_cld=(1/ctl_clt)*(ctl_rsds-(1-ctl_clt)*ctl_rsdscs)
    ctl_rsus_cld=(1/ctl_clt)*(ctl_rsus-(1-ctl_clt)*ctl_rsuscs)
    pert_rsut_cld=(1/pert_clt)*(pert_rsut-(1-pert_clt)*pert_rsutcs)
    pert_rsds_cld=(1/pert_clt)*(pert_rsds-(1-pert_clt)*pert_rsdscs)
    pert_rsus_cld=(1/pert_clt)*(pert_rsus-(1-pert_clt)*pert_rsuscs)

    # Mask these where values are unphysical   
    ctl_rsds_cld=MV.masked_where(ctl_rsds_cld > ctl_rsds,ctl_rsds_cld)
    ctl_rsus_cld=MV.masked_where(ctl_rsus_cld > ctl_rsus,ctl_rsus_cld)
    ctl_rsut_cld=MV.masked_where(ctl_rsut_cld < 0,ctl_rsut_cld)
    ctl_rsds_cld=MV.masked_where(ctl_rsds_cld < 0,ctl_rsds_cld)
    ctl_rsus_cld=MV.masked_where(ctl_rsus_cld < 0,ctl_rsus_cld)

    pert_rsds_cld=MV.masked_where(pert_rsds_cld > pert_rsds,pert_rsds_cld)
    pert_rsus_cld=MV.masked_where(pert_rsus_cld > pert_rsus,pert_rsus_cld)
    pert_rsut_cld=MV.masked_where(pert_rsut_cld < 0,pert_rsut_cld)
    pert_rsds_cld=MV.masked_where(pert_rsds_cld < 0,pert_rsds_cld)
    pert_rsus_cld=MV.masked_where(pert_rsus_cld < 0,pert_rsus_cld)    

    ## NOW THE FORMAL APRP CALCULATIONS:
    a_clr1,mu_clr1,ga_clr1,a_oc1,mu_cld1,ga_cld1 = \
        parameters(ctl_rsuscs,ctl_rsdscs,ctl_rsdt,ctl_rsutcs,ctl_rsus_cld,ctl_rsds_cld,ctl_rsut_cld) # control
    a_clr2,mu_clr2,ga_clr2,a_oc2,mu_cld2,ga_cld2 = \
        parameters(pert_rsuscs,pert_rsdscs,pert_rsdt,pert_rsutcs,pert_rsus_cld,pert_rsds_cld,pert_rsut_cld) # perturbed    

    ## Taylor et al. (2007) Eqn. 12b:
    A_1=albedo(ctl_clt,a_clr1,a_oc1,mu_clr1,mu_cld1,ga_clr1,ga_cld1)
    A_2=albedo(pert_clt,a_clr2,a_oc2,mu_clr2,mu_cld2,ga_clr2,ga_cld2)

    dA_amt_cld=     0.5*(albedo(pert_clt,a_clr1,a_oc1,mu_clr1,mu_cld1,ga_clr1,ga_cld1)-A_1) + \
                0.5*(A_2-albedo(ctl_clt,a_clr2,a_oc2,mu_clr2,mu_cld2,ga_clr2,ga_cld2))
    dA_a_clr=       0.5*(albedo(ctl_clt,a_clr2,a_oc1,mu_clr1,mu_cld1,ga_clr1,ga_cld1)-A_1) + \
                0.5*(A_2-albedo(pert_clt,a_clr1,a_oc2,mu_clr2,mu_cld2,ga_clr2,ga_cld2))
    dA_a_oc=        0.5*(albedo(ctl_clt,a_clr1,a_oc2,mu_clr1,mu_cld1,ga_clr1,ga_cld1)-A_1) + \
                0.5*(A_2-albedo(pert_clt,a_clr2,a_oc1,mu_clr2,mu_cld2,ga_clr2,ga_cld2))
    dA_abs_noncld=  0.5*(albedo(ctl_clt,a_clr1,a_oc1,mu_clr2,mu_cld1,ga_clr1,ga_cld1)-A_1) + \
                0.5*(A_2-albedo(pert_clt,a_clr2,a_oc2,mu_clr1,mu_cld2,ga_clr2,ga_cld2))
    dA_abs_cld=     0.5*(albedo(ctl_clt,a_clr1,a_oc1,mu_clr1,mu_cld2,ga_clr1,ga_cld1)-A_1) + \
                0.5*(A_2-albedo(pert_clt,a_clr2,a_oc2,mu_clr2,mu_cld1,ga_clr2,ga_cld2))
    dA_scat_noncld= 0.5*(albedo(ctl_clt,a_clr1,a_oc1,mu_clr1,mu_cld1,ga_clr2,ga_cld1)-A_1) + \
                0.5*(A_2-albedo(pert_clt,a_clr2,a_oc2,mu_clr2,mu_cld2,ga_clr1,ga_cld2))
    dA_scat_cld=    0.5*(albedo(ctl_clt,a_clr1,a_oc1,mu_clr1,mu_cld1,ga_clr1,ga_cld2)-A_1) + \
                0.5*(A_2-albedo(pert_clt,a_clr2,a_oc2,mu_clr2,mu_cld2,ga_clr2,ga_cld1))

    # if the cld fraction is less than 2%, set fields to be zero
    dA_amt_cld=MV.where(ctl_clt<0.02,0.,dA_amt_cld)
    dA_amt_cld=MV.where(pert_clt<0.02,0.,dA_amt_cld)
    dA_a_oc=MV.where(ctl_clt<0.02,0.,dA_a_oc)
    dA_a_oc=MV.where(pert_clt<0.02,0.,dA_a_oc)
    dA_abs_cld=MV.where(ctl_clt<0.02,0.,dA_abs_cld)
    dA_abs_cld=MV.where(pert_clt<0.02,0.,dA_abs_cld)
    dA_scat_cld=MV.where(ctl_clt<0.02,0.,dA_scat_cld)
    dA_scat_cld=MV.where(pert_clt<0.02,0.,dA_scat_cld)

    dA_a=dA_a_clr+dA_a_oc
    dA_cld=dA_abs_cld+dA_scat_cld+dA_amt_cld
    dA_noncld=dA_abs_noncld+dA_scat_noncld

    ## TOA SW Anomalies due to Surface Albedo Anomalies
    sfc_alb=-dA_a*ctl_rsdt
    sfc_alb_clr=-dA_a_clr*ctl_rsdt   
    sfc_alb_oc=-dA_a_oc*ctl_rsdt

    ## TOA SW Anomalies due to Cloud Anomalies
    cld=-dA_cld*ctl_rsdt
    cld_amt=-dA_amt_cld*ctl_rsdt    
    cld_scat=-dA_scat_cld*ctl_rsdt
    cld_abs=-dA_abs_cld*ctl_rsdt

    ## TOA SW Anomalies due to Non-cloud Anomalies
    noncld=-dA_noncld*ctl_rsdt
    noncld_scat=-dA_scat_noncld*ctl_rsdt  
    noncld_abs=-dA_abs_noncld*ctl_rsdt

    # set fields to zero when incoming solar radiation is zero
    dA_a=MV.where(ctl_rsdt<0.1,0.,dA_a)
    dA_noncld=MV.where(ctl_rsdt<0.1,0.,dA_noncld)
    sfc_alb=MV.where(ctl_rsdt<0.1,0.,sfc_alb)
    cld=MV.where(ctl_rsdt<0.1,0.,cld)
    noncld=MV.where(ctl_rsdt<0.1,0.,noncld)
    sfc_alb_clr=MV.where(ctl_rsdt<0.1,0.,sfc_alb_clr)
    sfc_alb_oc=MV.where(ctl_rsdt<0.1,0.,sfc_alb_oc)
    cld_amt=MV.where(ctl_rsdt<0.1,0.,cld_amt)
    cld_scat=MV.where(ctl_rsdt<0.1,0.,cld_scat)
    cld_abs=MV.where(ctl_rsdt<0.1,0.,cld_abs)
    noncld_scat=MV.where(ctl_rsdt<0.1,0.,noncld_scat)
    noncld_abs=MV.where(ctl_rsdt<0.1,0.,noncld_abs)

    # give fields axes attributes
    sfc_alb.setAxisList(pert_clt.getAxisList())
    sfc_alb_clr.setAxisList(pert_clt.getAxisList())
    sfc_alb_oc.setAxisList(pert_clt.getAxisList())
    cld.setAxisList(pert_clt.getAxisList())
    cld_amt.setAxisList(pert_clt.getAxisList())
    cld_scat.setAxisList(pert_clt.getAxisList())
    cld_abs.setAxisList(pert_clt.getAxisList())
    noncld.setAxisList(pert_clt.getAxisList())
    noncld_scat.setAxisList(pert_clt.getAxisList())
    noncld_abs.setAxisList(pert_clt.getAxisList())

    return (sfc_alb,sfc_alb_clr,sfc_alb_oc,cld,cld_amt,cld_scat,cld_abs,noncld,noncld_scat,noncld_abs)
###########################################################################

# Demonstration routine follows:

# Load in data:
direc='/work/zelinka1/git/aprp/data/'
exps = ['sstClim','sstClimAerosol']
variables = ['clt','rsdt','rsut','rsutcs','rsds','rsus','rsdscs','rsuscs']
for var in variables:
    for exp in exps:
        filename = direc+var+'_Amon_MPI-ESM-LR_'+exp+'_r1i1p2_185001-185212.nc'
        f=cdms.open(filename)
        data = f(var)
        f.close()
                           
        # Compute climatological annual cycle:
        avgdata = cdutil.ANNUALCYCLE.climatology(data) #(12, LAT, LON)
        
        if exp == 'sstClim':
            exec('ctl_'+var+' = avgdata')
        else:
            exec('pert_'+var+' = avgdata')

(sfc_alb,sfc_alb_clr,sfc_alb_oc,cld,cld_amt,cld_scat,cld_abs,noncld,noncld_scat,noncld_abs) = \
    APRP(ctl_clt,ctl_rsdt,ctl_rsut,ctl_rsutcs,ctl_rsds,ctl_rsus,ctl_rsdscs,ctl_rsuscs,\
         pert_clt,pert_rsdt,pert_rsut,pert_rsutcs,pert_rsds,pert_rsus,pert_rsdscs,pert_rsuscs)
         

# Plot Maps
lons=sfc_alb.getLongitude()[:]
lats=sfc_alb.getLatitude()[:]
LON, LAT = np.meshgrid(lons,lats)
bounds = np.arange(-12,14,2)
cmap = pl.cm.RdBu_r
bounds2 = np.append(np.append(-500,bounds),500) # This is only needed for norm if colorbar is extended
norm = mpl.colors.BoundaryNorm(bounds2, cmap.N) # make sure the colors vary linearly even if the desired color boundaries are at varied intervals

# Surface albedo and non-cloud terms
fig=pl.figure(figsize=(18,12)) # this creates and increases the figure size
names = ['sfc_alb','noncld','noncld_scat','noncld_abs']#'sfc_alb_clr','sfc_alb_oc'
for n,name in enumerate(names):
    ax=pl.subplot(2,2,n+1)
    m = Basemap(projection='robin',lon_0=0)
    m.drawmapboundary(fill_color='0.3')
    exec('DATA = MV.average('+name+',0)')
    im1 = m.contourf(LON,LAT,DATA,bounds,shading='flat',cmap=cmap,norm=norm,latlon=True,extend='both')
    m.drawcoastlines(linewidth=1.5)
    avgDATA = cdutil.averager(DATA, axis='xy', weights='weighted')
    pl.title(name+' ['+str(np.round(avgDATA,3))+']',fontsize=14)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.15)
    cb = pl.colorbar(im1, cax=cax,drawedges=True,ticks=bounds)
    cb.set_label('W/m$^2$') 
pl.subplots_adjust(hspace=-0.3) 
pl.suptitle('sstClimAerosol minus sstClim\nMPI-ESM-LR',fontsize=16,y=0.85)  
pl.savefig('/work/zelinka1/figures/APRP_noncld_example_maps.png', bbox_inches='tight')         


# Cloud terms
fig=pl.figure(figsize=(18,12)) # this creates and increases the figure size
names=['cld','cld_amt','cld_scat','cld_abs']
for n,name in enumerate(names):
    ax=pl.subplot(2,2,n+1)
    m = Basemap(projection='robin',lon_0=0)
    m.drawmapboundary(fill_color='0.3')
    exec('DATA = MV.average('+name+',0)')
    im1 = m.contourf(LON,LAT,DATA,bounds,shading='flat',cmap=cmap,norm=norm,latlon=True,extend='both')
    m.drawcoastlines(linewidth=1.5)
    avgDATA = cdutil.averager(DATA, axis='xy', weights='weighted')
    pl.title(name+' ['+str(np.round(avgDATA,3))+']',fontsize=14)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.15)
    cb = pl.colorbar(im1, cax=cax,drawedges=True,ticks=bounds)
    cb.set_label('W/m$^2$')    
pl.subplots_adjust(hspace=-0.3) 
pl.suptitle('sstClimAerosol minus sstClim\nMPI-ESM-LR',fontsize=16,y=0.85)  
pl.savefig('/work/zelinka1/figures/APRP_cld_example_maps.png', bbox_inches='tight')        

pl.show()

