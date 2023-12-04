# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 11:05:09 2022

@author: xins
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 15:37:11 2022

@author: xins
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 12:36:50 2021

@author: xins
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import flopy
import time

#%% Functions
def lin_func(x, intercept, slope):
    return intercept + (x * slope)

def make_figure(gwf, head, conc, spdis, sl, t, top0, slopels, csalt=35.,
                                xbar=None, vectors=False,pc=None):

    plt.rcParams['savefig.facecolor'] = '1.0'

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 5)
    ax.set_aspect('auto')

    pxs = flopy.plot.PlotCrossSection(
        model=gwf, line={'row': 0}, ax=ax)
    pxs.plot_grid(alpha=1., lw=0.1)

    # draw a patch collection, with the water table as the top of the cell

    pc.set_array(np.ma.masked_invalid(conc).flatten())
    ax.add_collection(pc)

    cbar = plt.colorbar(pc, shrink=0.5)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel(
        'CONCENTRATION, IN GRAMS PER LITER', rotation=90)

    atoplot = np.ma.masked_where(
        head < gwf.dis.botm.array, conc)
    pxs.contour_array(atoplot, levels=np.array([.5]) * csalt,
                      linestyles=['-'], colors=['k'])

    if vectors:
        pxs.plot_vector(qx,qy,qz, head=head, color='white',
                                    kstep=1, hstep=1, pivot='mid', minlength=0.5,
                                    scale=1., width=0.001, headwidth=3., headlength=5.0,
                                    headaxislength=4.5, minshaft=0.01,
                                    zorder=10, alpha=0.50)

    xsl = (sl - ztop0) / slopetop
    xbar = x_bar[0]
    if xbar is not None:
        xsl = np.min([xsl, xbar])
        sea_poly = [[0, sl], [xsl, sl], [xsl, lin_func(xsl, ztop0, slopetop)], [
            0, ztop0], [0, sl]]
    else:
        sea_poly = [[0, sl], [xsl, sl], [0, top0], [0, sl]]

    sea_patch = matplotlib.patches.Polygon(
        sea_poly, closed=True, edgecolor='k', facecolor='red')
    ax.add_patch(sea_patch)


    ax.set_title('TIME = {:.2f} days'.format(t))
    ax.set_xlabel('DISTANCE, IN METERS')
    ax.set_ylabel('ELEVATION, IN METERS')

    plt.tight_layout()
    
    fig.savefig('{}/{}.png'.format(modeldir,model_name), dpi=300)
    plt.close(fig)
    plt.show()
    return pxs, ax

#%%
# Model parameters
all_start = time.time()
nlayrowcol = [20, 1, 200]
Lx = 1e3

ztop0 = -2
zbot0 = ztop0 - 50
recharge = 200 / 1000. / 365.  # convert mm/yr to m/yr to m/day
porosity = 0.2
nyr = 100.
time_end = nyr*365.
periods_per_year = 365.
nper = 100
#nper = np.floor(nyr*periods_per_year).astype(int)
slopetops = [.005, .007, .01]
layers_ndeep = np.hstack([np.arange(2, 20, 2)])
nmodels = np.shape(slopetops)

slopebot = 0
nstp = 12  # every 2 month
nlay, nrow, ncol = nlayrowcol
delr = float(Lx)/float(ncol)
# Use cell center for elevation
xcol_centers = np.arange(delr*0.5, delr*ncol, delr)
slr_rate = 3.6/1000/365  # m/day NOAA report, global average slr rate https://www.climate.gov/news-features/understanding-climate/climate-change-global-sea-level
#slr_rate = 0 #spinup
# %%
# package setting parameters
save_flows = True

converge_continue = False
rewet = True
wetfct = 0.5
iwetit = 1
wetdry = 1e-4
ihdwet = 1
Khs = [0.01, 0.1, 1,10]

use_uzf = True

Kv = None
cfresh = 0
csalt = 35
rho0 = 1000. #density
drhodc = 0.7

hstrt = None
cstrt = None

relax = 0.97

barrier_type = 'dam'
sl_xs = [0, 30, 50, 100, 150, 200, 300]
# %%

# Input to model
for slopetop in slopetops:  # looping through slopetop
    # ----------- Region directory information -----------
    for sl_x in sl_xs:

        for Kh in Khs:

            for ndeep in layers_ndeep:

                model_start = time.time()

                wdir = r'E:\Postdoc\EESLR\data\model'
                #wdir = r'F:\OneDrive - University of Arkansas\advising\group\XinSu\paper1\xsec_uzf'
                main_model_dir = os.path.join(wdir, 'model')
                #modelname = 'SLR_3.6V3'
                
                if use_uzf:
                    modelname = 'SLR_uzf_update_0730'
                else:
                    modelname = 'SLR_3.6V4'
                
                modelname2 = 'spinup_uzf_update_0729_check'
                model_name_fmt = '{0}_{1}_{2}'
                region_dir = os.path.join(main_model_dir, modelname)
                results_dir = os.path.join(region_dir, 'output')
                model_name = model_name_fmt.format(slopetop, sl_x, ndeep)
                model_in_dir = os.path.join(region_dir, model_name)
                model_out_dir = os.path.join(model_in_dir, 'output')
                modeldir = os.path.join(
                    wdir, modelname, barrier_type, 'Kh_{}'.format(Kh), model_name)
                spinupdir = os.path.join(
                    wdir, modelname2, barrier_type, 'Kh_{}'.format(Kh), model_name)
                pre_fig = os.path.join(results_dir, 'Figures')

                m_info_dict = {'workspace': model_in_dir,
                               'model_name': model_name}
                
                print(
                    '------------- Model {} of {} -------------'.format(slopetop, nmodels))
                print('Model: {}, Kh = {} m / day'.format(model_name, Kh))
               

                model_dict1 = {'ws': modeldir, 'name': model_name, 'nlayrowcol': nlayrowcol,
                               'Lx': Lx, 'slopetop': slopetop, 'slopebot': 0,
                               'ztop0': ztop0, 'zbot0': zbot0, 'sl_start': -1.5, 'slr_rate': slr_rate,  # m/day
                               'time_end': time_end,
                               'Kh': Kh, 'ss': 0.008,
                               'porosity': porosity, 'recharge': recharge,
                               'alpha_l': 1., 'alpha_t': 0.1,
                               'nper': nper, 'nstp': 1,
                               
                               'cnst_c': True,
                               'drn_bool': True, 'barrier': None}  # drn_bool - True: free drainage away from cell; False: fresh water ponded above top of cell to equal sea level

                # %%MODFLOW 6 flopy simulation object (sim)

                # Model names
                gwfname = 'gwf_{}'.format(model_name)
                gwtname = 'gwt_{}'.format(model_name)
                # %%build MODFLOW 6 file
                sim = flopy.mf6.MFSimulation(sim_name=modelname, version='mf6',
                                             exe_name='mf6.exe',
                                             verbosity_level=1,
                                             sim_ws=modeldir)
                if converge_continue:
                    sim.name_file.continue_ = True
                    
                
                # %%Spatial discretization
                nlay, nrow, ncol = nlayrowcol
                delr = float(Lx)/float(ncol)

                # Use cell center for elevation
                xcol_centers = np.arange(delr*0.5, delr*ncol, delr)
                ztop = (ztop0 + xcol_centers*slopetop)[None, :]  # make 2D
                # Set layer thickness by total thickness and nlay
                botz = (zbot0 + xcol_centers*slopebot)[None, :]  # make 2D
                
                thick = ztop-botz
                laythick = thick/nlay  # thickness of each layer cell as function of x
                
                
                botm = ztop[0,0]-np.dot((laythick[0,0]*np.arange(1,nlay+1))[:,None],np.ones_like(ztop))
                botm = botm[:, None, :]  # add row (y) dimension
                # Vectorize delr,delc
                delr = np.ones((ncol)) * delr
                delc = 1
                delc = np.ones((nrow)) * delc

      #%%          # Time discretization
                perlen = nper * [time_end / nper]
                tdis_rc = []
                for i in range(nper):
                    tdis_rc.append((perlen[i], nstp, 1.0))

                # create tdis package
                tdis = flopy.mf6.ModflowTdis(sim, time_units='DAYS',
                                                 nper=nper, perioddata=tdis_rc)
                # %%
                newton = False
                if newton:
                    newtonoptions = 'newton'
                    no_ptc = 'ALL'
                    complexity = 'complex'
                else:
                    newtonoptions = None
                    no_ptc = None
                    complexity = 'simple'
                gwf = flopy.mf6.ModflowGwf(
                    sim, modelname=gwfname, newtonoptions='NEWTON UNDER_RELAXATION')
                # %%iterative model solution (IMS)
                n_inout = [500, 100]
                hr_close = [1e-2, 1e-2]
                ims_defaults = {'print_option': 'ALL', 'under_relaxation': 'DBD',
                                'linear_acceleration': 'BICGSTAB', 'scaling_method': 'NONE',
                                'reordering_method': 'NONE', }
                ims_dict = ims_defaults.copy()
                ims_dict.update({'outer_dvclose': hr_close[0], 'outer_maximum': n_inout[1],
                                 'inner_dvclose': hr_close[0], 'inner_maximum': n_inout[0],
                                 'rcloserecord': hr_close[1], 'relaxation_factor': relax,
                                 'filename': '{}.ims'.format(gwfname)})
                imsgwf = flopy.mf6.ModflowIms(
                    sim, **ims_dict, complexity='COMPLEX', no_ptcrecord='ALL')
                single_matrix = False
                if single_matrix:
                    sim.register_ims_package(imsgwf, [gwfname, gwtname])
                else:
                    sim.register_ims_package(imsgwf, [gwfname])
                # %% model domain
                #barrier = None

                

                barrier = None
                sl = lin_func(nyr*periods_per_year,
                              model_dict1['sl_start'], 0)
                x_sl = (sl - ztop0) / model_dict1['slopetop']
                x_sl = x_sl + sl_x
                sl_ind = (xcol_centers <= x_sl).nonzero()[0][-1]
                sl_ind = sl_ind + 2  # shift seaward (negative) or landward (+)
                barrier = np.zeros(nlayrowcol, dtype=bool)
              
                nwide = 1
                #barrier[:ndeep, :, sl_ind+1:sl_ind+1+nwide] = True #cutoffwall
                barrier[ndeep:, :, sl_ind+1:sl_ind+1+nwide] = True #dam
                x_bar = (sl_ind+1)*delr

                idomain = np.ones((nlay, nrow, ncol))
                # need to set idomain to 0 for impermeable cells
                if barrier is not None:
                    # barrier is bool mask [nlay,nrow,ncol] with True where impermeable
                    idomain[barrier] = 0
                    # first column in top layer that is part of the barrier
                    max_ghb_ind = barrier.nonzero()[-1][0]
                else:
                    max_ghb_ind = ncol

                # Make impermeable barrier


                dis = flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol,
                                              delr=delr, delc=delc,
                                              top=ztop, botm=botm, idomain=idomain)

                # %%Initial conditions
                hdsfile = os.path.join(spinupdir,'gwf_{}.hds'.format(model_name))
                if os.path.isfile(hdsfile):
                    #hdsfile = r'E:/Postdoc/EESLR/data/model/spinup_V1/dam/Kh_0.01/0.005_0_2/gwf_0.005_0_2.hds'
                    hdobj = flopy.utils.binaryfile.HeadFile(hdsfile, precision='double')
                    head0 = hdobj.get_data(idx = -1)
                    head0 [head0 > 1e4] = 0
                    
                    
                    if hstrt is None:
                        hstrt = head0  # start at initial sea level if not provided
                else:
                    hstrt = 0
                
                sl_start = model_dict1['sl_start']
                ic = flopy.mf6.ModflowGwfic(gwf, strt=hstrt)
                # %%
                # Node property flow
                ss = 0.008
                hhformulation_rhs = False

                npf = flopy.mf6.ModflowGwfnpf(gwf, xt3doptions=False,
                                              save_flows=save_flows,
                                              save_specific_discharge=True,
                                              icelltype=1,
                                              k=Kh, wetdry=1e-4, rewet_record=None,
                                              )
               
                #
                sto = flopy.mf6.ModflowGwfsto(
                    gwf, iconvert=1, sy=porosity, ss=ss, steady_state = None, transient={0:True})

                # NEW SEAWAT package
                buy_pd = [0, drhodc, 0, gwtname, 'SALINITY']  #drhodc is the real value that defines the slope o the density-concentration line 
                buy = flopy.mf6.ModflowGwfbuy(gwf, nrhospecies=len(buy_pd), denseref=rho0,  # concentration=cfresh,
                                              packagedata=buy_pd, hhformulation_rhs=hhformulation_rhs,)

                # %% GHB package (general head boundary package)
                # 2 GHB package, 1st seaward boundary, 2nd model top boundary

                # 1st GHB for seaward boundary cells
                delt = time_end / nper
                # uniform layer thickness by col only
                delv_top = ztop[0, 0]-botm[0, 0, 0]
                cond1 = delc[0] * int(delv_top) * Kh/(0.5*delr[0])
                timeseries1 = [(0.0, sl_start), (time_end+delt,
                                                 (sl_start + (time_end + delt)*slr_rate))]
               
                dens_salt = rho0 + drhodc * csalt

                ghblist1 = []
                for k in range(nlay):
                    ghblist1.append(
                        [(k, 0, 0), 'sealevel', cond1, csalt, dens_salt])
                ghb1 = flopy.mf6.ModflowGwfghb(gwf,
                                               stress_period_data=ghblist1,
                                               print_input=True,
                                               print_flows=True,
                                               save_flows=save_flows,
                                               pname='GHB-1',
                                               auxiliary=[
                                                   'SALINITY', 'DENSITY'],
                                               timeseries=timeseries1,
                                               filename='{}.ghb'.format(gwfname))
                ghb1.ts.time_series_namerecord = 'sealevel'
                ghb1.ts.interpolation_methodrecord = 'linearend'

                # %% 2nd top boundary condition of the study domain
                # 
                ghbdict2 = {}
                cncdict = {}
                uzfdict = {}
                if Kv is None:
                    Kv = Kh
                cond2 = delc[0] * delr[0] * Kv / (0.5 * delv_top)

                x_vals0 = np.cumsum(delr)-delr[0]  # cell left edge
                x_vals1 = x_vals0.copy() + delr[0]  # cell right edge
                tp0 = ztop0 + slopetop * x_vals0
                tp1 = ztop0 + slopetop * x_vals1
                
                    
                drn_bool = True
                uzf_pd1dict = {}
                
                uzf_pd2dict = {} 
                
                for iper, t in enumerate(np.arange(0, time_end + delt, delt)):

                    sl = sl_start + t * slr_rate
                    print(t, sl)
                    xsl = (sl - ztop0) / slopetop
                    ghblist2 = []
                    #rchlist = []
                    cnclist = []
                    #drnlist = []
                    uzf_pd1 = []  #packagedata
                    uzf_pd2 = []  #perioddata
                    
                    if iper == 0:
                        uzflist = []
                        nuzfcells = []
                    # x = delr[0] / 2.
                    for j in range(ncol):
                        if idomain[0, 0, j] != 0:  # not inactive cells

                            if sl < tp0[j] and sl < tp1[j]:
                                # sea level is lower than the top of the model cell
                                fcond = 0.  # only fresh

                            elif sl > tp0[j] and sl > tp1[j] and j > xsl:
                                # sea level higher than both edges of the model cell
                                fcond = 1.  # only saline

                            elif j < xsl:
                                # sea level intersects sloping top of cell so calculate
                                # fraction to use for conductance multiplier
                                fcond = (sl - tp0[j]) / \
                                    (tp1[j] - tp0[j])  # mixture
                            else:
                                fcond = 0.  # inside barrier, only fresh recharge

                            if fcond > 0.:
                                ghblist2.append(
                                    [(0, 0, j), sl, fcond * cond2, csalt, dens_salt, sl])
                            frch = 1. - fcond
                            temp_drn_elev = tp0[j]+(fcond*(tp1[j]-tp0[j]))
                            if frch > 0.:
                               # uzflist.append(
                                   # [(0, 0, j), recharge * frch, 0.])
                                if iper==0:
                                    nuzfcells.append ( (j)) # which columns in stress period are uzf
                                    
                                    
                               
                            else:
                               
                                if temp_drn_elev < sl:
                                        # flooded with fresh water to sea level
                                        ghblist2.append(
                                            [(0, 0, j), sl, cond2, cfresh, rho0, temp_drn_elev])

                            if fcond == 1.:
                                cnclist.append([(0, 0, j), csalt])
                            
                            # x += delr[0]
                    if len(ghblist2) > 0:
                        ghbdict2[iper] = ghblist2
                    
                        
                    #if len(rchlist) > 0:
                        #rchdict[iper] = rchlist
                    #if len(drnlist) > 0:
                        #drndict[iper] = drnlist
                    #if len(uzflist)  > 0:
                       # uzfdict[iper] = uzflist
                        
                    if len(cnclist) > 0:
                        cncdict[iper] = cnclist
                        

                
                if len(ghbdict2) > 0:
                    ghb2 = flopy.mf6.ModflowGwfghb(gwf,
                                                   stress_period_data=ghbdict2,
                                                   print_input=True,
                                                   print_flows=True,
                                                   save_flows=save_flows,
                                                   auxiliary=[
                                                       'SALINITY', 'DENSITY', 'ELEVATION'],
                                                   pname='GHB-2',
                                                   filename='{}2.ghb'.format(gwfname))

               
                drndict = [] # uncomment to turn off drains
                        
                if len(nuzfcells) > 0:
                    surfdep = delr[0] *slopetop 
                    icount = -1             
                    for icell,icol in enumerate(nuzfcells):
                        nlay_uzf = 1 # only to top layer
                        for ilay in np.arange(nlay_uzf): # first nlay_uzf layers could be dry
                            icount+=1
                            if ilay == 0:
                                if nlay_uzf ==1:
                                    ivcont = 0
                                else:
                                    ivcont = icount+1
                                uzf_pd1.append ((icount,(ilay,0,icol),1,ivcont,surfdep,Kv,porosity/2,porosity,porosity/2,4) )
                                uzf_pd2.append (( icount, recharge, 0, 0, 0,0,0,0,cfresh, rho0,))
                            elif ilay<nlay_uzf-1:
                                uzf_pd1.append ((icount,(ilay,0,icol),0,icount+1,surfdep,Kv,porosity/2,porosity,porosity/2,4) )
                                uzf_pd2.append (( icount, 0, 0, 0, 0,0,0,0,cfresh, rho0,))
                            else: # bottom layer
                                uzf_pd1.append ((icount,(ilay,0,icol),0,0,surfdep,Kv,porosity/2,porosity,porosity/2,4) )
                                uzf_pd2.append ((icount, 0, 0, 0, 0,0,0,0,cfresh, rho0,))
                    
                    if len(uzf_pd1)  > 0:
                        uzf_pd1dict[0] = uzf_pd1
                    if len(uzf_pd2)  >0:
                        uzf_pd2dict[0] = uzf_pd2
               
                if len(uzf_pd1) > 0 :
                    uzf = flopy.mf6.ModflowGwfuzf(gwf,
                                                  simulate_et = False,
                                                  simulate_gwseep = True,
                                                  packagedata = uzf_pd1,
                                                  perioddata = uzf_pd2,
                                                  pname='UZF-1',
                                                  auxiliary=[
                                                      'SALINITY', 'DENSITY'],
                                                  nuzfcells = icount+1,#len(nuzfcells),
                                                  print_input=True,
                                                  print_flows=True,
                                                  save_flows = save_flows,
                                                  nwavesets = 150,
                                                  ntrailwaves = 50,
                                                  )    
                    
                  
                
                if len(drndict) > 0:
                    drn = flopy.mf6.ModflowGwfdrn(gwf,
                                                  stress_period_data=drndict,
                                                  pname='DRN-1',
                                                  auxiliary=[
                                                      'SALINITY', 'ELEVATION'],
                                                  print_input=True,
                                                  print_flows=True,
                                                  save_flows=save_flows)
                    
                     
                # %%
                # output control
                saverecord = {0: [('HEAD', 'ALL'), ('BUDGET', 'ALL')],
                              1: [('HEAD', 'ALL'), ('BUDGET', 'ALL')],
                              nper - 1: [('HEAD', 'ALL'), ('BUDGET', 'ALL')]}
                printrecord = {0: [('HEAD', 'ALL'), ('BUDGET', 'LAST')],
                               1: None,
                               nper - 1: [('HEAD', 'ALL'), ('BUDGET', 'LAST')]}
                oc = flopy.mf6.ModflowGwfoc(gwf,
                                            budget_filerecord='{}.cbc'.format(
                                                gwfname),
                                            head_filerecord='{}.hds'.format(
                                                gwfname),
                                            headprintrecord=[
                                                ('COLUMNS', 10, 'WIDTH', 15,
                                                 'DIGITS', 6, 'GENERAL')],
                                            saverecord=saverecord,
                                            printrecord=printrecord)

                # %%create transport model
                gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname,
                                        model_nam_file='{}.nam'.format(gwtname))
                #gwt.write()
                
                # %%
                if not single_matrix:
                    imsgwt = flopy.mf6.ModflowIms(sim, print_option='ALL', complexity=complexity, no_ptcrecord='FIRST',
                                                  outer_dvclose=hr_close[0],
                                                  outer_maximum=n_inout[1],
                                                  under_relaxation='NONE',
                                                  inner_maximum=n_inout[0],
                                                  inner_dvclose=hr_close[0], rcloserecord=hr_close[1],
                                                  linear_acceleration='BICGSTAB',
                                                  scaling_method='NONE',
                                                  reordering_method='NONE',
                                                  relaxation_factor=relax,
                                                  filename='{}.ims'.format(gwtname))
                    sim.register_ims_package(imsgwt, [gwt.name])

                dis = flopy.mf6.ModflowGwtdis(gwt, nlay=nlay, nrow=nrow, ncol=ncol,
                                              delr=delr, delc=delc,
                                              top=ztop, botm=botm, idomain=idomain)

                # %% initial conditions
                cbcfile = os.path.join(spinupdir,'gwt_{}.ucn'.format(model_name))
                if os.path.isfile(cbcfile):
                    #cbcfile = r'E:/Postdoc/EESLR/data/model/spinup_V1/dam/Kh_0.01/0.005_0_2/gwt_0.005_0_2.ucn'
                    cnobj = flopy.utils.HeadFile(cbcfile, precision='double', text='CONCENTRATION')
                    conc0 = cnobj.get_data(idx= -1)
                    conc0[conc0 > 1e4]=0
                    if cstrt is None:
                        cstrt = conc0
                else:
                    cstrt = cfresh
                ic = flopy.mf6.ModflowGwtic(gwt, strt=cstrt)

                # %%
                # advection
                adv = flopy.mf6.ModflowGwtadv(gwt, scheme='UPSTREAM')

                # dispersion
                dsp = flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=None, xt3d_rhs=None, diffc=0.,
                                              alh=1, ath1=0.1)

                # mass storage and transfer
                mst = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)

                # sources
                sourcerecarray = [('GHB-1', 'AUX', 'SALINITY'),
                                  #('RCH-1', 'AUX', 'SALINITY'),
                                  ('UZF-1', 'AUX', 'SALINITY'),]
                if len(ghbdict2) > 0:
                    sourcerecarray.append(('GHB-2', 'AUX', 'SALINITY'))

                ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray,
                                              filename='{}.ssm'.format(gwtname))

                fmi = flopy.mf6.ModflowGwtfmi(
                    gwt, flow_imbalance_correction=True)

                # constant concentration
                cnst_c = False
                if cnst_c:
                    cnclist = []
                    for k in range(1, nlay):
                        cnclist.append([(k, 0, 0), csalt])
                    cnc1 = flopy.mf6.ModflowGwtcnc(
                        gwt, stress_period_data=cnclist)

                    if len(cncdict) > 0:
                        cnc2 = flopy.mf6.ModflowGwtcnc(gwt, stress_period_data=cncdict,
                                                       pname='CNC-2',
                                                       filename='{}2.cnc'.format(gwtname))
                # %%
                # output control
                saverecord = {0: [('CONCENTRATION', 'ALL'), ('BUDGET', 'LAST')],
                              1: [('CONCENTRATION', 'ALL')],
                              nper - 1: [('CONCENTRATION', 'ALL'), ('BUDGET', 'LAST')]}
                printrecord = {0: [('CONCENTRATION', 'ALL'), ('BUDGET', 'LAST')],
                               1: None,
                               nper - 1: [('CONCENTRATION', 'ALL'), ('BUDGET', 'LAST')]}
                oc = flopy.mf6.ModflowGwtoc(gwt,
                                            budget_filerecord='{}.cbc'.format(
                                                gwtname),
                                            concentration_filerecord='{}.ucn'.format(
                                                gwtname),
                                            concentrationprintrecord=[
                                                ('COLUMNS', 10, 'WIDTH', 15,
                                                 'DIGITS', 6, 'GENERAL')],
                                            saverecord=saverecord,
                                            printrecord=printrecord)
                # %%
                # GWF GWT exchange
                gwfgwt = flopy.mf6.ModflowGwfgwt(sim, exgtype='GWF6-GWT6',
                                                 exgmnamea=gwfname, exgmnameb=gwtname,
                                                 filename='{}.gwfgwt'.format(modelname))
                # %%
                sim.write_simulation()
                # %%
                success, buff = sim.run_simulation()
                if not success:
                    raise Exception("MODFLOW did not terminate normally.")

                gwf = sim.get_model('gwf_{}'.format(modelname))
                istp = -1

                # %% Getting results
                modelnames = list(sim.model_names)
                gwfname = modelnames[0]
                gwtname = modelnames[1]
                gwf = sim.get_model(gwfname)
                gwt = sim.get_model(gwtname)

                fname_h = os.path.join(modeldir, gwfname + '.hds')
                hdobj = flopy.utils.HeadFile(fname_h, precision='double')
                head = hdobj.get_alldata()
                head[head > 3e2] = np.nan

                fname_c = os.path.join(modeldir, gwtname + '.ucn')
                cnobj = flopy.utils.HeadFile(
                    fname_c, precision='double', text='CONCENTRATION')
                conc = cnobj.get_alldata()
                conc[conc > 1e2] = np.nan
                times = cnobj.get_times()

                fname_bu = os.path.join(modeldir, gwfname + '.cbc')
                budobj = flopy.utils.CellBudgetFile(
                    fname_bu, precision='double')
                spdis = budobj.get_data(text='DATA-SPDIS')[-1]
                qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf,head = head[-1])
    # %%  Plotting
                ec = None
                fc = None
                zorder = 2
                xv, yv, zv = gwf.modelgrid.xyzvertices
                patches = []
                for k in range(gwf.modelgrid.nlay):
                    for j in range(gwf.modelgrid.ncol):
                        x0 = xv[0, j]
                        x1 = xv[0, j + 1]
                        y0 = zv[k, 0, j]
                        y1 = zv[k + 1, 0, j]
                        if head is not None:
                            y0 = min(y0, head[-1][k, 0, j])
                            y0 = max(y0, gwf.modelgrid.botm[k, 0, j])
                        poly = [[x0, y0], [x1, y0], [
                            x1, y1], [x0, y1], [x0, y0]]
                        # print(poly)
                        patch = matplotlib.patches.Polygon(
                            poly, closed=True, edgecolor=ec, facecolor=fc)
                        patches.append(patch)
                pc = matplotlib.collections.PatchCollection(
                    patches, zorder=zorder, match_original=True, cmap='jet')



                # %%plotting
                x_bar = (sl_ind+1)*delr
                px, axf = make_figure(gwf, head[istp], conc[istp], spdis,
                                      sl, times[istp], ztop0, model_dict1['slopetop'], xbar=x_bar, vectors=True,pc=pc)
# %%
print('Total time for modeling: {0:4.1f} min'.format(
    (time.time()-all_start)/60.))

