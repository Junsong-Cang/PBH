try:
    from main import *
except:
    from src.main import *

def HyRec_Path(idx = 0):
    r = main_path + 'src/HyRec/'
    return r

def Get_Lbh_Spectra(Lbh, E):
    '''
    dEtot/dEdt
    --inputs--
    E : energy in eV
    Lbh : bolometric luminosity in eV/s
    --output--
    dEtot/dEdt in 1/s
    '''

    # Emin = 2e3
    # Emax = 1e4
    # E1 = 1e5
    # E = np.linspace(Emin, Emax, 10000)
    # S1 = ((E/Emin)**(-0.5)) * np.exp(-E/E1)
    # Normalisation = np.trapz(x = E, y = S1)
    # print(Normalisation)

    Normalisation = 4681.372338022844
    Emin = 2e3
    Emax = 1e4
    E1 = 1e5
    if E<Emin or E>Emax:
        r = 0
    else:
        Shape = ((E/Emin)**(-0.5)) * np.exp(-E/E1)
        L0 = Lbh/Normalisation
        r = L0 * Shape
        # print(Lbh, '--', L0, '--', Shape)
    return r

def Get_I(Lbh, E, m0, fbh0, z):
    '''
    output : dN/dEdVdt in 1/(Mpc^3 s eV)
    '''
    OmC = cosmo['OmC']
    h = cosmo['h']
    RhoCr = 2.775e11 * h**2 # msun/Mpc^3
    RhoC = OmC * RhoCr * (1+z)**3
    nbh = fbh0 * RhoC/m0 # Mpc^-3
    dNdEdt = Get_Lbh_Spectra(Lbh, E)/E

    r = nbh * dNdEdt
    return r
    
def Get_dEdVdt_Kernel(m0, fbh0, Lbh_vec, z_vec, channel):
    '''
    Get dEdVdt_dep as a function of z
    ----inputs----
    Lbh_vec : Lbh in eV
    channel : deposition channel [1, 3, 4] = [HIon, LyA, Heat]
    '''
    
    Mpc_cm = Constants['Mpc']*100
    TransferData = Interpolation_Table['TransferData']
    cid = channel - 1
    T_axis = TransferData['T']
    z_axis = TransferData['z_axis']
    Ek_axis = TransferData['Ek_axis']
    Hubble_axis = TransferData['Hubble_axis']
    nz = len(z_axis)
    ne = len(Ek_axis)
    
    z_vec_ = z_vec[::-1]
    Lbh_vec_ = Lbh_vec[::-1]

    Lbh_axis = np.interp(x = np.log(z_axis), xp = np.log(z_vec_), fp = Lbh_vec_, left = 0, right = 0)

    OmC = cosmo['OmC']
    h = cosmo['h']
    RhoCr = 2.775e11 * h**2
    RhoC = RhoCr * OmC *(1+z_axis)**3 # mcun/Mpc^3
    nbh = fbh0 * RhoC/m0 # Mpc^3
    
    dEdVdt_tot = nbh * Lbh_axis/ Mpc_cm**3
    I_Factor = np.empty((nz, ne))

    for z_idx in np.arange(0, nz):
        for e_idx in np.arange(0, ne):
            z = z_axis[z_idx]
            E = Ek_axis[e_idx]
            I_Factor[z_idx, e_idx] = Get_I(
                Lbh = Lbh_axis[z_idx],
                E = E,
                m0 = m0,
                fbh0 = fbh0,
                z = z)
    
    Kernel_E = np.linspace(0, 1, ne)
    Kernel_Z = np.linspace(0, 1, nz)
    r = np.linspace(0, 1, nz)
    
    for zd_idx in np.arange(0, nz):
        for zi_idx in np.arange(0, nz):
            for e_idx in np.arange(0, ne):
                
                E = Ek_axis[e_idx]
                T = T_axis[cid, zd_idx, e_idx, zi_idx]
                I = I_Factor[zi_idx, e_idx]
                Kernel_E[e_idx] = E * T * I
            
            Int_E = np.trapz(x = Ek_axis, y = Kernel_E)
            Hp = Hubble_axis[zi_idx]
            zp = z_axis[zi_idx]
            
            Kernel_Z[zi_idx] = Int_E/(Hp * (1+zp)**4)
            if z_axis[zi_idx] < z_axis[zd_idx]:
                Kernel_Z[zi_idx] = 0
        
        H = Hubble_axis[zd_idx]
        z = z_axis[zd_idx]
        Int_zi = np.trapz(x = z_axis, y = Kernel_Z)
        r[zd_idx] = H * (1+z)**3 * Int_zi
    
    dEdVdt_dep = r /Mpc_cm**3 # convert to eV/Mpc^3/s unit
    fc = dEdVdt_dep/dEdVdt_tot # Deposition efficiency
    r = {'z': z_axis, 'dEdVdt_dep': dEdVdt_dep, 'fc' : fc}
    
    return r

def Print_dEdVdt(dEdVdt_dataset, File):
    HIon = dEdVdt_dataset['dEdVdt_HIon']
    LyA = dEdVdt_dataset['dEdVdt_LyA']
    Heat = dEdVdt_dataset['dEdVdt_Heat']
    nz = len(LyA)
    F=open(File,'w')
    for idx in np.arange(0,nz):
        print("{0:.8E}".format(HIon[idx]), file=F)
    for idx in np.arange(0,nz):
        print("{0:.8E}".format(LyA[idx]), file=F)
    for idx in np.arange(0,nz):
        print("{0:.8E}".format(Heat[idx]), file=F)
    F.close()

def Get_dEdVdt_Lite(
        fbh0 = 1e-5,
        m0 = 1e3,
        zmin = 8,
        zmax = 1000,
        nz = 10000,
        Use_Halo_Boost = False,
        Use_EoR = False,
        Use_OmB_FeedBack = True,
        Fix_mdot = False,
        Fixed_mdot = 10,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Show_Extrap_MSG = False,
        Ignore_Evo = False
    ):
    '''
    Get dEdVdt_dep and fc at all 3 channels, dEdVdt_dep are in eV/s/cm^3
    also get dEdVdt_inj at entier evo history
    '''
    r0 =  Get_Evo_Lite(
            fbh0 = fbh0,
            m0 = m0,
            OmB0 = cosmo['OmB'],
            zmin = zmin,
            zmax = zmax,
            nz = nz,
            Use_Halo_Boost = Use_Halo_Boost,
            Use_EoR = Use_EoR,
            Use_OmB_FeedBack = Use_OmB_FeedBack,
            Fix_mdot = Fix_mdot,
            Fixed_mdot = Fixed_mdot,
            Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
            Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
            Use_mdot_Interp = Use_mdot_Interp,
            Show_Extrap_MSG = Show_Extrap_MSG,
            Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
            Ignore_Evo = Ignore_Evo)
    
    r1 = Get_dEdVdt_Kernel(
        m0 = m0,
        fbh0 = fbh0,
        Lbh_vec = r0['Lbh_vec_eV'],
        z_vec = r0['z'],
        channel = 1)
    r3 = Get_dEdVdt_Kernel(
        m0 = m0,
        fbh0 = fbh0,
        Lbh_vec = r0['Lbh_vec_eV'],
        z_vec = r0['z'],
        channel = 3)
    r4 = Get_dEdVdt_Kernel(
        m0 = m0,
        fbh0 = fbh0,
        Lbh_vec = r0['Lbh_vec_eV'],
        z_vec = r0['z'],
        channel = 4)
    
    z = r1['z']
    dEdVdt_HIon = r1['dEdVdt_dep']
    dEdVdt_LyA = r3['dEdVdt_dep']
    dEdVdt_Heat = r4['dEdVdt_dep']
    fc_HIon = r1['fc']
    fc_LyA = r3['fc']
    fc_Heat = r4['fc']

    # Let's also print OmB to check consistency with background HyRec cosmology
    z_vec = r0['z'][::-1]
    OmB_ratio_vec = r0['OmB_ratio'][::-1]
    OmB_ratio = np.interp(x = z, xp = z_vec, fp = OmB_ratio_vec)

    r = {
        'z' : z,
        'dEdVdt_HIon' : dEdVdt_HIon,
        'dEdVdt_LyA' : dEdVdt_LyA,
        'dEdVdt_Heat' : dEdVdt_Heat,
        'fc_HIon' : fc_HIon,
        'fc_LyA' : fc_LyA,
        'fc_Heat' : fc_Heat,
        'OmB_ratio' : OmB_ratio
        }
    
    return r

def Get_dEdVdt(
        fbh0 = 1e-5,
        mc = 1e3,
        sbh = 0.5,
        Use_Ps_domain = False,
        PsA = 0.1,
        kc = 1e5,
        SigmaPs = 0.5,
        zmin = 10,
        zmax = 500,
        nz = 500,
        nm = 500,
        DeltaC = 0.45,
        Use_EoR = False,
        dmdm0_method = 0,
        Use_Halo_Boost = False,
        sbh_range = [-3, 5],
        show_status = True,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Break_OmB_Lost_Ratio = 1.1,
        Fix_mdot = False,
        Fixed_mdot = 10,
        Show_Extrap_MSG = False,
        Ignore_Evo = False
    ):
    '''
    Get dEdVdt_dep and fc. If not Use_Ps_domain and sbh < 0, will use monochromatic results for m0 = mc, fbh0 = fbh0
    '''
    if show_status:
        t1 = TimeNow()
        
    if not Use_Ps_domain and sbh < 0:
        r= Get_dEdVdt_Lite(
            fbh0 = fbh0,
            m0 = mc,
            zmin = zmin,
            zmax = zmax,
            nz = nz,
            Use_Halo_Boost = Use_Halo_Boost,
            Use_EoR = Use_EoR,
            Use_OmB_FeedBack = True,
            Fix_mdot = Fix_mdot,
            Fixed_mdot = Fixed_mdot,
            Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
            Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
            Use_mdot_Interp = Use_mdot_Interp,
            Show_Extrap_MSG = Show_Extrap_MSG,
            Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
            Ignore_Evo = Ignore_Evo)
        if show_status:
            print('Time Used for Get_dEdVdt: ', time.time() - t1)
        return r
    
    Evo =  Get_Evo(
        mc = mc,
        sbh = sbh,
        fbh0 = fbh0,
        Use_Ps_domain = Use_Ps_domain,
        PsA = PsA,
        kc = kc,
        SigmaPs = SigmaPs,
        DeltaC = DeltaC,
        OmB0 = cosmo['OmB'],
        Use_EoR = Use_EoR,
        dmdm0_method = dmdm0_method,
        Use_Halo_Boost = Use_Halo_Boost,
        zmin = zmin,
        zmax = zmax,
        Use_HZ = True,
        nz = nz,
        nm = nm,
        sbh_range = sbh_range,
        show_status = show_status,
        Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
        Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
        Use_mdot_Interp = Use_mdot_Interp,
        Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
        Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
        Break_OmB_Lost_Ratio = Break_OmB_Lost_Ratio,
        Fix_mdot = Fix_mdot,
        Show_Extrap_MSG = Show_Extrap_MSG,
        Fixed_mdot = Fixed_mdot,
        Ignore_Evo = Ignore_Evo)
    
    if show_status:
        print('Evolution hisotry acquired, computing deposition rate')
    MF0 = Evo['MF'][0,:]
    m0_vec = Evo['m_vec'][0,:]
    f0 = Evo['fbh_vec'][0] # could be different from input if UsePs
    Lbh = Evo['Lbh_vec_eV'] # index : z, m
    z_vec = Evo['z_vec']
    
    z_axis = Interpolation_Table['TransferData']['z_axis']
    nzt = len(z_axis)

    d1_lite = np.zeros((nzt, nm))
    d3_lite = np.zeros((nzt, nm))
    d4_lite = np.zeros((nzt, nm))

    for mid in np.arange(0, nm):
        if show_status:
            print('Get_dEdVdt status @ m:', mid/nm)
        
        m0 = m0_vec[mid]
        Lbh_vec = Lbh[:,mid]

        r1 = Get_dEdVdt_Kernel(m0 = m0, fbh0 = f0, Lbh_vec = Lbh_vec, z_vec = z_vec, channel = 1)['dEdVdt_dep']
        r3 = Get_dEdVdt_Kernel(m0 = m0, fbh0 = f0, Lbh_vec = Lbh_vec, z_vec = z_vec, channel = 3)['dEdVdt_dep']
        r4 = Get_dEdVdt_Kernel(m0 = m0, fbh0 = f0, Lbh_vec = Lbh_vec, z_vec = z_vec, channel = 4)['dEdVdt_dep']
        
        for zid in np.arange(0, nzt):
            d1_lite[zid, mid] = r1[zid]
            d3_lite[zid, mid] = r3[zid]
            d4_lite[zid, mid] = r4[zid]
    
    d1 = np.linspace(0, 1, nzt)
    d3 = np.linspace(0, 1, nzt)
    d4 = np.linspace(0, 1, nzt)
    
    # Integrate over initial mass function to get total dEdVdt
    
    lnm0 = np.log(m0_vec)
    for zid in np.arange(0, nzt):

        if show_status:
            print('Get_dEdVdt status @ z:', zid/nzt)
        
        F1 = d1_lite[zid,:] * MF0
        F3 = d3_lite[zid,:] * MF0
        F4 = d4_lite[zid,:] * MF0

        d1[zid] = np.trapz(x = lnm0, y = F1)
        d3[zid] = np.trapz(x = lnm0, y = F3)
        d4[zid] = np.trapz(x = lnm0, y = F4)
    
    # Now get fc, it's ok to have NaN because we don't use z>1000 in Evo
    dEdVdt_inj_vec = Evo['dEdVdt']
    dEdVdt_inj = np.interp(x = z_axis, xp = z_vec[::-1], fp = dEdVdt_inj_vec[::-1], left = 0, right = 0)
    
    f1 = d1/dEdVdt_inj
    f3 = d3/dEdVdt_inj
    f4 = d4/dEdVdt_inj
    # Let's also print OmB to check consistency with background HyRec cosmology
    z_vec = Evo['z_vec'][::-1]
    OmB_ratio_vec = Evo['omb_ratio'][::-1]
    OmB_ratio = np.interp(x = z_axis, xp = z_vec, fp = OmB_ratio_vec)

    if show_status:
        print('Time Used for Get_dEdVdt: ', time.time() - t1)

    r = {
        'z' : z_axis,
        'dEdVdt_HIon' : d1,
        'dEdVdt_LyA' : d3,
        'dEdVdt_Heat' : d4,
        'fc_HIon' : f1,
        'fc_LyA' : f3,
        'fc_Heat' : f4,
        'dEdVdt_inj' : dEdVdt_inj,
        'OmB_ratio' : OmB_ratio
        }
    
    return r


def Call_HyRec(
        fbh0 = 1e-8,
        mc = 1e2,
        sbh = 1,
        Use_Ps_domain = False,
        PsA = 0.1,
        kc = 1e5,
        SigmaPs = 0.5,
        zmin = 10,
        zmax = 500,
        nz = 2000,
        nm = 1000,
        DeltaC = 0.45,
        Use_EoR = False,
        dmdm0_method = 0,
        Use_Halo_Boost = False,
        sbh_range = [-3, 5],
        show_status = True,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Break_OmB_Lost_Ratio = 1.1,
        Fix_mdot = False,
        Fixed_mdot = 10,
        Show_Extrap_MSG = False,
        Ignore_Evo = False
    ):

    dEdVdt_Data = Get_dEdVdt(
        fbh0 = fbh0,
        mc = mc,
        sbh = sbh,
        Use_Ps_domain = Use_Ps_domain,
        PsA = PsA,
        kc = kc,
        SigmaPs = SigmaPs,
        zmin = zmin,
        zmax = zmax,
        nz = nz,
        nm = nm,
        DeltaC = DeltaC,
        Use_EoR = Use_EoR,
        dmdm0_method = dmdm0_method,
        Use_Halo_Boost = Use_Halo_Boost,
        sbh_range = sbh_range,
        show_status = show_status,
        Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
        Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
        Use_mdot_Interp = Use_mdot_Interp,
        Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
        Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
        Break_OmB_Lost_Ratio = Break_OmB_Lost_Ratio,
        Fix_mdot = Fix_mdot,
        Show_Extrap_MSG = Show_Extrap_MSG,
        Fixed_mdot = Fixed_mdot,
        Ignore_Evo = Ignore_Evo)
    
    # passwd is for mpi, 难道我真的是个天才
    passwd_length = 10
    passwd = 'tmp_'
    for idx in np.arange(0, passwd_length):
        passwd = passwd + str(np.random.randint(10))
    
    IF = passwd + '_in.txt'
    OF = passwd + '_out.txt'
    IF_Full = HyRec_Path() + IF
    OF_Full = HyRec_Path() + OF
    Print_dEdVdt(dEdVdt_dataset = dEdVdt_Data, File = IF_Full)
    cmd1 = 'cd ' + HyRec_Path() + ';'
    cmd2 = './hyrec <'
    cmd3 = IF + '>' + OF
    cmd = cmd1 + cmd2 + cmd3
    os.system(cmd)

    # Loading data
    HyRec_Data = np.loadtxt(OF_Full)
    
    # HyRec may fail for very high injection rate, I have tunned it to return NaN for this scenario
    NaN_check = np.isnan(HyRec_Data)
    if True in NaN_check:
        z = np.logspace(0,4,100)-1
        xe = np.nan * z
        t = np.nan * z
        warnings.warn('HyRec collapsed, injection rate too high, I am returning NaN')
    else:
        z = HyRec_Data[:,0]
        xe = HyRec_Data[:,1]
        t = HyRec_Data[:,2]
    
    # clean up
    if os.path.exists(OF_Full) and os.path.exists(IF_Full):
        os.remove(OF_Full)
        os.remove(IF_Full)
    
    OmB_z_axis = Interpolation_Table['TransferData']['z_axis']

    r = {
        'z':z, 
        'xe': xe, 
        't': t, 
        'OmB_ratio' : dEdVdt_Data['OmB_ratio'],
        'zb' : OmB_z_axis}

    return r

def Find_FbhMax_EDGES(
        T21 = -100,
        mc = 1e0,
        sbh = -1,
        zmax = 500,
        nz = 2000,
        nm = 1000,
        LgF_Min = -20,
        dmdm0_method = 0,
        Use_Halo_Boost = False,
        sbh_range = [-3, 5],
        show_status = True,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Fix_mdot = False,
        Fixed_mdot = 10,
        Show_Extrap_MSG = False,
        Ignore_Evo = False
    ):
    
    # Only search in this range
    LgF_Max = np.log10(0.9)
    
    # Note that we are using z=17
    TkMax = 1783.65/(36.324 - T21)
    
    def Tk_dif(LgF):
        RecHistory = Call_HyRec(
            fbh0 = 10**LgF,
            mc = mc,
            sbh = sbh,
            Use_Ps_domain = False,
            zmin = 16,
            zmax = zmax,
            nz = nz,
            nm = nm,
            dmdm0_method = dmdm0_method,
            Use_Halo_Boost = Use_Halo_Boost,
            sbh_range = sbh_range,
            show_status = False,
            Use_EoR = False,
            Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
            Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
            Use_mdot_Interp = Use_mdot_Interp,
            Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
            Break_OmB_Lost_Ratio = 1.1,
            Fix_mdot = Fix_mdot,
            Fixed_mdot = Fixed_mdot,
            Show_Extrap_MSG = Show_Extrap_MSG,
            Ignore_Evo = Ignore_Evo)
        
        z = RecHistory['z'][::-1]
        tk = RecHistory['t'][::-1]
        NaN_Check = np.isnan(tk)
        if True in NaN_Check:
            t17 = 1e5
        else:
            t17 = np.interp(x = 17, xp = z, fp = tk)
        dif = t17 - TkMax
        return dif

    try:
        LfMax = Solve(F = Tk_dif, Xmin = LgF_Min, Xmax = LgF_Max, Precision = 1e-2, show_status = show_status)
    except:
        # Normally this happens for LowM, but HiM might also be susceptible
        dif = Tk_dif(LgF_Min)
        if dif > 0:
            LfMax = LgF_Min
        else:
            LfMax = 0
    r = 10**LfMax

    return r

def Find_PsMax_EDGES(
        T21 = -100,
        kc = 1e5,
        SigmaPs = 0.5,
        LgF_Min = -30,
        zmax = 500,
        nz = 1000,
        nm = 500,
        DeltaC = 0.45,
        dmdm0_method = 0,
        Use_Halo_Boost = False,
        show_status = True,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Fix_mdot = False,
        Fixed_mdot = 10,
        Show_Extrap_MSG = False,
        Ignore_Evo = False
    ):
    '''
    Speed : 23m for nz = 1500, nm = 500
    '''
    
    # Note that we are using z=17
    TkMax = 1783.65/(36.324 - T21)
    
    LgF_Max = np.log10(0.7)
    # Smarter way to find scan region for A
    
    def Find_fbh_for_A(A):
        MF = Phi_Profile(
            A = A,
            Sigma = SigmaPs,
            kc = kc,
            DeltaC = DeltaC,
            map_nx = 500
        )
        if not np.isnan(MF['Width']):
            return MF['fbh']
        else:
            return 1e-280

    def Find_A(fbh):
        # each call takes about 10s
        dif = lambda la: Find_fbh_for_A(10**la) - fbh
        try:
            r = Solve(F = dif, Xmin = -4, Xmax = 2, Precision = 1e-3)
        except:
            # Expand search
            r = Solve(F = dif, Xmin = -8, Xmax = 6, Precision = 1e-3)
        r = 10**r
        return r
    
    A_Search_Min = Find_A(10**LgF_Min)
    A_Search_Max = Find_A(10**LgF_Max)
    LA_Search_Min = np.log10(A_Search_Min)
    LA_Search_Max = np.log10(A_Search_Max)
    
    def Tk_dif(LgA):
        RecHistory =  Call_HyRec(
            Use_Ps_domain = True,
            PsA = 10**LgA,
            kc = kc,
            SigmaPs = SigmaPs,
            zmin = 16,
            zmax = zmax,
            nz = nz,
            nm = nm,
            DeltaC = DeltaC,
            Use_EoR = False,
            dmdm0_method = dmdm0_method,
            Use_Halo_Boost = Use_Halo_Boost,
            show_status = False,
            Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
            Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
            Use_mdot_Interp = Use_mdot_Interp,
            Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
            Break_OmB_Lost_Ratio = 1.1,
            Fix_mdot = Fix_mdot,
            Fixed_mdot = Fixed_mdot,
            Show_Extrap_MSG = Show_Extrap_MSG,
            Ignore_Evo = Ignore_Evo)
        
        z = RecHistory['z'][::-1]
        tk = RecHistory['t'][::-1]
        NaN_Check = np.isnan(tk)
        if True in NaN_Check:
            t17 = 1e5
        else:
            t17 = np.interp(x = 17, xp = z, fp = tk)
        dif = t17 - TkMax
        return dif
    
    try:
        LgAmax = Solve(F = Tk_dif, Xmin = LA_Search_Min, Xmax = LA_Search_Max, Precision = 1e-2, show_status = show_status)
    except:
        # Normally this happens for LowA, but HiA might also be susceptible
        dif = Tk_dif(LA_Search_Min)
        if dif > 0:
            LgAmax = LA_Search_Min
        else:
            # Actually we can use this, this would be the constraint from fbh0 < 1
            LgAmax = LA_Search_Max
    
    Amax = 10**LgAmax
    Pmax = Amax/(np.sqrt(2 * np.pi)*SigmaPs)

    return Pmax

'''
t1 = TimeNow()
a = Find_PsMax_EDGES(Fix_mdot = 0)
print(a)        
Timer(t1)
'''
