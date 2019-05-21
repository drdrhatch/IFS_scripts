from ParIO import *
import numpy as np

def init_read_parameters(suffix):

    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict

    return pars

def otherRef(suffix, pars):

    qref = 1.6E-19
    c  = 1.
    m_kg = 1.673E-27
    Bref = pars['Bref']
    Tref = pars['Tref']
    nref = pars['nref']
    Lref = pars['Lref']
    mref = pars['mref']
    nref = nref * 1.E19
    Tref = Tref * qref * 1.E03
    mref = mref * m_kg
    pref = nref * Tref
    cref = np.sqrt(Tref / mref)
    Omegaref = qref * Bref / mref / c
    rhoref = cref / Omegaref
    rhorefStar = rhoref / Lref
    
    Apar_norm = Bref * Lref * rhorefStar ** 2
    phi_norm = Tref / qref * rhorefStar

    Gamma_gb = cref * nref * rhorefStar ** 2
    Qheat_gb = cref * pref * rhorefStar ** 2

    print 'Tref = ', Tref, 'J'
    print 'nref = ', nref, 'm^-3'
    print 'pref = ', pref, 'N / m^2'
    print 'cref = ', cref, 'm / s'
    print 'Omegaref = ', Omegaref, 'Hz'
    print 'rhoref = ', rhoref, 'm'
    print 'rhorefStar = ', rhorefStar, 'dimensionless'
    print 'phiref = ',phi_norm, 'V'
    print 'Aparref = ', Apar_norm, 'Tesla * m'
    print 'Gamma_gb = ', Gamma_gb, '/ m^2 / s'
    print 'Q_gb =', Qheat_gb, 'J / m^2 / s'

    return pref, cref, Omegaref, rhoref, rhorefStar, Apar_norm, Gamma_gb, Qheat_gb
