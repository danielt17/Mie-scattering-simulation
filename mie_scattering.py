# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 19:36:22 2022

@author: danie
"""

# %% Imports

import numpy as np
import matplotlib.pyplot as plt

# %% Functions

def mie_scattering(x,refrel,nang):
    '''
    This code is based upon Craig F. Bohren, Donald R. Huffman - Absorption and 
    Scattering of Light by Small Particles (1983, Wiley). The code calculates
    the scattering efficiencies of a homogneous sphere. From fortarn to python,
    taken from Appendix A of the book (Page 477 - 482).
    Inputs:
        x       - Size parameter = k*radius = 2pi/lambda * radius 
                  (lambda is the wavelength in the medium around the scatterers)
        refrel  - Refraction index (n in complex form for example:  1.5+0.02*1j)
        nang    - Number of angles for S1 and S2 function in range from 0 to pi/2
    outputs:
        S1, S2  - Funtion which correspond to the (complex) phase functions
        Qext    - Extinction efficiency (extintion is absorbed + scattered)
        Qsca    - Scattering efficiency 
        Qback   - Backscatter efficiency
        gsca    - Asymmetry parameter
    '''
    s1_1=np.zeros(nang,dtype=np.complex128); s1_2=np.zeros(nang,dtype=np.complex128)
    s2_1=np.zeros(nang,dtype=np.complex128); s2_2=np.zeros(nang,dtype=np.complex128)
    pi=np.zeros(nang); tau=np.zeros(nang)
    mxnang = 1000
    if (nang > mxnang):
        raise Exception('Error: maximum number of angles is: ' + str(mxnang) + ', but got: ' + str(nang) + ' angles.')
    if (nang < 2):
        raise Exception('Require number of angles > 1 in order to calculate scattering intensities')
    pii = 4.*np.arctan(1)
    dx = x
    drefrl = refrel
    y = x*drefrl
    ymod = abs(y)
    # Series expansion terminated after NSTOP terms
    # Logarithmic derivatives calculated from NMX on down
    xstop = x + 4*np.cbrt(x) + 2
    nmx = max(xstop,ymod) + 15.0
    nmx= np.fix(nmx)
    nstop = int(xstop)
    # Maximum number of recursion calls
    nmxx = 1e6 
    if (nmx > nmxx):
        raise Exception("Error maximum number of recursion calls is: " + str(int(nmxx)) + ", but current parameters require: "  + str(int(ymod)) + ' calls.')
    dang = 0.5*pii/ (nang-1)
    amu=np.arange(0.0,nang,1)
    amu=np.cos(amu*dang)
    pi0=np.zeros(nang); pi1=np.ones(nang)
    # Logarithmic derivative D(J) calculated by downward recurrence. 
    # beginning with initial value (0.,0.) at J=NMX
    nn = int(nmx)-1
    d=np.zeros(nn+1,dtype=np.complex128)
    for n in range(0,nn):
        en = nmx - n
        d[nn-n-1] = (en/y) - (1/ (d[nn-n]+en/y))
    # Riccati-Bessel functions with real argument X calculated by upward recurrence
    psi0 = np.cos(dx); psi1 = np.sin(dx)
    chi0 = -np.sin(dx); chi1 = np.cos(dx)
    xi1 = psi1-chi1*1j
    qsca = 0
    gsca = 0
    p = -1
    an = 0
    bn = 0
    for n in range(0,nstop):
        en = n+1.0
        fn = (2.*en+1.)/(en* (en+1.))
        # for given N, PSI  = psi_n        CHI  = chi_n
        # PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
        # PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
        # Calculate psi_n and chi_n
        psi = (2.*en-1.)*psi1/dx - psi0
        chi = (2.*en-1.)*chi1/dx - chi0
        xi = psi-chi*1j
        # Store previous values of AN and BN for use
        # in computation of g=<cos(theta)>
        if (n > 0):
            an1 = an
            bn1 = bn
        # Compute AN and BN:
        an = (d[n]/drefrl+en/dx)*psi - psi1
        an = an/ ((d[n]/drefrl+en/dx)*xi-xi1)
        bn = (drefrl*d[n]+en/dx)*psi - psi1
        bn = bn/ ((drefrl*d[n]+en/dx)*xi-xi1)
        # Augment sums for Qsca and g=<cos(theta)>
        qsca += (2.*en+1.)* (abs(an)**2+abs(bn)**2)
        gsca += ((2.*en+1.)/ (en* (en+1.)))*( np.real(an)* np.real(bn)+np.imag(an)*np.imag(bn))
        if (n > 0):
            gsca += ((en-1.)* (en+1.)/en)*( np.real(an1)* np.real(an)+np.imag(an1)*np.imag(an)+np.real(bn1)* np.real(bn)+np.imag(bn1)*np.imag(bn))
        # Now calculate scattering intensity pattern
        # First do angles from 0 to 90
        pi=0+pi1    # 0+pi1 because we want a hard copy of the values
        tau=en*amu*pi-(en+1.)*pi0
        s1_1 += fn* (an*pi+bn*tau)
        s2_1 += fn* (an*tau+bn*pi)
        # Now do angles greater than 90 using PI and TAU from
        # angles less than 90.
        # P=1 for N=1,3,...% P=-1 for N=2,4,...
        # remember that we have to reverse the order of the elements
        # of the second part of s1 and s2 after the calculation
        p = -p
        s1_2+= fn*p* (an*pi-bn*tau)
        s2_2+= fn*p* (bn*pi-an*tau)
        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1 = psi1-chi1*1j
        # Compute pi_n for next value of n
        # For each angle J, compute pi_n+1
        # from PI = pi_n , PI0 = pi_n-1
        pi1 = ((2.*en+1.)*amu*pi- (en+1.)*pi0)/ en
        pi0 = 0+pi   # 0+pi because we want a hard copy of the values
    # Have summed sufficient terms.
    # Now compute QSCA,QEXT,QBACK,and GSCA
    # we have to reverse the order of the elements of the second part of s1 and s2 
    s1=np.concatenate((s1_1,s1_2[-2::-1])); s2=np.concatenate((s2_1,s2_2[-2::-1]))
    gsca = 2.*gsca/qsca
    qsca = (2./ (dx*dx))*qsca
    qext = (4./ (dx*dx))* np.real(s1[0])
    qback = 4*(abs(s1[2*nang-2])/dx)**2    
    return s1,s2,qext,qsca,qback,gsca

# %% Main

if __name__ == '__main__':
    radius = 0.000001                                                       # radius of spherical particle
    center_wave_length = radius                                             # center wave length
    lower_wave_length = center_wave_length * 1e-2                           # start of wave length vector
    upper_wave_length = center_wave_length * 1e2                            # end of wave length vector
    wave_lengths = np.linspace(lower_wave_length,upper_wave_length,100000)  # wave length of EM wave
    xs = (2*np.pi*radius)/wave_lengths                                      # size parameter
    G = np.pi * radius ** 2                                                 # Particle cross section area for a sphere
    refrel = 3+0.02j                                                      # reflective coefficent
    nang = 180                                                              # number of angles to calculate at
    qext_ar,qsca_ar,qback_ar,gsca_ar = [],[],[],[]
    for x in xs:
        _,_,qext,qsca,qback,gsca = mie_scattering(x,refrel,nang)
        qext_ar.append(qext)
        qsca_ar.append(qsca)
        qback_ar.append(qback)
        gsca_ar.append(gsca)
    qabs_ar = list(np.array(qext_ar) - np.array(qsca_ar))
    plt.figure()
    plt.subplot()
    plt.semilogx(wave_lengths*100,qext_ar,'-',linewidth=3,label='$Q_{extinction}$')
    plt.semilogx(wave_lengths*100,qsca_ar,'-',linewidth=3,label='$Q_{scattering}$')
    plt.semilogx(wave_lengths*100,qabs_ar,'-',linewidth=3,label='$Q_{absorption}$')
    plt.semilogx(wave_lengths*100,qback_ar,'-',linewidth=3,label='$Q_{backscatter}$')
    plt.title('Mie scattering from a homogeneous sphere ($radius: r = ' + str(radius*100) + ' [cm]$, refractive index: $n_{sphere} = ' + str(refrel) + ' $), scattering efficiencies.')
    plt.xlabel('Wavelength [cm]')
    plt.ylabel('Efficiency')
    plt.xlim(lower_wave_length*100*1e1, upper_wave_length*1e-1*100)
    plt.legend(loc='upper right')
    plt.figure()
    plt.semilogx(wave_lengths*100,gsca_ar,linewidth=2)
    plt.title('Mie scattering from a homogeneous sphere ($radius: r = ' + str(radius*100) + ' [cm]$, refractive index: $n_{sphere} = ' + str(refrel) + ' $), Asymmetry coefficent.')
    plt.xlabel('Wavelength [cm]')
    plt.ylabel('Efficiency')
    plt.xlim(lower_wave_length*100*1e1, upper_wave_length*1e-1*100)