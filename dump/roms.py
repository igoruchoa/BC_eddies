
import numpy as np


def zlev(ncw):

    h = ncw['h'][:]
    hc = ncw['hc'][:]
    theta_s,theta_b,N = 6,0,30
    type='r'
    vtransform=2
    zeta = ncw['zeta'][:]

    comp = np.arange(0,30)

    for d in np.arange(0,zeta.shape[0]):
        z = zlevs(h,zeta[d,:,:],theta_s,theta_b,hc,N,type,vtransform);
        i,j,k = z.shape
        zn = z.reshape(i,j*k)

        for zi in zn.T:
            I = np.argsort(zi)
            ok = (I == comp);
            if any(~ok):
                print('DEU PAU LIXO')

#             print('dia %s funcionou!'%d)
#         print('mes funcionou!')
    return z



def zlevs(h,zeta,theta_s,theta_b,hc,N,type='r',vtransform=2):
    """
    ################################################################
    #
    #  function z = zlevs(h,zeta,theta_s,theta_b,hc,N,type);
    #
    #  this function compute the depth of rho or w points for ROMS
    #
    #  On Input:
    #
    #    type    'r': rho point 'w': w point
    #
    #  On Output:
    #
    #    z       Depths (m) of RHO- or W-points (3D matrix).
    #
    #  Further Information:
    #  http://www.brest.ird.fr/Roms_tools/
    #
    #  This file is part of ROMSTOOLS
    #
    #  ROMSTOOLS is free software; you can redistribute it and/or modify
    #  it under the terms of the GNU General Public License as published
    #  by the Free Software Foundation; either version 2 of the License,
    #  or (at your option) any later version.
    #
    #  ROMSTOOLS is distributed in the hope that it will be useful, but
    #  WITHOUT ANY WARRANTY; without even the implied warranty of
    #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #  GNU General Public License for more details.
    #
    #  You should have received a copy of the GNU General Public License
    #  along with this program; if not, write to the Free Software
    #  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
    #  MA  02111-1307  USA
    #
    #  Copyright (c) 2002-2006 by Pierrick Penven
    #  e-mail:Pierrick.Penven@ird.fr
    #  Translated to Python by Rafael Soutelino, rsoutelino@gmail.com
    #  Last Modification: Aug, 2010
    ################################################################
    """

    M, L = h.shape

    # Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.

    sc_r=np.zeros(N)
    Cs_r=np.zeros(N)
    sc_w=np.zeros(N+1)
    Cs_w=np.zeros(N+1)

    if (vtransform == 2):
        ds=1./N;
        if type=='w':

            sc_w[0] = -1.0
            sc_w[-1] =  0
            Cs_w[0] = -1.0
            Cs_w[-1] =  0

            sc_w[1:N] = ds*(np.arange(1,N)-N)
            Cs_w = csf(sc_w, theta_s,theta_b)
            N+=1
    #    Still in matlab code
    #    disp(['===================================='])
    #    for k=N:-1:1
    #        disp(['Niveau S=',num2str(k),' Cs=',num2str( Cs_w(k), '%8.7f')])
    #    end
    #    disp(['===================================='])

        else:

            sc = ds*(np.arange(1,N+1)-N-0.5)
            Cs_r = csf(sc, theta_s,theta_b)
            sc_r = sc.copy()
    #    Still in matlab code
    #    disp(['===================================='])
    #    for k=N:-1:1
    #        disp(['Niveau S=',num2str(k),' Cs=',num2str( Cs_r(k), '%8.7f')])
    #    end
    #    disp(['===================================='])

    else:
        cff1 = 1.0 / np.sinh( theta_s )
        cff2 = 0.5 / np.tanh( 0.5*theta_s )
        ##
        if type=='w':
            sc = ( np.arange(0,N+1) - N ) / N
            N  = N + 1
        else:
            sc = ( np.arange(1,N+1) - N - 0.5 ) / N
        ##

        Cs = (1 - theta_b) * cff1 * np.sinh( theta_s * sc ) + \
            theta_b * ( cff2 * np.tanh(theta_s *(sc + 0.5) ) - 0.5 )

    #   Still in matlab code
    #    disp(['===================================='])
    #    for k=N:-1:1
    #        disp(['Niveau S=',num2str(k),' Cs=',num2str( Cs(k), '%8.7f')])
    #    end
    #    disp(['===================================='])

    # ------------------------------------------------------------------ #

    # Create S-coordinate system: based on model topography h(i,j),
    # fast-time-averaged free-surface field and vertical coordinate
    # transformation metrics compute evolving depths of of the three-
    # dimensional model grid. Also adjust zeta for dry cells.

    Dcrit = 0.2 # min water depth in dry cells
    h[h==0] = 1e-14
    zeta[zeta<(Dcrit-h)] = Dcrit - h[zeta<(Dcrit-h)]
    hinv = 1./h
    z = np.zeros([N,M,L])

    if (vtransform == 2):

        if type=='w':
            cff1 = Cs_w.copy()
            cff2 = sc_w.copy() + 1
            sc = sc_w.copy()
        else:
            cff1 = Cs_r.copy()
            cff2 = sc_r.copy() + 1
            sc = sc_r.copy()

        h2 = h + hc
        cff = hc*sc
        h2inv = 1./h2

        for k in np.arange(0,N):
            z0 = cff[k] + cff1[k]*h
            z[k,:,:] = (z0*h)/h2 + zeta*(1.+ z0*h2inv)

    else:

        cff1 = Cs.copy()
        cff2 = sc + 1
        cff = hc*(sc-Cs)

        for k in range(0, N):
            z0         = cff[k] + cff1[k]*h
            z[k,:,:] = z0 + zeta*( 1 + z0*hinv )

    return z


def csf(sc, theta_s,theta_b):
    #######################################################################
    #
    #  function h = csf(sc, theta_s,theta_b);
    #
    #  Further Information:
    #  http://www.brest.ird.fr/Roms_tools/
    #
    #  This file is part of ROMSTOOLS
    #
    #  ROMSTOOLS is free software; you can redistribute it and/or modify
    #  it under the terms of the GNU General Public License as published
    #  by the Free Software Foundation; either version 2 of the License,
    #  or (at your option) any later version.
    #
    #  ROMSTOOLS is distributed in the hope that it will be useful, but
    #  WITHOUT ANY WARRANTY; without even the implied warranty of
    #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #  GNU General Public License for more details.
    #
    #  You should have received a copy of the GNU General Public License
    #  along with this program; if not, write to the Free Software
    #  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
    #  MA  02111-1307  USA
    #
    ######################################################################

    if (theta_s > 0 ):
        csrf=(1-np.cosh(sc*theta_s))/(np.cosh(theta_s)-1)
    else:
        csrf=-sc**2

    if (theta_b > 0):
        h = (np.exp(theta_b*csrf)-1)/(1-np.exp(-theta_b))
    else:
        h  = csrf.copy();

    return h