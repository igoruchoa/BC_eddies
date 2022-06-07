import numpy as np

def vortoa(Xg,Yg,X,Y,U,V,corrlenx,corrleny,err,b=0):
    '''
    Vectoa is a vectorial objective analysis function.
    It interpolates a velocity field (U and V, east and north velocity components)
    into a streamfunction field (dimension MxN) solving the Laplace equation:
        $nabla^{2}Psi=0$.
    ======================================================================
    Input:
        Xg & Yg - Grid of of interpolation points (i.e. LON & LAT grid)
        X & Y - Arrays of observation points
        U & V - Arrays of observed east and north components of velocity
        corrlen & err - Correlation length scales (in x and y) and error for a
                    gaussian streamfunction covariance function (floats)
        b - Constant value that forces a correction in the data mean value.
            Unless it is defined, b=0 by default.
    ======================================================================
    Output:
        PSI - Streamfuction field matrix with MxN dimension.
              The dimension of the output is always the same of XC & YC
    ======================================================================
    PYTHON VERSION by:
       Iury Sousa and Hélio Almeida - 30 May 2016
       Laboratório de Dinâmica Oceânica - IOUSP
    ======================================================================'''
    # making sure that the input variables aren't changed
    xc,yc,x,y,u,v=Xg.copy(),Yg.copy(),X.copy(),Y.copy(),U.copy(),V.copy()

    corrlen = corrleny
    xc = xc*( corrleny*1./corrlenx)
    x = x*(corrleny*1./corrlenx)

    n = len(x)
    # Joins all the values of velocity (u then v) in one column-wise array.
    # Being u and v dimension len(u)=y, then uv.shape = (2*y,1)
    uv=np.array([np.hstack((u,v))]).T    #data entered row-wise


    # CALCULATING angles and distance

    # pp is a join of two matrix calculating the distance between every point of observation and all others like the
    # example below
    # len(y) = M -> pp[0].shape = MxM being:

    #    pp[0][i] = y-y[i]
    pp = -np.tile(y,(n,1)).T+np.tile(y,(n,1)),-np.tile(x,(n,1)).T+np.tile(x,(n,1))

    #
    t = []
    for ii,jj in zip(pp[0].ravel(),pp[1].ravel()):
        t.append(np.math.atan2(ii,jj))
    t = np.array(t)
    t.shape = pp[0].shape
    # t end up to be angles and d2 the distances between every observation point and all the others
    d2=((np.tile(x,(n,1)).T-np.tile(x,(n,1)))**2+(np.tile(y,(n,1)).T-np.tile(y,(n,1)))**2)



    lambd = 1/(corrlen**2)
    bmo=b*err/lambd
    R=np.exp(-lambd*d2)        #%longitudinal
    S=R*(1-2*lambd*d2)+bmo #%transverse
    R=R+bmo

    A=np.zeros((2*n,2*n))
    A[0:n,0:n]=(np.cos(t)**2)*(R-S)+S
    A[0:n,n:2*n]=np.cos(t)*np.sin(t)*(R-S)
    A[n:2*n,0:n]=A[0:n,n:2*n]
    A[n:2*n,n:2*n]=(np.sin(t)**2)*(R-S)+S
    A=A+err*np.eye(A.shape[0])

    # angles and distances
    nv1,nv2 =xc.shape
    nv=nv1*nv2


    xc = xc.T.ravel()
    yc = yc.T.ravel()


    #% the same idea of pp but this time for interpolation points
    ppc = -np.tile(yc,(n,1)).T+np.tile(y,(nv,1)),-np.tile(xc,(n,1)).T+np.tile(x,(nv,1))
    tc = []
    for ii,jj in zip(ppc[0].ravel(),ppc[1].ravel()):
        tc.append(np.math.atan2(ii,jj))
    tc = np.array(tc)
    tc.shape = ppc[0].shape
    d2=((np.tile(xc,(n,1)).T-np.tile(x,(nv,1)))**2+(np.tile(yc,(n,1)).T-np.tile(y,(nv,1)))**2)
    R=np.exp(-lambd*d2)+bmo;
    S=R*(1-2*lambd*d2)


    C=np.zeros((2*n,2*n))
    C[0:nv,0:n]=(np.cos(tc)**2)*(R-S)+S
    C[0:nv,n:2*n]=np.cos(tc)*np.sin(tc)*(R-S)
    C[nv:2*nv,0:n]=np.cos(tc)*np.sin(tc)*(R-S)
    C[nv:2*nv,n:2*n]=(np.sin(tc)**2)*(R-S)+S
    
    R=-4*lambd*(2-d2*lambd)*R
    P=np.zeros((nv,2*n))
    P[:,:n]=np.sin(tc)*np.sqrt(d2)*R
    P[:,n:2*n]=-np.cos(tc)*np.sqrt(d2)*R

    pv=np.dot(P,np.linalg.solve(A,uv)).reshape(nv2,nv1).T
    vc=None#np.dot(C,np.linalg.solve(A,uv)).reshape(nv2,nv1).T

    return pv,vc


def vectoa(Xg,Yg,X,Y,U,V,corrlenx,corrleny,err,b=0):
    '''
    Vectoa is a vectorial objective analysis function.
    It interpolates a velocity field (U and V, east and north velocity components)
    into a streamfunction field (dimension MxN) solving the Laplace equation:
        $nabla^{2}Psi=0$.
   ======================================================================
    Input:
        Xg & Yg - Grid of of interpolation points (i.e. LON & LAT grid)
        X & Y - Arrays of observation points
        U & V - Arrays of observed east and north components of velocity
        corrlen & err - Correlation length scales (in x and y) and error for a
                    gaussian streamfunction covariance function (floats)
        b - Constant value that forces a correction in the data mean value.
            Unless it is defined, b=0 by default.
   ======================================================================
    Output:
        PSI - Streamfuction field matrix with MxN dimension.
              The dimension of the output is always the same of XC & YC
   ======================================================================
    PYTHON VERSION by:
       Iury Sousa and Hélio Almeida - 30 May 2016
       Laboratório de Dinâmica Oceânica - IOUSP
   ======================================================================'''
   # making sure that the input variables aren't changed
    xc,yc,x,y,u,v=Xg.copy(),Yg.copy(),X.copy(),Y.copy(),U.copy(),V.copy()

    corrlen = corrleny
    xc = xc*( corrleny*1./corrlenx)
    x = x*(corrleny*1./corrlenx)

    n = len(x)
    # Joins all the values of velocity (u then v) in one column-wise array.
    # Being u and v dimension len(u)=y, then uv.shape = (2*y,1)
    uv=np.array([np.hstack((u,v))]).T    #data entered row-wise


    # CALCULATING angles and distance

    # pp is a join of two matrix calculating the distance between every point of observation and all others like the
    # example below
    # len(y) = M -> pp[0].shape = MxM being:

    #    pp[0][i] = y-y[i]
    pp = -np.tile(y,(n,1)).T+np.tile(y,(n,1)),-np.tile(x,(n,1)).T+np.tile(x,(n,1))

    #
    t = []
    for ii,jj in zip(pp[0].ravel(),pp[1].ravel()):
        t.append(np.math.atan2(ii,jj))
    t = np.array(t)
    t.shape = pp[0].shape
    # t end up to be angles and d2 the distances between every observation point and all the others
    d2=((np.tile(x,(n,1)).T-np.tile(x,(n,1)))**2+(np.tile(y,(n,1)).T-np.tile(y,(n,1)))**2)



    lambd = 1/(corrlen**2)
    bmo=b*err/lambd
    R=np.exp(-lambd*d2)        #%longitudinal
    S=R*(1-2*lambd*d2)+bmo #%transverse
    R=R+bmo

    A=np.zeros((2*n,2*n))

    A[0:n,0:n]=(np.cos(t)**2)*(R-S)+S
    A[0:n,n:2*n]=np.cos(t)*np.sin(t)*(R-S)
    A[n:2*n,0:n]=A[0:n,n:2*n]
    A[n:2*n,n:2*n]=(np.sin(t)**2)*(R-S)+S
    A=A+err*np.eye(A.shape[0])

    # angles and distances
    nv1,nv2 =xc.shape
    nv=nv1*nv2


    xc = xc.T.ravel()
    yc = yc.T.ravel()


    #% the same idea of pp but this time for interpolation points
    ppc = -np.tile(yc,(n,1)).T+np.tile(y,(nv,1)),-np.tile(xc,(n,1)).T+np.tile(x,(nv,1))
    tc = []
    for ii,jj in zip(ppc[0].ravel(),ppc[1].ravel()):
        tc.append(np.math.atan2(ii,jj))
    tc = np.array(tc)
    tc.shape = ppc[0].shape
    d2=((np.tile(xc,(n,1)).T-np.tile(x,(nv,1)))**2+(np.tile(yc,(n,1)).T-np.tile(y,(nv,1)))**2)
    R=np.exp(-lambd*d2)+bmo;


    P=np.zeros((nv,2*n))
    # streamfunction-velocity covariance
    P[:,0:n]=np.sin(tc)*np.sqrt(d2)*R;
    P[:,n:2*n]=-np.cos(tc)*np.sqrt(d2)*R;

    PSI=np.dot(P,np.linalg.solve(A,uv))   # solvi the linear system
    PSI=PSI.reshape(nv2,nv1).T

    return PSI


def scaloa(xc, yc, x, y, t=None, corrlenx=None,corrleny=None, err=None, zc=None):
    """
    Scalar objective analysis. Interpolates t(x, y) into tp(xc, yc)
    Assumes spatial correlation function to be isotropic and Gaussian in the
    form of: C = (1 - err) * np.exp(-d**2 / corrlen**2) where:
    d : Radial distance from the observations.
    Parameters
    ----------
    corrlen : float
    Correlation length.
    err : float
    Random error variance (epsilon in the papers).
    Return
    ------
    tp : array
    Gridded observations.
    ep : array
    Normalized mean error.
    Examples
    --------
    See https://ocefpaf.github.io/python4oceanographers/blog/2014/10/27/OI/
    Notes
    -----
    The funcion `scaloa` assumes that the user knows `err` and `corrlen` or
    that these parameters where chosen arbitrary. The usual guess are the
    first baroclinic Rossby radius for `corrlen` and 0.1 e 0.2 to the sampling
    error.
    """
    corrlen = corrleny
    xc = xc*( corrleny*1./corrlenx)
    x = x*(corrleny*1./corrlenx)

    n = len(x)
    x, y = np.reshape(x, (1, n)), np.reshape(y, (1, n))
    # Squared distance matrix between the observations.
    d2 = ((np.tile(x, (n, 1)).T - np.tile(x, (n, 1))) ** 2 +
    (np.tile(y, (n, 1)).T - np.tile(y, (n, 1))) ** 2)
    nv = len(xc)
    xc, yc = np.reshape(xc, (1, nv)), np.reshape(yc, (1, nv))
    # Squared distance between the observations and the grid points.
    dc2 = ((np.tile(xc, (n, 1)).T - np.tile(x, (nv, 1))) ** 2 +
    (np.tile(yc, (n, 1)).T - np.tile(y, (nv, 1))) ** 2)
    # Correlation matrix between stations (A) and cross correlation (stations
    # and grid points (C))
    A = (1 - err) * np.exp(-d2 / corrlen ** 2)
    C = (1 - err) * np.exp(-dc2 / corrlen ** 2)
    if 0: # NOTE: If the parameter zc is used (`scaloa2.m`)
        A = (1 - d2 / zc ** 2) * np.exp(-d2 / corrlen ** 2)
        C = (1 - dc2 / zc ** 2) * np.exp(-dc2 / corrlen ** 2)
    # Add the diagonal matrix associated with the sampling error. We use the
    # diagonal because the error is assumed to be random. This means it just
    # correlates with itself at the same place.
    A = A + err * np.eye(len(A))
    # Gauss-Markov to get the weights that minimize the variance (OI).
    tp = None
    ep = 1 - np.sum(C.T * np.linalg.solve(A, C.T), axis=0) / (1 - err)
    if t is not None:
        t = np.reshape(t, (n, 1))
        tp = np.dot(C, np.linalg.solve(A, t))
        #if 0: # NOTE: `scaloa2.m`
        #  mD = (np.sum(np.linalg.solve(A, t)) /
        #  np.sum(np.sum(np.linalg.inv(A))))
        #  t = t - mD
        #  tp = (C * (np.linalg.solve(A, t)))
        #  tp = tp + mD * np.ones(tp.shape)
        return tp, ep

    else:
        print("Computing just the interpolation errors.")
        #Normalized mean error. Taking the squared root you can get the
        #interpolation error in percentage.
        return ep


def uv2zeta(x,y,U,V):
    '''
     ZETA k-component of rotational by velocity field
       Returns the scalar field of vorticity from U and V
       velocity components. X and Y are longitude and latitude
       matrices, respectively, as U and V. All matrices have
       same dimension. Grid could be curvilinear, but the error
       will increase.
       --> !!! IMPORTANT !!!  <-----------------------------------------
       The grid indexes IM and JM must have origin on the left-inferior
       corner, increasing to right and to up, respectively.
                                  GRID
                :            :            :            :
                |            |            |            |
       (JM = 2) o ---------- o ---------- o ---------- o --- ...
                |            |            |            |
                |            |            |            |
       (JM = 1) o ---------- o ---------- o ---------- o --- ...
                |            |            |            |
                |            |            |            |
       (JM = 0) o ---------- o ---------- o ---------- o --- ...
            (IM = 0)     (IM = 1)     (IM = 2)     (IM = 3)
       -----------------------------------------------------------------
      USAGE:  u,v = psi2uv(x,y,psi)
      INPUT:
         x     = decimal degrees (+ve E, -ve W) [-180..+180]
         y     = decimal degrees (+ve N, -ve S) [- 90.. +90]
         u     = velocity zonal component [m s^-1]
         v     = velocity meridional component [m s^-1]
      OUTPUT:
         ZETA  = Relative vorticity field [s^-1]
      EXAMPLE:
        #Subset function for quiver plot
        st = 2
        sub = lambda P: P[::st,::st]
        #Number of points
        q = 30
        #Grid
        X,Y = np.meshgrid(np.linspace(-10,10,q),np.linspace(-20,20,q))
        #Length in km
        ly,lx = 40*60*1852,10*60*1852
        #Velocity field
        U     = np.cos(60*1852*2*np.pi*Y/ly)
        V     = np.cos((np.pi/2)+60*1852*2*np.pi*X/lx)
        #Correction to longitudes
        coslat = np.cos(np.pi*Y/180)
        #Analytical Vorticity
        Va = -(2*np.pi/(lx*coslat))*np.sin((np.pi/2)+60*1852*2*np.pi*X/lx)+\
            (2*np.pi/ly)*np.sin(60*1852*2*np.pi*Y/ly)
        #Numerical Vorticity
        Vn = rel_vort(X,Y,U,V)
        #Scale factor to apply on Vn
        fac = Va.max()/Vn.max()
        #Figure
        #args for zeta plot
        kw = {'cmap':'RdBu_r','alpha':0.8}
        fig,(a1,a2,a3) = plt.subplots(1,3,figsize=(15,8))
        a1.pcolormesh(X,Y,Va,**kw)
        a1.quiver(sub(X),sub(Y),sub(U),sub(V),scale=10,linewidths=0.5)
        a1.set_title('ANALYTICAL')
        C = a2.contour(X,Y,(Va-Vn)*100/Va.max(),np.arange(0,100,5),colors='k')
        a2.clabel(C,fmt='%2i')
        a2.set_title('ERROR [%]')
        A = a3.pcolormesh(X,Y,Vn*fac,**kw)
        fig.colorbar(A)
        a3.quiver(sub(X),sub(Y),sub(U),sub(V),scale=10,linewidths=0.5)
        a3.set_title('NUMERIC')
      AUTHOR:
       Wandrey Watanabe e Iury Sousa    - 29 Jun 2016
       Laboratório de Dinâmica Oceânica - IOUSP
       ======================================================================
    '''
    #função para cálculo de ângulos
    angcalc = lambda dy,dx: np.math.atan2(dy,dx)

    #funções para os diferenciais com borda
    dX      = lambda d: np.hstack([(d[:,1]-d[:,0])[:,None],
                        d[:,2:]-d[:,:-2],(d[:,-1]-d[:,-2])[:,None]])
    dY      = lambda d: np.vstack([(d[1,:]-d[0,:])[None,:],
                        d[2:,:]-d[:-2,:],(d[-1,:]-d[-2,:])[None,:]])

    # x = coordenada natural i
    # y = coordenada natural j
    # X = coordenada zonal cartesiana
    # Y = coordenada meridional cartesiana

    dyX = dX(y)*60*1852*np.cos(y*np.pi/180)
    dyY = dY(y)*60*1852

    dxX = dX(x)*60*1852*np.cos(y*np.pi/180)
    dxY = dY(x)*60*1852

    dsx = np.sqrt(dxX**2 + dyX**2)
    dsy = np.sqrt(dxY**2 + dyY**2)

    #Derivadas naturais de U
    dUy = dY(U)
    dUx = dX(U)

    #Derivadas naturais de V
    dVy = dY(V)
    dVx = dX(V)

    #a função map aplicará a função lambda angcalc em todos os elementos
    #indicados depois da virgula e armazenará em uma lista
    angX = list(map(angcalc,dyX.ravel(),dxX.ravel()))
    angY = list(map(angcalc,dyY.ravel(),dxY.ravel()))

    #fazendo rechape dos ângulos para o formato das matrizes de distâncias
    angX = np.reshape(angX,dsx.shape)
    angY = np.reshape(angY,dsy.shape)

    #calculando as componentes X e Y das derivadas
    dv1,dv2 =  (dVy/dsy)*np.cos(angY),  (dVx/dsx)*np.cos(angX)
    du1,du2 = -(dUy/dsy)*np.sin(angY), -(dUx/dsx)*np.sin(angX)


    # zerando valores aos quais a derivada tendeu ao infinito:
    # isto só acontece se uma das dimensões da grade for paralalela a um
    # eixo cartesiano
    dv1[np.isinf(dv1)] = 0
    dv2[np.isinf(dv2)] = 0
    du1[np.isinf(du1)] = 0
    du2[np.isinf(du2)] = 0

    #somando as componentes
    ZETA = dv1+dv2+du1+du2

    return ZETA
