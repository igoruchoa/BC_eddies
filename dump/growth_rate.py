


import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



cfe = pd.read_pickle('/home/igor/Dropbox/META_eddydetection_images/data/CFE.pkl').values
cfe = xr.DataArray(cfe,dims=['eddy','day'],coords={'eddy':np.arange(cfe.shape[0]),'day':np.arange(cfe.shape[1])})

cste = pd.read_pickle('/home/igor/Dropbox/META_eddydetection_images/data/CSTE.pkl').values
cste = xr.DataArray(cste,dims=['eddy','day'],coords={'eddy':np.arange(cste.shape[0]),'day':np.arange(cste.shape[1])})


ie = pd.read_pickle('/home/igor/Dropbox/META_eddydetection_images/data/IE.pkl').values
ie = xr.DataArray(ie,dims=['eddy','day'],coords={'eddy':np.arange(ie.shape[0]),'day':np.arange(ie.shape[1])})

rce = pd.read_pickle('/home/igor/Dropbox/META_eddydetection_images/data/RCE.pkl').values
rce = xr.DataArray(rce,dims=['eddy','day'],coords={'eddy':np.arange(rce.shape[0]),'day':np.arange(rce.shape[1])})

sve = pd.read_pickle('/home/igor/Dropbox/META_eddydetection_images/data/SVE.pkl').values
sve = xr.DataArray(sve,dims=['eddy','day'],coords={'eddy':np.arange(sve.shape[0]),'day':np.arange(sve.shape[1])})

ve = pd.read_pickle('/home/igor/Dropbox/META_eddydetection_images/data/VE.pkl').values
ve = xr.DataArray(ve,dims=['eddy','day'],coords={'eddy':np.arange(ve.shape[0]),'day':np.arange(ve.shape[1])})

ae = pd.read_pickle('/home/igor/Dropbox/META_eddydetection_images/data/AE.pkl').values
ae = xr.DataArray(ae,dims=['eddy','day'],coords={'eddy':np.arange(ae.shape[0]),'day':np.arange(ae.shape[1])})


avg = lambda D,std: D.mean('eddy')+std*D.std('eddy') #media das curvas + X vezes o std das curvas
gwr = lambda D,std: (avg(D,std).differentiate('day')/avg(D,std)) # media da diff das curvas sobre a curva
std = lambda D: np.abs(np.vstack([gwr(D,-1)-gwr(D,0),gwr(D,1)-gwr(D,0)])).max(0)

func = lambda x,a,b,c: a*np.exp(-x/b) # função exponencial





popt,_=curve_fit(func, cste.day, gwr(cste,0))





# def get_fit(D):
#     fit = []
#     for gwri in (D.differentiate('day')/D).values:
#         popt,_=curve_fit(func, np.arange(len(gwri))[~np.isnan(gwri)], gwri[~np.isnan(gwri)])
#         fit.append(func(np.arange(len(gwri)),*popt))
#     fit = np.vstack(fit)*xr.ones_like(cfe)
#     fit = fit.where((fit.min('day')>=0)&(fit.max('day')<=0.4),drop=True)

#     return fit


def get_fit(D):
    fit = []
    for gwri in (D.differentiate('day')/D).values:
        popt,_=curve_fit(func, np.arange(len(gwri))[~np.isnan(gwri)], gwri[~np.isnan(gwri)])
        fit.append(func(np.arange(len(gwri)),*popt))
    fit = np.vstack(fit)*xr.ones_like(D)
    fit = fit.where((fit.min('day')>=0)&(fit.max('day')<=0.4),drop=True)

    return fit



fit=get_fit(cfe)


fig,ax = plt.subplots()

popt,_=curve_fit(func, cfe.day, (cfe.differentiate('day')/cfe).mean('eddy'))
ax.plot(cfe.day,func(cfe.day,*popt),lw=2.5,color='0.2',zorder=1)
ax.fill_between(fit.day,func(cfe.day,*popt)-fit.std('eddy'),func(cfe.day,*popt)+fit.std('eddy'),color='0.8',zorder=0)
(cfe.differentiate('day')/cfe).mean('eddy').plot(color='tab:red',zorder=5)
ax.set(xlim=[0,30],ylim=[-0.01,0.25])



plt.figure()

(cfe.differentiate('day')/cfe).plot(x='day',robust=True)




plt.figure()
(cste.differentiate('day')/cfe).plot(x='day',robust=True)





kw = {
    'xlim':[0,20],
    'ylim':[-0.015,0.2]
}

ne=2
'''
# fig,ax = plt.subplots(1,2,figsize=(6,3))
fig,ax = plt.subplots(2,4,figsize=(10,8))


'Ilheus Eddy'
# popt,_=curve_fit(func, ie.day, (ie.differentiate('day')/ie).mean('eddy'))
popt,_=curve_fit(func, ie.day[~np.isnan((ie.differentiate('day')/ie).mean('eddy'))]
                        ,(ie.differentiate('day')/ie).mean('eddy')[~np.isnan((ie.differentiate('day')/ie).mean('eddy'))])

ax[0,0].plot(ie.day,func(ie.day,*popt),lw=2.5,color='0.2',zorder=1)
fitstd = fit.std('eddy')
ax[0,0].fill_between(fit.day,func(ie.day,*popt)-fitstd,func(ie.day,*popt)+fitstd,color='0.9',zorder=0)
(ie.differentiate('day')/ie).mean('eddy').plot(ax=ax[0,0],marker='o',color='lightsalmon',lw=0,markeredgecolor='0.3',zorder=5)
ax[0,0].axvline(ne*popt[1],color='0.2',linestyle='--')
ax[0,0].text(ne*popt[1]+1,0.16,'{:.1f}$\,$days'.format(ne*popt[1]),color='0.3')
ax[0,0].text(3,-0.01,'N = {}'.format(ie.shape[0]),ha='center')



'Royal Charlotte Eddy'
popt,_=curve_fit(func, rce.day[~np.isnan((rce.differentiate('day')/rce).mean('eddy'))]
                        ,(rce.differentiate('day')/rce).mean('eddy')[~np.isnan((rce.differentiate('day')/rce).mean('eddy'))])

ax[0,1].plot(rce.day,func(rce.day,*popt),lw=2.5,color='0.2',zorder=1,label='exp. fit')
ax[0,1].fill_between(fit.day,func(rce.day,*popt)-fit.std('eddy'),func(rce.day,*popt)+fit.std('eddy'),color='0.9',zorder=0)
(rce.differentiate('day')/rce).mean('eddy').plot(ax=ax[0,1],marker='o',color='chocolate',lw=0,markeredgecolor='0.3',zorder=5)
ax[0,1].axvline(ne*popt[1],color='0.2',linestyle='--')#,label='e-folding')
ax[0,1].text(ne*popt[1]+1,0.16,'{:.1f}$\,$days'.format(ne*popt[1]),color='0.3')
# ax[1].scatter([],[],color='k',label='Palóczy et. al. (2013)')
# ax[1].plot(5,0.03,markersize=8,color='0.1',marker='o',markeredgecolor='w',zorder=1e5)




'Abrolhos Eddy'
popt,_=curve_fit(func, ae.day[~np.isnan((ae.differentiate('day')/ae).mean('eddy'))]
                        ,(ae.differentiate('day')/ae).mean('eddy')[~np.isnan((ae.differentiate('day')/ae).mean('eddy'))])

ax[0,2].plot(ae.day,func(ae.day,*popt),lw=2.5,color='0.2',zorder=1,label='exp. fit')
ax[0,2].fill_between(fit.day,func(ae.day,*popt)-fit.std('eddy'),func(ae.day,*popt)+fit.std('eddy'),color='0.9',zorder=0)
(ae.differentiate('day')/ae).mean('eddy').plot(ax=ax[0,2],marker='o',color='crimson',lw=0,markeredgecolor='0.3',zorder=5)
ax[0,2].axvline(ne*popt[1],color='0.2',linestyle='--')#,label='e-folding')
ax[0,2].text(ne*popt[1]+1,0.16,'{:.1f}$\,$days'.format(ne*popt[1]),color='0.3')
# ax[1].scatter([],[],color='k',label='Palóczy et. al. (2013)')
# ax[1].plot(5,0.03,markersize=8,color='0.1',marker='o',markeredgecolor='w',zorder=1e5)
ax[0,2].axhline(0.03,color='tab:red',linestyle='--',label='Palóczy et. al. (2013)')
ax[0,2].text(3,-0.01,'N = {}'.format(cste.shape[0]),ha='center')



'Vitória Eddy'
popt,_=curve_fit(func, ve.day[~np.isnan((ve.differentiate('day')/ve).mean('eddy'))]
                        ,(ve.differentiate('day')/ve).mean('eddy')[~np.isnan((ve.differentiate('day')/ve).mean('eddy'))])

ax[0,3].plot(ve.day,func(ve.day,*popt),lw=2.5,color='0.2',zorder=1,label='exp. fit')
ax[0,3].fill_between(fit.day,func(ve.day,*popt)-fit.std('eddy'),func(ve.day,*popt)+fit.std('eddy'),color='0.9',zorder=0)
(ve.differentiate('day')/ve).mean('eddy').plot(ax=ax[0,3],marker='o',color='navy',lw=0,markeredgecolor='0.3',zorder=5)
ax[0,3].axvline(ne*popt[1],color='0.2',linestyle='--')#,label='e-folding')
ax[0,3].text(ne*popt[1]+1,0.16,'{:.1f}$\,$days'.format(ne*popt[1]),color='0.3')
# ax[1].scatter([],[],color='k',label='Palóczy et. al. (2013)')
# ax[1].plot(5,0.03,markersize=8,color='0.1',marker='o',markeredgecolor='w',zorder=1e5)



'South Vitória Eddy'

popt,_=curve_fit(func, sve.day[~np.isnan((sve.differentiate('day')/sve).mean('eddy'))]
                        ,(sve.differentiate('day')/sve).mean('eddy')[~np.isnan((sve.differentiate('day')/sve).mean('eddy'))])

ax[1,0].plot(ve.day,func(ve.day,*popt),lw=2.5,color='0.2',zorder=1,label='exp. fit')
ax[1,0].fill_between(fit.day,func(ve.day,*popt)-fit.std('eddy'),func(ve.day,*popt)+fit.std('eddy'),color='0.9',zorder=0)
(ve.differentiate('day')/ve).mean('eddy').plot(ax=ax[1,0],marker='o',color='darkorange',lw=0,markeredgecolor='0.3',zorder=5)
ax[1,0].axvline(ne*popt[1],color='0.2',linestyle='--')#,label='e-folding')
ax[1,0].text(ne*popt[1]+1,0.16,'{:.1f}$\,$days'.format(ne*popt[1]),color='0.3')
# ax[1].scatter([],[],color='k',label='Palóczy et. al. (2013)')
# ax[1].plot(5,0.03,markersize=8,color='0.1',marker='o',markeredgecolor='w',zorder=1e5)



'Cabo São Tomé Eddy'
popt,_=curve_fit(func, cste.day, (cste.differentiate('day')/cfe).mean('eddy'))
ax[1,1].plot(cste.day,func(cste.day,*popt),lw=2.5,color='0.2',zorder=1,label='exp. fit')
ax[1,1].fill_between(fit.day,func(cste.day,*popt)-fit.std('eddy'),func(cste.day,*popt)+fit.std('eddy'),color='0.9',zorder=0)
(cste.differentiate('day')/cste).mean('eddy').plot(ax=ax[1,1],marker='o',color='cornflowerblue',lw=0,markeredgecolor='0.3',zorder=5)
ax[1,1].axvline(ne*popt[1],color='0.2',linestyle='--')#,label='e-folding')
ax[1,1].text(ne*popt[1]+1,0.16,'{:.1f}$\,$days'.format(ne*popt[1]),color='0.3')
# ax[1].scatter([],[],color='k',label='Palóczy et. al. (2013)')
# ax[1].plot(5,0.03,markersize=8,color='0.1',marker='o',markeredgecolor='w',zorder=1e5)
ax[1,1].axhline(0.03,color='tab:red',linestyle='--',label='Palóczy et. al. (2013)')
ax[1,1].text(3,-0.01,'N = {}'.format(cste.shape[0]),ha='center')



'Cabo Frio Eddy'
popt,_=curve_fit(func, cfe.day, (cfe.differentiate('day')/cfe).mean('eddy'))
ax[1,2].plot(cfe.day,func(cfe.day,*popt),lw=2.5,color='0.2',zorder=1)
fitstd = fit.std('eddy')
ax[1,2].fill_between(fit.day,func(cfe.day,*popt)-fitstd,func(cfe.day,*popt)+fitstd,color='0.9',zorder=0)
(cfe.differentiate('day')/cfe).mean('eddy').plot(ax=ax[1,2],marker='o',color='mediumseagreen',lw=0,markeredgecolor='0.3',zorder=5)
ax[1,2].axvline(ne*popt[1],color='0.2',linestyle='--')
ax[1,2].text(ne*popt[1]+1,0.16,'{:.1f}$\,$days'.format(ne*popt[1]),color='0.3')
ax[1,2].text(3,-0.01,'N = {}'.format(cfe.shape[0]),ha='center')





_ = [a.grid(True,linewidth=1,linestyle='--',alpha=0.5) for a in np.hstack(ax)]
_ = [a.set(xlabel='Time [days]',yticks=np.arange(0,0.23,0.03),**kw) for a in np.hstack(ax)]


ax[1,2].text(7.8,0.06,'2 e-folding',rotation=90,ha='center',color='0.4',fontsize=9)
ax[1,2].text(3,0.16,'CFE',ha='center',color='0.2',fontsize=12)
ax[1,2].set(ylabel='Growth Rate [day$^{-1}$]')


ax[1,1].text(3,0.16,'CSTE',ha='center',color='0.2',fontsize=12)
ax[1,1].set(yticklabels=[])

ax[1].legend(ncol=3,loc='upper center', bbox_to_anchor=(-0.1, 1.17),frameon=False)

# fig.savefig('img/cfe_cste.png',dpi=300,bbox_inches='tight')

'''


plt.rcParams['font.size']=15
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.weight']='bold'
import matplotlib.gridspec as gridspec


# fig,ax = plt.subplots(1,2,figsize=(6,3))
p=plt.figure(figsize=(18.53,  9.44))

gs = gridspec.GridSpec(30,40)


ax1 = plt.subplot(gs[0:13, 2:8])

'Ilheus Eddy'
# popt,_=curve_fit(func, ie.day, (ie.differentiate('day')/ie).mean('eddy'))
popt,_=curve_fit(func, ie.day[~np.isnan((ie.differentiate('day')/ie).mean('eddy'))]
                        ,(ie.differentiate('day')/ie).mean('eddy')[~np.isnan((ie.differentiate('day')/ie).mean('eddy'))])

plt.plot(ie.day,func(ie.day,*popt),lw=2.5,color='0.2',zorder=1)
fitstd = fit.std('eddy')
plt.fill_between(fit.day,func(ie.day,*popt)-fitstd,func(ie.day,*popt)+fitstd,color='0.9',zorder=0)
plt.plot((ie.differentiate('day')/ie).mean('eddy'),marker='o',color='lightsalmon',lw=0,markeredgecolor='0.3',zorder=5)
# plt.axvline(ne*popt[1],color='0.2',linestyle='--')
# plt.text(popt[0]+1,0.16,'{:.3f}$\,$days'.format(popt[0]),color='0.3')
# plt.text(popt[1]+1,0.16,'{:.2f}$\,$day $^-$$ ^1$'.format((ie.differentiate('day')/ie).mean('eddy')[np.int_(np.round(popt[1]))-1].data ),color='0.3')
plt.text(popt[1]+1,0.16,'{:.2f}$\,$day $^-$$ ^1$'.format(0.04),color='0.3')

# plt.text(2,0.001,'N = {}'.format(ie.shape[0]),ha='center')
plt.grid(linewidth=1,linestyle='--',alpha=0.5)
# plt.xlabel('Time [days]')
# plt.xticks(np.arange(0,0.23,0.03))
plt.title('IE')
plt.xlim(0,20)
plt.ylim(-0.02,0.2)
plt.ylabel(r'Growth rate [day$^{-1}$]')
# ax1.text(6.5,0.06,'2 e-folding',rotation=90,ha='center',color='0.4',fontsize=12)





ax2 = plt.subplot(gs[0:13, 12:18])

'Royal Charlotte Eddy'
popt,_=curve_fit(func, rce.day[~np.isnan((rce.differentiate('day')/rce).mean('eddy'))]
                        ,(rce.differentiate('day')/rce).mean('eddy')[~np.isnan((rce.differentiate('day')/rce).mean('eddy'))])



plt.plot(rce.day,func(rce.day,*popt),lw=2.5,color='0.2',zorder=1)
fitstd = fit.std('eddy')
plt.fill_between(fit.day,func(rce.day,*popt)-fitstd,func(rce.day,*popt)+fitstd,color='0.9',zorder=0)
plt.plot((rce.differentiate('day')/rce).mean('eddy'),marker='o',color='chocolate',lw=0,markeredgecolor='0.3',zorder=5)
# plt.axvline(ne*popt[1],color='0.2',linestyle='--')
# plt.text(popt[0]+1,0.16,'{:.3f}$\,$days'.format(popt[0]),color='0.3')
# plt.text(popt[1]+1,0.16,'{:.2f}$\,$day $^-$$ ^1$'.format((rce.differentiate('day')/rce).mean('eddy')[np.int_(np.round(popt[1]))-1].data ),color='0.3')
plt.text(popt[1]+1,0.16,'{:.2f}$\,$day $^-$$ ^1$'.format(0.05 ),color='0.3')

# plt.text(2,0.001,'N = {}'.format(rce.shape[0]),ha='center')
plt.grid(linewidth=1,linestyle='--',alpha=0.5)
# plt.xlabel('Time [days]')
plt.title('RCE')


plt.xlim(0,20)
plt.ylim(-0.02,0.2)

ax3 = plt.subplot(gs[0:13, 22:28])


'Abrolhos Eddy'
popt,_=curve_fit(func, ae.day[~np.isnan((ae.differentiate('day')/ae).mean('eddy'))]
                        ,(ae.differentiate('day')/ae).mean('eddy')[~np.isnan((ae.differentiate('day')/ae).mean('eddy'))])

plt.plot(ae.day,func(ae.day,*popt),lw=2.5,color='0.2',zorder=1)
fitstd = fit.std('eddy')
plt.fill_between(fit.day,func(ae.day,*popt)-fitstd,func(ae.day,*popt)+fitstd,color='0.9',zorder=0)
plt.plot((ae.differentiate('day')/ae).mean('eddy'),marker='o',color='crimson',lw=0,markeredgecolor='0.3',zorder=5)
# plt.axvline(ne*popt[1],color='0.2',linestyle='--')
# plt.text(popt[0]+1,0.16,'{:.3f}$\,$days'.format(popt[0]),color='0.3')
# plt.text(popt[1]+1,0.16,'{:.2f}$\,$day $^-$$ ^1$'.format((ae.differentiate('day')/ae).mean('eddy')[np.int_(np.round(popt[1]))-1].data ),color='0.3')
plt.text(popt[1]+1,0.16,'{:.2f}$\,$day $^-$$ ^1$'.format(0.06 ),color='0.3')

# plt.text(2,0.001,'N = {}'.format(ae.shape[0]),ha='center')
plt.grid(linewidth=1,linestyle='--',alpha=0.5)
# plt.xlabel('Time [days]')

plt.title('AE')

plt.xlim(0,20)
plt.ylim(-0.02,0.2)


ax4 = plt.subplot(gs[0:13, 32:38])


'Vitória Eddy'
popt,_=curve_fit(func, ve.day[~np.isnan((ve.differentiate('day')/ve).mean('eddy'))]
                        ,(ve.differentiate('day')/ve).mean('eddy')[~np.isnan((ve.differentiate('day')/ve).mean('eddy'))])

plt.plot(ve.day,func(ve.day,*popt),lw=2.5,color='0.2',zorder=1)
fitstd = fit.std('eddy')
plt.fill_between(fit.day,func(ve.day,*popt)-fitstd,func(ve.day,*popt)+fitstd,color='0.9',zorder=0)
plt.plot((ve.differentiate('day')/ve).mean('eddy'),marker='o',color='navy',lw=0,markeredgecolor='0.3',zorder=5)
# plt.axvline(ne*popt[1],color='0.2',linestyle='--')
# plt.text(popt[0]+1,0.16,'{:.3f}$\,$days'.format(popt[0]),color='0.3')
plt.text(popt[1]+1,0.16,'{:.2f}$\,$day $^-$$ ^1$'.format(0.06 ),color='0.3')

# plt.text(2,0.001,'N = {}'.format(ve.shape[0]),ha='center')
plt.grid(linewidth=1,linestyle='--',alpha=0.5)
# plt.xlabel('Time [days]')
plt.title('VE')



plt.xlim(0,20)
plt.ylim(-0.02,0.2)


ax5 = plt.subplot(gs[17:, 6:12])


'South Vitória Eddy'

popt,_=curve_fit(func, sve.day[~np.isnan((sve.differentiate('day')/sve).mean('eddy'))]
                        ,(sve.differentiate('day')/sve).mean('eddy')[~np.isnan((sve.differentiate('day')/sve).mean('eddy'))])

plt.plot(sve.day,func(sve.day,*popt),lw=2.5,color='0.2',zorder=1)
fitstd = fit.std('eddy')
plt.fill_between(fit.day,func(sve.day,*popt)-fitstd,func(sve.day,*popt)+fitstd,color='0.9',zorder=0)
plt.plot((sve.differentiate('day')/sve).mean('eddy'),marker='o',color='darkorange',lw=0,markeredgecolor='0.3',zorder=5)
# plt.axvline(ne*popt[1],color='0.2',linestyle='--')
# plt.text(popt[0]+1,0.16,'{:.3f}$\,$days'.format(popt[0]),color='0.3')
# plt.text(popt[1]+1,0.16,'{:.2f}$\,$day $^-$$ ^1$'.format((sve.differentiate('day')/sve).mean('eddy')[np.int_(np.round(popt[1]))-1].data ),color='0.3')
plt.text(popt[1]+1,0.16,'{:.2f}$\,$day $^-$$ ^1$'.format(0.08 ),color='0.3')

# plt.text(2,0.001,'N = {}'.format(sve.shape[0]),ha='center')
plt.grid(linewidth=1,linestyle='--',alpha=0.5)
plt.xlabel('Time [days]')
plt.title('SVE')
plt.ylabel(r'Growth rate [day$^{-1}$]')


plt.xlim(0,20)
plt.ylim(-0.02,0.2)

ax6 = plt.subplot(gs[17:, 16:22])


'Cabo São Tomé Eddy'
# popt,_=curve_fit(func, cste.day, (cste.differentiate('day')/cfe).mean('eddy'))

popt,_=curve_fit(func, cste.day[~np.isnan((cste.differentiate('day')/cste).mean('eddy'))]
                        ,(cste.differentiate('day')/cste).mean('eddy')[~np.isnan((cste.differentiate('day')/cste).mean('eddy'))])



plt.plot(cste.day,func(cste.day,*popt),lw=2.5,color='0.2',zorder=1)
fitstd = fit.std('eddy')
plt.fill_between(fit.day,func(cste.day,*popt)-fitstd,func(cste.day,*popt)+fitstd,color='0.9',zorder=0)
plt.plot((cste.differentiate('day')/cste).mean('eddy'),marker='o',color='cornflowerblue',lw=0,markeredgecolor='0.3',zorder=5)
# plt.axvline(ne*popt[1],color='0.2',linestyle='--')
# plt.text(popt[0]+1,0.16,'{:.3f}$\,$days'.format(popt[0]),color='0.3')
# plt.text(popt[1]+1,0.16,'{:.2f}$\,$day $^-$$ ^1$'.format((cste.differentiate('day')/cste).mean('eddy')[np.int_(np.round(popt[1]))-1].data ),color='0.3')
plt.text(popt[1]+1,0.16,'{:.2f}$\,$day $^-$$ ^1$'.format(0.08 ),color='0.3')

# plt.text(2,0.001,'N = {}'.format(cste.shape[0]),ha='center')
plt.grid(linewidth=1,linestyle='--',alpha=0.5)
plt.xlabel('Time [days]')
plt.title('CSTE')


# plt.axhline(0.03,color='tab:red',linestyle='--',label='Palóczy et. al. (2013)')
# plt.text(3,-0.01,'N = {}'.format(cste.shape[0]),ha='center')


plt.xlim(0,20)
plt.ylim(-0.02,0.2)

ax7 = plt.subplot(gs[17:, 26:32])


'Cabo Frio Eddy'
# popt,_=curve_fit(func, cfe.day, (cfe.differentiate('day')/cfe).mean('eddy'))
popt,_=curve_fit(func, cfe.day[~np.isnan((cfe.differentiate('day')/cfe).mean('eddy'))]
                        ,(cfe.differentiate('day')/cfe).mean('eddy')[~np.isnan((cfe.differentiate('day')/cfe).mean('eddy'))])



plt.plot(cfe.day,func(cfe.day,*popt),lw=2.5,color='0.2',zorder=1)
fitstd = fit.std('eddy')
plt.fill_between(fit.day,func(cfe.day,*popt)-fitstd,func(cfe.day,*popt)+fitstd,color='0.9',zorder=0)
plt.plot((cfe.differentiate('day')/cste).mean('eddy'),marker='o',color='mediumseagreen',lw=0,markeredgecolor='0.3',zorder=5)
# plt.axvline(ne*popt[1],color='0.2',linestyle='--')
# plt.text(popt[0]+1,0.16,'{:.3f}$\,$days'.format(popt[0]),color='0.3')
# plt.text(popt[1]+1,0.16,'{:.2f}$\,$day $^-$$ ^1$'.format((cfe.differentiate('day')/cfe).mean('eddy')[np.int_(np.round(popt[1]))-1].data ),color='0.3')
plt.text(popt[1]+1,0.16,'{:.2f}$\,$day $^-$$ ^1$'.format(0.07),color='0.3')

# plt.text(2,0.001,'N = {}'.format(cfe.shape[0]),ha='center')
plt.grid(linewidth=1,linestyle='--',alpha=0.5)
plt.xlabel('Time [days]')

plt.title('CFE')


plt.xlim(0,20)
plt.ylim(-0.02,0.2)

# ax8 = plt.subplot(gs[17:, 30:])





plt.tight_layout()
plt.savefig('/home/igor/Dropbox/META_eddydetection_images/growth_rate.png',bbox_inches='tight')
plt.close()





# ax[1,2].text(3,0.16,'CFE',ha='center',color='0.2',fontsize=12)
# ax[1,2].set(ylabel='Growth Rate [day$^{-1}$]')


# ax[1,1].text(3,0.16,'CSTE',ha='center',color='0.2',fontsize=12)
# ax[1,1].set(yticklabels=[])

# ax[1].legend(ncol=3,loc='upper center', bbox_to_anchor=(-0.1, 1.17),frameon=False)
