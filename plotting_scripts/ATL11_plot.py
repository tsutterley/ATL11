# -*- coding: utf-8 -*-
"""
Created on Wed May 16 21:05:50 2018

@author: ben
"""
import matplotlib.pyplot as plt
import numpy as np

class ATL11_plot:
    def __init__(self, D11):
        self.D11=D11

        fig=plt.figure()

        self.ax2=plt.axes([0.05,  0.05, 0.4, 0.375])
        self.ax3=plt.axes([0.55,  0.05, 0.4, 0.375])
        self.ax1=plt.axes([0.05, 0.525, 0.9, 0.375])

        self.h_errorbars=D11.plot()

        #fig.canvas.mpl_connect('pick_event', self.pick_event)
        fig.canvas.mpl_connect('button_press_event', self.pick_event)

        plt.show(block=True)
        return

    def pick_event(self, event):
        ONEYEAR=24*3600*365.25
        t0=1.5*ONEYEAR
        #xx=event.artist.get_xdata()
        xx=event.xdata
        #ii=event.ind
        D11=self.D11
        x0=D11.ref_surf.x_atc
        this=np.argmin(np.abs(x0-xx))
        #this=np.flatnonzero(x0==xx[ii])
        print("this is {0:d}".format(this))
        print("x0={0:d}".format(x0[this]))
        #inc are the selected cycles
        inc=np.flatnonzero(D11.cycle_stats.n_segs[this,:])
        #print "inc:"
        #print inc
        plt.sca(self.ax2)
        plt.cla()
        # plot the heights and
        yy=D11.cycle_stats.y_atc[this,:].ravel()
        hh=D11.corrected_h.h_corr[this,:].ravel()
        ss=D11.corrected_h.h_corr_sigma[this,:].ravel()

        hc=D11.corrected_h.h_corr[this,:]
        hc_sigma=D11.corrected_h.h_corr_sigma[this,:]
        t=D11.corrected_h.delta_time[this,:]
        good=np.logical_and(np.logical_and(hc_sigma<10 , np.isfinite(hc_sigma)), np.isfinite(t))
        G=np.ones((good.sum(), 2))
        G[:,1]=(t[good]-t0)/ONEYEAR
        cov_data=np.diag(hc_sigma[good]**2)
        cov_data_i=np.diag(1/hc_sigma[good]**2)
        GcG=(G.transpose().dot(cov_data_i).dot(G))
        Ginv=np.linalg.solve(GcG, G.transpose().dot(cov_data_i))
        m=Ginv.dot(hc[good])
        sigma_m=np.sqrt(np.diagonal(Ginv.dot(cov_data).dot(Ginv.transpose())))

        plt.errorbar(yy, hh, ss, fmt='o')
        plt.plot(yy, D11.cycle_stats.h_mean[this,:].ravel(), 'kx')
        plt.plot(yy[inc], hh[inc],'r*', markersize =12)
        plt.sca(self.ax3)
        plt.cla()
        #cycles=np.arange(yy.size)
        plt.errorbar(t.ravel()/ONEYEAR, hh.ravel(), ss.ravel(), fmt='o')
        plt.plot(t[inc]/ONEYEAR, hh[inc],'r*', markersize=12)
        tt=np.linspace(0, 3, 40)*ONEYEAR
        plt.plot(tt/ONEYEAR, m[0]+m[1]*(tt-t0)/ONEYEAR,'k-')
        eZ=np.sqrt(sigma_m[0]**2+(sigma_m[1]*(tt-t0)/ONEYEAR)**2)
        plt.plot(tt/ONEYEAR, m[0]+m[1]*(tt-t0)/ONEYEAR+eZ,'k--' )
        plt.plot(tt/ONEYEAR, m[0]+m[1]*(tt-t0)/ONEYEAR-eZ,'k--' )
        plt.show(block=False)
        plt.title('dhdt={0:f}'.format(m[1]))
        return
