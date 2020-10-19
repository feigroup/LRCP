#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 17:41:11 2019

@author: catherine
"""

import numpy as np
#from sklearn.model_selection import train_test_split
import random
from interval import Interval
#m is the number of anchors
#p is the dimension of every X_ti
#Theta is a low rank matrix of size m by p
#R is the rank of Theta
#omega is the covariance matrix of X
#T is the sample size
#signal is the coefficient for error term
#mu is the mean for error term
#sigma is the variance for error term 

def simu(sigma1,m,p,R,omega,T,signal,mu,df):
    Theta  = np.random.normal(mu,sigma1,size=(m,p))
    u, s, v= np.linalg.svd(Theta, full_matrices = False)
    s[R:] = 0
    smat  = np.diag(s)
    Theta = np.dot(u, np.dot(smat, v))
    X= np.random.multivariate_normal(np.zeros(p), omega, m*T)
    #X= np.random.multivariate_normal(np.zeros(m * p), omega, T)
    X= np.reshape(X,[T,m,p])
    Xnew = np.zeros([T, m, m *p])
    for j in range(T):
        for i in range(m):
            em = np.zeros(m)
            em[i] = 1
            Xnew[j, i, :] = np.kron(X[j,i,:], np.transpose(em))
    # Xnew is n, m, m*p
    E=np.random.standard_t(df,size=T*m)
    E=(E-np.mean(E))/np.std(E)
    E=np.reshape(E,[T,m])
    Y=np.zeros([T,m])
    for j in range(T):
        Y[j,:]=np.sum(np.multiply(Theta,X[j,:,:]),axis=1)+signal*E[j,:]
    return Xnew, Y

def prsm(pX, pY, pnew, nh, m, thetaini, rhoini, alpha, beta, lambdap, epstol, maxiter, diag):
    pX = pX.reshape(nh * m, pnew, order = 'F')
    pY = pY.reshape(nh * m, 1, order = 'F')
    thetay = thetaini
    rho = rhoini
    itr = 0
    diff = 1
    while itr <= maxiter and diff > epstol:
        Xinverse = np.linalg.inv((2 * np.matmul(pX.transpose(), pX)/(nh * m) + beta * np.eye(pnew)))
        XY = beta * thetay + rho + (2 * np.matmul(pX.transpose(), pY)/(nh * m))[:, 0]
        thetax = np.matmul(Xinverse, XY)
        rhohalf = rho - alpha * beta * (thetax - thetay)
        mtheta = (thetax - rhohalf).reshape(m, int(pnew/m), order = 'F')
        u, s, v= np.linalg.svd(mtheta, full_matrices=False)
        snew = np.clip(s - lambdap/beta, 0, None)
        thetay = np.reshape(np.dot(u, np.dot(np.diag(snew), v)), pnew, order = 'F')
        rho = rhohalf - alpha * beta * (thetax - thetay)
        diff = np.sum(np.abs(thetax - thetay))
        itr = itr + 1
        #if(diag == True):
        #    print(diff)
    return thetay, pX, pY, rho

#cancpt is the candicate of change points
#Kmax is the maximum number of segmentations
#U is the optimal information criterion for different number of change points without penalty term
#taumat is the location of change points

def dynProg(cancpt, Kmax, X, Y, pnew, m, thetaini, rhoini, alpha, beta, lambdap, epstol, maxiter):  
    Nr  = Kmax - 2
    npt = cancpt.shape[0]
    V=np.zeros((npt, npt))
    V[:]=np.nan
    for j1 in range(npt):
        for j2 in range((j1+1), npt):
            if (j1==0) or (j1==(npt-1)):
                start = cancpt[j1]
            else: 
                start = cancpt[j1]+1
            end = cancpt[j2]+1
            pX1 = X[start:end, :, :]
            pY1 = Y[start:end, :]
            thetal, pX1, pY1, rhol = prsm(pX1, pY1, pnew, pY1.shape[0], m, thetaini, rhoini, alpha, beta, lambdap, epstol, maxiter, True)
            error1 = pY1[:, 0] - np.matmul(pX1, thetal)
            ssrl = np.matmul(error1.transpose(), error1)
            V[j1, j2] = ssrl
    Vt=V[0:(npt-1),1:npt]
    nw=Vt.shape[0]
    U = []#np.zeros(Kmax)
    U.append(Vt[0,nw-1])
    D = Vt[:,nw-1].copy()
    Pos = np.zeros((nw, Nr))
    Pos[:] = np.nan
    Pos[nw-1,:] = nw
    taumat = np.zeros((Nr, Nr))
    taumat[:] = np.nan
    for k in range(Nr):
        for j in range(nw-1):
            dist = Vt[j,j:(nw-1)] + D[(j+1):nw]
            D[j] = np.min(dist)
            Pos[j,0] = int(np.where(dist == np.min(dist))[0][0] + j) +1
            if k > 0:
                Pos[j,1:(k+1)] = Pos[int(Pos[j,0]),range(k)]
        U.append(D[0])
        taumat[k,range(k+1)] = Pos[0,range(k + 1)]
    return taumat, U

def evaluate(Selected_Index,true_cpt_pos):
    cpt_dist=np.zeros([len(true_cpt_pos)])
    sel_dist=np.zeros([len(Selected_Index)])
    for i in range(len(true_cpt_pos)):
        dist=abs(Selected_Index-true_cpt_pos[i])
        cpt_dist[i]=np.min(dist)
    cpt_dist_max=np.max(cpt_dist)
    for j in range(len(Selected_Index)):
        dist=abs(Selected_Index[j]-true_cpt_pos)
        sel_dist[j]=np.min(dist)
    sel_dist_max=np.max(sel_dist)
    return cpt_dist_max,sel_dist_max

    
###Simulation Start
Itermax=20
  
m=8
p=8
signal=0.5
R1=3
df=5

omega = np.zeros((p, p))
for j in range(0, p):
    for i in range(0, p):
        omega[i, j] = 0.5**(np.abs(i - j))

T1=100
T2=75
T3=150
T4=75
T5=100
nall=T1+T2+T3+T4+T5
true_cpt_pos=np.array([T1-1,T1+T2-1,T1+T2+T3-1,T1+T2+T3+T4-1])

cpt_num=np.zeros([Itermax])
d_under=np.zeros([Itermax])
d_over=np.zeros([Itermax])
MSIC=np.zeros([Itermax])
h_selected=np.zeros([Itermax])

np.random.seed(2019)

test_number=5

for Iter in range(Itermax):
    print(Iter)
    X1,Y1=simu(1,m,p,R1,omega,T1,signal,0,df)
    X2,Y2=simu(1.5,m,p,R1,omega,T2,signal,0,df)
    X3,Y3=simu(2,m,p,R1,omega,T3,signal,0,df)
    X4,Y4=simu(1,m,p,R1,omega,T4,signal,0,df)
    X5,Y5=simu(1.5,m,p,R1,omega,T5,signal,0,df)
    Xall=np.vstack((X1, X2,X3,X4,X5))
    Yall=np.vstack((Y1, Y2,Y3,Y4,Y5))
    
    hlist=range(16*int((m/np.log(nall))**(0.5)),48*int((m/np.log(nall))**(0.5)),4*int((m/np.log(nall))**(0.5)))
    hn=len(hlist)
    MSE_test=np.zeros([hn])
    test_sample=random.sample(range(0,nall),test_number)
    for kk in range(test_number):
        #X,X_test,Y,Y_test=train_test_split(Xall,Yall,test_size=0.002,random_state=0)
        X=np.delete(Xall,test_sample[kk],axis=0)
        Y=np.delete(Yall,test_sample[kk],axis=0)
        X_test=Xall[test_sample[kk],:,:]
        Y_test=Yall[test_sample[kk],:]
        
        for jj in range(hn):
            h=int(hlist[jj])
            print(h)
            pnew = X.shape[2]
            n= X.shape[0]
            S = np.zeros(n)
            Xini = X.reshape(n * m, pnew, order = 'F')
            Yini = Y.reshape(n * m, 1, order = 'F')
            XXinv = np.matmul(Xini.transpose(), Xini)
            XXinv = np.linalg.pinv(XXinv)
            thetaini = np.matmul(XXinv, np.matmul(Xini.transpose(), Yini))[:, 0]
            rhoini = np.zeros(pnew)
            lambdap =  np.sqrt((m + pnew/m)/h)* m**(-1)  * 0.3
            rhoini[:] = lambdap
            maxiter = 1000
            epstol = 1e-3
            alpha = 0.9
            beta = 1
            for j in range(h+1, n - h):
                pX1 = X[(j - h):(j + 1), :, :]
                pY1 = Y[(j - h):(j + 1), :]
                pX2 = X[(j + 1): (j + h + 1), :, :]
                pY2 = Y[(j + 1): (j + h + 1), :]
                pX = X[(j - h): (j + h + 1), :, :]
                pY = Y[(j - h): (j + h + 1),  :]
                thetal, pX1, pY1, rhol = prsm(pX1, pY1, pnew, pY1.shape[0], m, thetaini, rhoini, alpha, beta, lambdap, epstol, maxiter, True)
                thetar, pX2, pY2, rhor = prsm(pX2, pY2, pnew, pY2.shape[0], m, thetaini, rhoini, alpha, beta, lambdap, epstol, maxiter, True)
                theta, pX, pY, rho = prsm(pX, pY, pnew, pY.shape[0], m, thetaini, rhoini, alpha, beta, lambdap, epstol, maxiter, True)
                error1 = pY1[:, 0] - np.matmul(pX1, thetal)
                error2 = pY2[:, 0] - np.matmul(pX2, thetar)
                error = pY[:, 0] - np.matmul(pX, theta)
                ssrl = np.matmul(error1.transpose(), error1)/(h * m)
                ssrr = np.matmul(error2.transpose(), error2)/(h * m)
                ssr = np.matmul(error.transpose(), error)/(h * m)
                S[j] = ssr - ssrl - ssrr
                # print(j)
            cpt = np.array([0])
    
            for j in range (h + 1, n-h):
                if S[j] == np.max(S[(j - h):(j + h + 1)]):
                    cpt = np.hstack((cpt, j))
        
            cpt_new=np.hstack((cpt,(n-1)));

            taumat, U = dynProg(cpt_new, len(cpt_new),X, Y, pnew, m, thetaini, rhoini, alpha, beta, lambdap, epstol, maxiter)
            temp = np.array(U) + (np.array(range(len(U)))) * pow(n,2.0/3)
            #MSIC[Iter]=np.min(temp)
            smindex = int(np.where(temp == np.min(temp))[0][0])-1
            Selected_Index=cpt_new[taumat[smindex,range(smindex+1)].astype(int)]
            print(len(Selected_Index))
            if len(Selected_Index)==0:
                #Selected_Index=np.array[(0,(n-1))]
                #cpt_dist_max,sel_dist_max=evaluate(Selected_Index,true_cpt_pos)
                #Ymean=np.mean(Y,axis=0)
                theta_all, X_all, Y_all, rho_all = prsm(X, Y, pnew, Y.shape[0], m, thetaini, rhoini, alpha, beta, lambdap, epstol, maxiter, True)
                error = Y_test - np.matmul(X_test, theta_all)
                MSE_test[jj]=MSE_test[jj]+sum(pow(error,2))
            else:
                for i in range(len(Selected_Index)+1):
                    left=np.hstack((0,Selected_Index))
                    right=np.hstack((Selected_Index,nall))
                    if (test_sample[kk] in Interval(left[i],right[i]))==True:
                        Yinter=Y[left[i]:right[i],:]
                        Xinter=X[left[i]:right[i],:,:]
                        theta_interval, X_interval, Y_interval, rho_interval = prsm(Xinter, Yinter, pnew, Yinter.shape[0], m, thetaini, rhoini, alpha, beta, lambdap, epstol, maxiter, True)
                        error = Y_test - np.matmul(X_test, theta_interval)
                        MSE_test[jj]=MSE_test[jj]+sum(pow(error,2))
                        #Ymean=np.mean(Y[left[i]:right[i],:],axis=0)
                        #MSE_test[jj]=MSE_test[jj]+sum(pow((Y_test-Ymean),2))
            
    print(MSE_test)
    position=np.where(MSE_test==np.min(MSE_test))[0]
    position_selected=position[len(position)-1] 
    h_selected[Iter]=hlist[position_selected]
    X=Xall
    Y=Yall
    hh=int(h_selected[Iter])
    pnew = X.shape[2]
    n= X.shape[0]
    S = np.zeros(n)
    Xini = X.reshape(n * m, pnew, order = 'F')
    Yini = Y.reshape(n * m, 1, order = 'F')
    XXinv = np.matmul(Xini.transpose(), Xini)
    XXinv = np.linalg.pinv(XXinv)
    thetaini = np.matmul(XXinv, np.matmul(Xini.transpose(), Yini))[:, 0]
    rhoini = np.zeros(pnew)
    lambdap =  np.sqrt((m + pnew/m)/hh)* m**(-1)  * 0.3
    rhoini[:] = lambdap
    maxiter = 1000
    epstol = 1e-3
    alpha = 0.9
    beta = 1
    for j in range(hh+1, n - hh):
        pX1 = X[(j - hh):(j + 1), :, :]
        pY1 = Y[(j - hh):(j + 1), :]
        pX2 = X[(j + 1): (j + hh + 1), :, :]
        pY2 = Y[(j + 1): (j + hh + 1), :]
        pX = X[(j - hh): (j + hh + 1), :, :]
        pY = Y[(j - hh): (j + hh + 1),  :]
        thetal, pX1, pY1, rhol = prsm(pX1, pY1, pnew, pY1.shape[0], m, thetaini, rhoini, alpha, beta, lambdap, epstol, maxiter, True)
        thetar, pX2, pY2, rhor = prsm(pX2, pY2, pnew, pY2.shape[0], m, thetaini, rhoini, alpha, beta, lambdap, epstol, maxiter, True)
        theta, pX, pY, rho = prsm(pX, pY, pnew, pY.shape[0], m, thetaini, rhoini, alpha, beta, lambdap, epstol, maxiter, True)
        error1 = pY1[:, 0] - np.matmul(pX1, thetal)
        error2 = pY2[:, 0] - np.matmul(pX2, thetar)
        error = pY[:, 0] - np.matmul(pX, theta)
        ssrl = np.matmul(error1.transpose(), error1)/(hh * m)
        ssrr = np.matmul(error2.transpose(), error2)/(hh * m)
        ssr = np.matmul(error.transpose(), error)/(hh * m)
        S[j] = ssr - ssrl - ssrr
                # print(j)
    cpt = np.array([0])
    
    for j in range (hh + 1, n-hh):
        if S[j] == np.max(S[(j - hh):(j + hh + 1)]):
            cpt = np.hstack((cpt, j))
    cpt_new=np.hstack((cpt,(n-1)));

    taumat, U = dynProg(cpt_new, len(cpt_new),X, Y, pnew, m, thetaini, rhoini, alpha, beta, lambdap, epstol, maxiter)
    temp = np.array(U) + (np.array(range(len(U)))) * pow(n,2.0/3)
    MSIC[Iter]=np.min(temp)
    smindex = int(np.where(temp == np.min(temp))[0][0])-1
    Selected_Index=cpt_new[taumat[smindex,range(smindex+1)].astype(int)]
    print(len(Selected_Index))
    if len(Selected_Index)==0:
        Selected_Index=np.array[(0,(n-1))]
        cpt_dist_max,sel_dist_max=evaluate(Selected_Index,true_cpt_pos)
    else:
        cpt_dist_max,sel_dist_max=evaluate(Selected_Index,true_cpt_pos)
    print(cpt_dist_max)
    print(sel_dist_max)
    cpt_num[Iter]=len(Selected_Index)
    d_under[Iter]=sel_dist_max
    d_over[Iter]=cpt_dist_max

np.save('MULcovTcpt_num1',cpt_num)
np.save('MULcovTd_under1',d_under)
np.save('MULcovTd_over1',d_over)
    
        
                        
            
        
    