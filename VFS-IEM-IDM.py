# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 17:30:11 2021

@author: 龚文武
"""

import numpy as np
np.set_printoptions(threshold=np.inf)


#--------------------------------计算灾害等级----------------------------------#
def matrix(X): 
    M = np.zeros((X.shape[0],X.shape[1]))
    for r in range(X.shape[0]):
        for l in range(X.shape[1]):
            ab_rl = X[r][l]
            a_rl = ab_rl[0]
            b_rl = ab_rl[1]
            M_rl = (X.shape[1] - (l+1))*a_rl/(X.shape[1] - 1) + ((l+1) - 1)*b_rl/(X.shape[1] - 1)
            M[r][l] = round(M_rl,1) 
    return M

def mu_a0(u,ab,cd): 
    beta = 1 
    if (u>33): 
        M_0 = matrix(ab)[0]
        ab_0 = ab[0]
        cd_0 =cd[0]
    elif (10<u<=33):
        M_0 = matrix(ab)[1]
        ab_0 = ab[1]
        cd_0 =cd[1]
    else:
        M_0 = matrix(ab)[2]
        ab_0 = ab[2]
        cd_0 =cd[2]
    mu_A = [0 for _ in range(len(M_0))]
    for l in range(len(M_0)):
        ab_l = ab_0[l] 
        cd_l = cd_0[l]
        if (cd_l[0]<=u<ab_l[0]):
            mu_l = 0.5*(1 - ((u -ab_l[0])/(cd_l[0] - ab_l[0]))**beta) 
        elif (M_0[l]<=u<ab_l[1]):
            mu_l = 0.5*(1 + ((u -ab_l[1])/(M_0[l] - ab_l[1]))**beta) 
        elif (ab_l[0]<=u<M_0[l]):
            mu_l = 0.5*(1 + ((u -ab_l[0])/(M_0[l] - ab_l[0]))**beta)
        elif (ab_l[1]<=u<cd_l[1]):
            mu_l = 0.5*(1 - ((u -ab_l[1])/(cd_l[1] - ab_l[1]))**beta)
        else:
            mu_l = 0 
        mu_A[l] = round(mu_l,3)
        #mu_A[l] = mu_l
    return mu_A

def mu_a1(u,ab,cd): 
    beta = 1
    if (u>33): 
        M_1 = matrix(ab)[0]
        ab_1 = ab[0]
        cd_1 =cd[0]
    elif (10<u<=33):
        M_1 = matrix(ab)[1]
        ab_1 = ab[1]
        cd_1 =cd[1]
    else:
        M_1 = matrix(ab)[2]
        ab_1 = ab[2]
        cd_1 =cd[2]
    #M_0 = matrix(ab)[0] 
    #M_1 = matrix(ab)[1] 
    #M_2 = matrix(ab)[2]
    mu_A = [0 for _ in range(len(M_1))]
    for l in range(len(M_1)):
        ab_l = ab_1[l]
        cd_l = cd_1[l]
        if (ab_l[0]<=u<M_1[l]):
            mu_l = 0.5*(1 + ((u -ab_l[0])/(M_1[l] - ab_l[0]))**beta)
        elif (cd_l[0]<=u<ab_l[0]):
                mu_l = 0.5*(1 - ((u -ab_l[0])/(cd_l[0] - ab_l[0]))**beta)
        elif (M_1[l]<=u<ab_l[1]):
            mu_l = 0.5*(1 + ((u -ab_l[1])/(M_1[l] - ab_l[1]))**beta)
        elif (ab_l[1]<=u<cd_l[1]):
            mu_l = 0.5*(1 - ((u -ab_l[1])/(cd_l[1] - ab_l[1]))**beta)
        else:
            mu_l = 0
        mu_A[l] = round(mu_l,3)
        #mu_A[l] = mu_l
    return mu_A

def mu_a2(u,ab,cd): 
    beta = 1
    if (u>33): 
        M_2 = matrix(ab)[0]
        ab_2 = ab[0]
        cd_2 =cd[0]
    elif (10<u<=33):
        M_2 = matrix(ab)[1]
        ab_2 = ab[1]
        cd_2 =cd[1]
    else:
        M_2 = matrix(ab)[2]
        ab_2 = ab[2]
        cd_2 =cd[2]
    #M_0 = matrix(ab)[0] 
    #M_1 = matrix(ab)[1]
    #M_2 = matrix(ab)[2] 
    mu_A = [0 for _ in range(len(M_2))]
    for l in range(len(M_2)):
        ab_l = ab_2[l]
        cd_l = cd_2[l]
        if (ab_l[0]<=u<M_2[l]):
            mu_l = 0.5*(1 + ((u -ab_l[0])/(M_2[l] - ab_l[0]))**beta) 
        elif (ab_l[1]<=u<cd_l[1]):
            mu_l = 0.5*(1 - ((u -ab_l[1])/(cd_l[1] - ab_l[1]))**beta) 
        elif (cd_l[0]<=u<ab_l[0]):
                mu_l = 0.5*(1 - ((u -ab_l[0])/(cd_l[0] - ab_l[0]))**beta)
        elif (M_2[l]<=u<ab_l[1]):
            mu_l = 0.5*(1 + ((u -ab_l[1])/(M_2[l] - ab_l[1]))**beta)
        else:
            mu_l = 0 
        mu_A[l] = round(mu_l,3)
        #mu_A[l] = mu_l
    return mu_A

def f(U): 
    f = np.zeros((U.shape[0],U.shape[1]))
    for r in range(U.shape[0]):
        U_r = sum(U[r])
        for l in range(U.shape[1]):
            u_rl = U[r][l]
            f_rl = u_rl/U_r 
            f[r][l] = f_rl
    return f
    
def ListMul(A,B,k):
    if len(A) == len(B):
        mul = map(lambda a,b: (a * b)**k, A,B)
        res = sum(mul)
        return res
    return ('之前的列表输入有误！')

def Ln(f): 
    Ln = np.zeros((f.shape[0],f.shape[1]))
    for i in range(f.shape[0]):
        for j in range(f.shape[1]):
            if (f[i][j] != 0):
                Ln_ij = np.log(f[i][j])
            else:
                Ln_ij = 0
            Ln[i][j] = Ln_ij
    return Ln

def h_list(f):
    h = np.zeros(f.shape[0])
    F = np.zeros(f.shape[0])
    lnf = Ln(f) 
    for r in range(f.shape[0]):
        F_r = ListMul(lnf[r],f[r],1)
        h_r = (-1/np.log(f.shape[1]))*F_r
        F[r] = F_r
        h[r] = h_r
    print(F) 
    return h

def weight(f): 
    h = h_list(f) 
    print(h)
    w = np.zeros(f.shape[0]) 
    for r in range(f.shape[0]):
        h_r = h[r]
        w_r = (1-h_r)/(f.shape[0]-sum(h))
        w[r] = round(w_r,2)
    #print(w)
    return w

def H_value(X,U,Y,a,b): 
    w = weight(f(U)) 
    mu = np.zeros(X.shape[1])
    if X.shape[1] == len(Y):
        for l in range(X.shape[1]):
            x_r = X[:,l]
            x = 1 - x_r          
            mu_l = (1 + (ListMul(w,x,b)/ListMul(w,x_r,b))**(a/b))**(-1) 
            mu[l] = round(mu_l, 2)
        print(mu) 
        return ListMul(mu,Y,1) 
    return ('分类有误！')

def normal_1(X,u,h): 
    mu = 0
    for i in range(X.shape[0]):
        mu = 0
        x_i = X[i]
        mu_i = np.exp((-(x_i-u)**2)/(2*h**2)) 
        mu += mu_i
    return mu

def Q_1(X,U,h):
    Q = np.zeros(U.shape[0]) 
    for j in range(U.shape[0]):
        u_j = U[j]
        q_j = normal_1(X,u_j,h)
        Q[j] = round(q_j,2) 
    return Q

#--------------------------------估计概率分布----------------------------------#

def normal_2(X,u,v,h1,h2): 
    mu = 0
    for i in range(X.shape[0]): 
        x_i = X[i]
        x1 = x_i[0]
        x2 = x_i[1]
        mu_i = np.exp(-((x1-u)**2)/(2*h1**2) - ((x2-v)**2)/(2*h2**2))
        mu += mu_i
    return mu

def Q_2(X,U,V,h1,h2): 
    Q = np.zeros((U.shape[0],V.shape[0]))
    for j in range(U.shape[0]): 
        u_j = U[j]
        for k in range(V.shape[0]):
            v_k = V[k]
            q_jk = normal_2(X,u_j,v_k,h1,h2)
            Q[j][k] = round(q_jk,2) 
    print("得到二维累计信息扩散矩阵：")
    return Q

def P_2(Q): 
    P = np.zeros((Q.shape[0],Q.shape[1]))
    for j in range(Q.shape[0]):
        for k in range(Q.shape[1]):
            q_jk = Q[j][k]
            p_jk = q_jk/sum(map(sum,Q))
            P[j][k] = round(p_jk,3) 
    print("得到致灾因子联合概率分布：")
    return P

def P_con(Q): 
    P_c = np.zeros((Q.shape[0],Q.shape[1]))
    for j in range(Q.shape[0]):
        for k in range(Q.shape[1]):
            q_jk = Q[j][k]
            p_cj = q_jk/sum(Q[j,:])
            P_c[j][k] = round(p_cj,3) 
    print("得到致灾因子条件概率分布：")
    return P_c

#-------------------------------估计脆弱性曲线---------------------------------#

def normal_3(X,u,v,o,h1,h2,h3): 
    mu = 0
    for i in range(X.shape[0]): 
        x_i = X[i]
        x1 = x_i[0]
        x2 = x_i[1]
        x3 = x_i[2]
        mu_i = np.exp(-((x1-u)**2)/(2*h1**2) - ((x2-v)**2)/(2*h2**2) - ((x3-o)**2)/(2*h3**2))
        mu += mu_i
    return mu

def Q_3(X,U,V,O,h1,h2,h3): 
    Q = np.zeros((U.shape[0],V.shape[0],O.shape[0]))
    for j in range(U.shape[0]): 
        u_j = U[j]
        for k in range(V.shape[0]):
            v_k = V[k]
            for l in range(O.shape[0]):
                o_l = O[l]
                q_jkl = normal_3(X,u_j,v_k,o_l,h1,h2,h3) 
                Q[j][k][l] = round(q_jkl,2) 
    print("得到三维累计扩散矩阵：")
    return Q

def S_max(Q): 
    S = np.zeros(Q.shape[2])
    for l in range(Q.shape[2]):
        Q_l = Q[:,:,l]
        s_l = Q_l.max()
        S[l] = s_l
    return S

def R_f(Q): 
    R = np.zeros((Q.shape[0],Q.shape[1],Q.shape[2]))
    for j in range(Q.shape[0]): 
        for k in range(Q.shape[1]):
            for l in range(Q.shape[2]):
                Q_jk = Q[:,:,l]
                q_jkl = Q_jk[j][k]
                s_l = S_max(Q)[l]
                r_jkl = q_jkl/s_l 
                R[j][k][l] = round(r_jkl,2) 
    print("得到R_f矩阵：")
    return R
    
def V_surface(X1,U,X2,V,O,R): 
    f_vul = np.zeros((R.shape[0],R.shape[1]))
    h1 = 10; h2 = 0.19 
    Q_J = Q_1(X1,U,h1)
    Q_K = Q_1(X2,V,h2)
    for j in range(Q_J.shape[0]):
        for k in range(Q_K.shape[0]):
            mu_B = np.zeros(R.shape[2])
            for l in range(R.shape[2]):
                R_l = R[:,:,l]
                r_jk = R_l[j][k]
                m1 = min(Q_J[j],r_jk)
                m2 = min(Q_K[k],r_jk)
                mu_l = max(m1,m2)
                mu_B[l] = mu_l
            f_vul_jk = ListMul(mu_B,O,1)
            f_vul[j][k] = round(f_vul_jk,2)
    print("得到的脆弱性曲线：")
    return f_vul

def ListSum(A,B): 
    if A.shape[0] == B.shape[0]:
        res = np.zeros(A.shape[0])
        for i in range(A.shape[0]):
            res[i] = ListMul(A[i],B[i],1)
        print("得到动态期望风险大小为")
        return res
    return ('矩阵维度不一样')

##-----------------------------需要输入以下参数-------------------------------##
ab = np.array([[[0,50],[50,100],[100,150],[150,250]],
              [[8,10.8],[10.8,17.2],[17.2,23.6],[23.6,30]],
              [[0,2],[2,5],[5,8],[8,11]]])
    
cd = np.array([[[0,100],[0,150],[50,250],[100,250]],
              [[8,17.2],[8,23.6],[10.8,30],[17.2,30]],
              [[0,5],[0,8],[2,11],[5,11]]]) # 数据集矩阵


Sample1 = [[176,2.72],[198,3.00],[254,3.74],[203,2.32],[251,2.49],[261,2.74],[173,1.93],
           [269,2.72],[179,2.31],[203,3.95],[226,2.56],[164,1.94],[181,1.99],[211,1.53],
           [223,2.13],[261,3.06],[197,1.83],[255,2.48],[232,2.92],[273,2.96],[211,3.68],
           [227,1.88],[287,2.28],[290,3.11],[161,3.67],[202,2.11],[232,2.46],[236,3.20],
           [242,3.03],[285,2.48],[155,2.47],[197,1.48],[220,2.45],[255,3.93],[182,1.37],
           [210,3.02],[233,2.90],[241,1.80]]

Sample2 = [[2.72,0.3819],[3.00,1.1520],[3.74,1.0750],[2.32,0.2571],
           [2.49,0.1450],[2.74,0.2983],[1.93,0.0765],[2.72,0.4013],
           [2.31,0.1895],[3.95,2.4800],[2.56,0.3648],[1.94,0.1527],
           [1.99,0.1984],[1.53,0.0452],[2.13,0.1423],[3.06,0.9351],
           [1.83,0.0841],[2.48,0.1682],[2.92,1.0410],[2.96,0.8352],
           [3.68,2.4521],[1.88,0.0251],[2.28,0.0362],[3.11,0.7341],
           [3.67,2.3580],[2.11,0.2461],[2.46,1.3100],[3.20,1.6130],
           [3.03,0.8872],[2.48,0.5902],[2.47,0.6952],[1.48,0.0267],
           [2.45,0.5241],[3.93,2.4260],[1.37,0.0528],[3.02,0.3182],
           [2.90,0.5931],[1.8,0.0725]]

Sample3 = [[176,2.72,0.6819],[198,3.00,1.3520],[254,3.74,1.3750],[203,2.32,0.2571],
           [251,2.49,0.4450],[261,2.74,0.9831],[173,1.93,0.0765],[269,2.72,0.4013],
           [179,2.31,0.2895],[203,3.95,2.3800],[226,2.56,0.7648],[164,1.94,0.1527],
           [181,1.99,0.1984],[211,1.53,0.0452],[223,2.13,0.1423],[261,3.06,1.2351],
           [197,1.83,0.0841],[255,2.48,0.7682],[232,2.92,0.7410],[273,2.96,0.8352],
           [211,3.68,2.1521],[227,1.88,0.0251],[287,2.28,0.2362],[290,3.11,0.9341],
           [161,3.67,2.0580],[202,2.11,0.2461],[232,2.46,1.3100],[236,3.20,1.6130],
           [242,3.03,1.8872],[285,2.48,0.5902],[155,2.47,0.6952],[197,1.48,0.0267],
           [220,2.45,0.5241],[255,3.93,2.2260],[182,1.37,0.0528],[210,3.02,0.8182],
           [233,2.90,0.8931],[241,1.80,0.0725]]

sample1 = [176,198,254,203,251,261,173,269,179,203,226,164,181,211,223,261,197,255,232,273,
           211,227,287,290,161,202,232,236,242,285,155,197,220,255,182,210,233,241]

sample2 = [2.72,3.74,2.32,2.49,2.74,1.93,2.72,2.31,3.95,2.56,1.94,1.99,1.53,2.13,3.06,1.83,2.48,
           2.92,2.96,3.68,1.88,2.28,3.11,3.67,2.11,2.46,3.2,3.03,2.48,2.47,1.48,2.45,3.93,1.37,3.02,2.9,1.8]

sample3 = [0.3819,1.1520,1.0750,0.2571,0.1450,0.2983,0.0765,0.4013,0.1895,2.4800,0.3648,0.1527,0.1984,0.0452,
           0.1423,0.9351,0.0841,0.1682,1.0410,0.8352,2.4521,0.0251,0.0362,0.7341,2.3580,0.2461,1.3100,1.6130,
           0.8872,0.5902,0.6952,0.0267,0.5241,2.4260,0.0528,0.3182,0.5931,0.0725]
