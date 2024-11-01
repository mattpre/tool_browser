import numpy as np


def compute_Ea(A,C,q,phi,delta,gamma,tau,z,nk):

    F = np.array([A[0],z])

    vect = np.array([-1,np.tan(tau-phi)])
    vect /= np.linalg.norm(vect)

    Ea = []
    P = []
    B = []
    G = 0
    Q = 0

    pt0 = A
    for kpt in range(len(C)):
        Eak = []
        dpt = (C[kpt]-pt0)/nk
        for kk in range(nk):
            pt00 = pt0+kk*dpt
            pt01 = pt0+(kk+1)*dpt
            area = 0.5*abs(np.cross(pt00-F,pt01-F))
            g = area*gamma
            G += g
            Q += np.cross(q[kpt],pt01-pt00)
            Bk = (G+Q)*np.array([np.cos(phi),np.sin(phi)])
            Pk = (G+Q)*np.array([(np.sin(phi)+np.cos(phi)*np.tan(tau-phi))/((pt01[1]-F[1])/(pt01[0]-F[0])+np.tan(tau-phi)),
                             (np.sin(phi)+np.cos(phi)*np.tan(tau-phi))/((pt01[1]-F[1])/(pt01[0]-F[0])+np.tan(tau-phi))*(pt01[1]-F[1])/(pt01[0]-F[0])])
            P.append(Pk)
            B.append(Bk)
            Eak.append(np.dot(Pk-Bk,vect))
        Ea.append(max(Eak))
        pt0 = C[kpt]

    return max(Ea),Ea.index(max(Ea)),P,B
