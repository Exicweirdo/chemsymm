import numpy as np
import numpy.linalg as la
# this module contain several algorithms to calculate multipole moments
# the input xyzm here is a 4×n array contains x y z and mass, and charge is also compatible here
#inertia moment
def inertia(xyzm: np.ndarray, I: np.ndarray) -> int:
    for i in range(3):
        for j in range(3):
            I[i, j] = 0
    for iatom in range(xyzm.shape[1]):
        for i in range(3):
            I[i, i] += np.dot(xyzm[:3, iatom], xyzm[:3, iatom]) * xyzm[3, iatom]
            for j in range(3):
                I[i, j] -= xyzm[i, iatom]*xyzm[j, iatom] * xyzm[3, iatom]
    return 0
#mass center
def m_cent(xyzm: np.ndarray):
    return np.sum(xyzm[:3,:]*xyzm[3,:], axis = 1)/np.sum(xyzm[3,:])

# l=2,3 spherical harmonics
def Y_2(theta, phi, m):
    if abs(m) == 2:
        return 1/4 * np.sqrt(3/2) * np.sin(theta)**2 *np.exp(m*(0+1.j)*phi)
    elif abs(m) == 1:
        return 1/2 *(-m/abs(m))* np.sqrt(3/2) * np.sin(theta) * np.cos(theta) *np.exp(m*(0+1.j)*phi)
    elif m == 0:
        return 1/4 * (3*np.cos(theta)**2 - 1)
def Y_3(theta, phi, m):
    if abs(m) == 3:
        return (-m/abs(m)) * 1/8 * np.sqrt(5) * np.exp(m*(0+1.j)*phi) * np.sin(theta)**3
    elif abs(m) == 2:
        return 1/4 * np.sqrt(15/2) * np.exp(m*(0+1.j)*phi) * np.sin(theta)**2 * np.cos(theta)
    elif abs(m) == 1:
        return (-m/abs(m)) * 1/8 * np.sqrt(3) * np.exp(m*(0+1.j)*phi) * np.sin(theta) * (5*np.cos(theta)**2 - 1)
    elif m == 0:
        return 1/4 * (5*np.cos(theta)**3 - 3*np.cos(theta)) 

#gravitational quadrupole and octopole moment
#Q_lm = Σ_i r_i^l * Y_lm(θ_i, φ_i) * m_i
def G2moment(xyzm:np.ndarray):
    xyz = xyzm[:-1, :]
    rtp = np.zeros(xyz.shape)
    rtp[0,:] = la.norm(xyz, axis=0)
    rtp[1,:] = np.arctan2(la.norm(xyz[0:2,:], axis=0), xyz[2,:])
    rtp[2,:] = np.arctan2(xyz[1,:], xyz[0,:])
    moment = np.zeros(5)
    for m in range(0,3):
        z1 = np.sum(xyzm[3,:] * Y_2(rtp[1,:], rtp[2,:],  m)*rtp[0,:]**2)
        z2 = np.sum(xyzm[3,:] * Y_2(rtp[1,:], rtp[2,:], -m)*rtp[0,:]**2)
        moment[2-m] = ((z1 - z2) / 2 * (-(0+1.j)**(m-1))).real
        moment[2+m] = ((z1 + z2) / 2 * (-(0+1.j)**m)).real
    return moment

def G3moment(xyzm:np.ndarray):
    xyz = xyzm[:-1, :]
    rtp = np.zeros(xyz.shape)
    rtp[0,:] = la.norm(xyz, axis=0)
    rtp[1,:] = np.arctan2(la.norm(xyz[0:2,:], axis=0), xyz[2,:])
    rtp[2,:] = np.arctan2(xyz[1,:], xyz[0,:])
    moment = np.zeros(7)
    for m in range(0,4):
        z1 = np.sum(xyzm[3,:] * Y_3(rtp[1,:], rtp[2,:],  m)*rtp[0,:]**3)
        z2 = np.sum(xyzm[3,:] * Y_3(rtp[1,:], rtp[2,:], -m)*rtp[0,:]**3)
        moment[3-m] = ((z1 - z2) / 2 * (-(0+1.j)**(m-1))).real
        moment[3+m] = ((z1 + z2) / 2 * (-(0+1.j)**m)).real
        #moment[m+3] = np.sum(xyzm[3,:] * Y_3(rtp[1,:], rtp[2,:], m)*rtp[0,:]**3)
    return moment

# gravitational moments of pow_n basis
# example: pow_n=1 [x, y, z]; pow_n=2 [x2, y2, z2, xy, xz, yz]...
# the moment of which is eigenvalue of {Σ_n m_i^n t_i^n t_j^n}_ij of basis {t_i}
# due to its rank, it has C(2+pow_n, 2) non-zero eigenvalue
def multipole(xyzm:np.ndarray, pow_n:int):
    xyz = xyzm[:3,:]
    m_half = (xyzm[3,:]**0.5)
    pole_half = m_half.reshape(1,-1)
    for i in range(pow_n):
        pole_half = np.einsum("ik,jk->ijk", pole_half, xyz).reshape(3**(i+1),-1)
    pole_mat = pole_half @ pole_half.T
    return la.eigh(pole_mat)[0][::-1][:(pow_n+1)*(pow_n+2)//2]