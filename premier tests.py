import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
import pandas as pd
import math
import sys

np.set_printoptions(threshold=sys.maxsize)

cmap_inferno = matplotlib.cm.get_cmap('inferno')

K =  1000000#facteur de convection
T0 = 381
a=0.002

Nr =  10     #nombre de mailles selon e_r
Ntheta = 20  #nombre de mailles selon e_theta
Nphi = 2    #nombre de mailles selon e_phi

rho = 8.9e3 #Kg.m^-3
Cp = 0.385e3 #J.Kg^-1
Lambda = 200

Tinfty = 450 #Température a l'infini

alpha = Lambda / (rho * Cp)

maillage = np.zeros((Nr,Ntheta,Nphi,4)) # x : ensemble des points du maillage
                                # point d'indice i selon er, j selon e_theta et k selon e_phi 
                                #   -> maillage[i][j][k]

maillage[:,:,:,3] = T0          #set every node temperature to T_0

meshes = []

def h(theta):
    return K*(1+math.cos(theta))




for i in range(Nr):
    for j in range(Ntheta):
            for k in range(Nphi):

                maillage[i,j,k,:3] = ((i)*a/(Nr), (j+0.5)*math.pi/Ntheta, 2*k*math.pi/Nphi)


def polar_to_cartesian(polar):
    cart = np.copy(polar)
    cart[:,:,:,0] = polar[:,:,:,0] * np.sin(polar[:,:,:,1])*np.cos(polar[:,:,:,2])
    cart[:,:,:,1] = polar[:,:,:,0] * np.sin(polar[:,:,:,1])*np.sin(polar[:,:,:,2])
    cart[:,:,:,2] = polar[:,:,:,0] * np.cos(polar[:,:,:,1])
    return cart


DoIprint = False




def compute_next_step(p_mesh, dt):
    Nr, Ntheta, Nphi, _ = np.shape(p_mesh)
    n_mesh = np.copy(p_mesh)

    Delta_r = a / Nr
    Delta_theta = math.pi / Ntheta
    Delta_phi = 2 * math.pi / Nphi
    for i in range(0, Nr):
        # on passe plusieurs fois sur le node 0,0,0 (le centre de la sphère)
        if i == 0:

            #r = (i) * Delta_r
            for j in range(Ntheta):
                #theta = (j+1.5) * Delta_theta
                for k in range(Nphi):
                    #phi = (k+1) * Delta_phi
                    sumtheta = 0
                    sumthetaplus = 0
                    sumthetaplusplus = 0
                    for m in range(Ntheta):
                        sumtheta += p_mesh[1,m,k,3]
                        sumthetaplus += p_mesh[2,m,k,3]
                        sumthetaplusplus += p_mesh[3,m,k,3]
                        
                        
                    #n_mesh[i,j,k,3] = 450
                    #n_mesh[i,j,k,3] = alpha*dt*(sumthetaplusplus + 2*sumthetaplus + 4*sumtheta-7*Ntheta*p_mesh[i,j,k,3])/(Delta_r**2) + p_mesh[i,j,k,3]
                    n_mesh[i,j,k,3] = alpha*dt*(sumtheta  - Ntheta*p_mesh[i,j,k,3])/(Delta_r**2) + p_mesh[i,j,k,3]
                    if k == 0 and j == 0 and DoIprint:
                        print(i,j,k)
                        print(n_mesh[i,j,k,3])



        elif i == Nr - 1:
            r = (i+1) * Delta_r
            for j in range(Ntheta):
                theta = (j+1.5) * Delta_theta
                for k in range(Nphi):
                    phi = (k+1) * Delta_phi
                    n_mesh[i,j,k,3] = (h(theta) + Lambda/Delta_r)**-1 * (h(theta) * Tinfty + Lambda * p_mesh[i-1,j,k, 3]/Delta_r)
                    if k == 0 and DoIprint:
                        print(i,j,k)
                        print(n_mesh[i,j,k,3])
        else:

            r = (i) * Delta_r
            for j in range(Ntheta):
                theta = (j+1.5) * Delta_theta
                if j == 0:
                    for k in range(Nphi):
                        phi = (k+1) * Delta_phi
                        n_mesh[i,j,k,3] = alpha * dt * ((p_mesh[i+1,j,k,3]- 2 * p_mesh[i,j,k,3]  + p_mesh[i-1,j,k,3])/(Delta_r**2) 
                                                        + 1/r * (p_mesh[i+1,j,k,3] - p_mesh[i-1,j,k,3])/Delta_r
                                                        + 1/r**2 * (p_mesh[i,(j+1)%Ntheta,k,3] - 2 * p_mesh[i,j,k,3] + p_mesh[i,j,k,3])/(Delta_theta**2)
                                                        + 1/(r**2 * math.tan(theta)) * (abs(p_mesh[i,(j+1),k,3] - p_mesh[i,j,k,3])) / (2* Delta_theta) 
                                                        ) + p_mesh[i,j,k,3]
                        if k == 0 and i == 1  and DoIprint:    
                            print(i,j,k)
                            print(n_mesh[i,j,k,3])
                            print((p_mesh[i+1,j,k,3]- 2 * p_mesh[i,j,k,3]  + p_mesh[i-1,j,k,3]))
                            print(p_mesh[i+1,j,k,3] - p_mesh[i-1,j,k,3])
                            print(p_mesh[i,(j+1)%Ntheta,k,3] - 2 * p_mesh[i,j,k,3] + p_mesh[i,j,k,3])
                            print(p_mesh[i,(j+1)%Ntheta,k,3] - p_mesh[i,j,k,3])

                elif j == Ntheta - 1:
                    for k in range(Nphi):
                        phi = (k+1) * Delta_phi
                        n_mesh[i,j,k,3] = alpha * dt * ((p_mesh[i+1,j,k,3]- 2 * p_mesh[i,j,k,3]  + p_mesh[i-1,j,k,3])/(Delta_r**2) 
                                                        + 1/r * (p_mesh[i+1,j,k,3] - p_mesh[i-1,j,k,3])/Delta_r
                                                        + 1/r**2 * (p_mesh[i,j,k,3] - 2 * p_mesh[i,j,k,3] + p_mesh[i,(j-1)%Ntheta,k,3])/(Delta_theta**2)
                                                        + 1/(r**2 * math.tan(theta)) * abs(p_mesh[i,j,k,3] - p_mesh[i,(j-1)%Ntheta,k,3]) / (2* Delta_theta) 
                                                        ) + p_mesh[i,j,k,3]
                        
                        if k == 0 and i==1 and DoIprint:
                            print(i,j,k)
                            print(n_mesh[i,j,k,3])
                        

                else:
                    for k in range(Nphi):
                        phi = (k+1) * Delta_phi
                        n_mesh[i,j,k,3] = alpha * dt * ((p_mesh[i+1,j,k,3]- 2 * p_mesh[i,j,k,3]  + p_mesh[i-1,j,k,3])/(Delta_r**2) 
                                                        + 1/r * (p_mesh[i+1,j,k,3] - p_mesh[i-1,j,k,3])/Delta_r
                                                        + 1/r**2 * (p_mesh[i,(j+1)%Ntheta,k,3] - 2 * p_mesh[i,j,k,3] + p_mesh[i,(j-1)%Ntheta,k,3])/(Delta_theta**2)
                                                        + 1/(r**2 * math.tan(theta)) * abs(p_mesh[i,(j+1)%Ntheta,k,3] - p_mesh[i,(j-1)%Ntheta,k,3]) / (2* Delta_theta) 
                                                        ) + p_mesh[i,j,k,3]
                        
                        if k == 0 and i==1 and DoIprint:
                            print(i,j,k)
                            print(n_mesh[i,j,k,3])                
                        
    return n_mesh


def simulation():
    dt = 0.000001
    meshes = []
    meshes.append(maillage)
    for i in range(100000):
        print(i)
        meshes.append(compute_next_step(meshes[i], dt))
    return meshes

simulation_meshes = simulation()


k = 0 
norm = matplotlib.colors.Normalize(vmin=380, vmax=451)
def update_graph(num):
    factor = 100
    global k
    mesh_c = polar_to_cartesian((simulation_meshes[num*factor]))
    graph._offsets3d = mesh_c[:,:,:,0].ravel(), mesh_c[:,:,:,1].ravel(), mesh_c[:,:,:,2].ravel()
    
    graph.set_facecolor(cmap_inferno(norm(mesh_c[:,:,:,3].ravel())))
        
    title.set_text('3D Test, num={}'.format(num*factor))
    return title, graph, 


fig = plt.figure()
fig.set_figheight(15)
fig.set_figwidth(15)
ax = fig.add_subplot(111, projection='3d')
ax.axes.set_xlim3d(left=-0.004, right=0.004) 
ax.axes.set_ylim3d(bottom=-0.004, top=0.004)
ax.axes.set_zlim3d(bottom=-0.004, top=0.004)


title = ax.set_title('3D Test')

mesh_c = polar_to_cartesian((simulation_meshes[99]))
print(mesh_c)
print(cmap_inferno(norm(mesh_c[:,:,:,3].ravel())))



graph = ax.scatter(mesh_c[:,:,:,0].ravel(), mesh_c[:,:,:,1].ravel(), mesh_c[:,:,:,2].ravel())
#graph = ax.scatter(mesh_c[:,:,:,0].ravel(), mesh_c[:,:,:,1].ravel(), mesh_c[:,:,:,2].ravel(),c = mesh_c[:,:,:,3], cmap= "inferno")
ani = matplotlib.animation.FuncAnimation(fig, update_graph,999, blit=False)

plt.show()
ani.save("girophare.gif", fps=1)
#print(simulation_meshes)