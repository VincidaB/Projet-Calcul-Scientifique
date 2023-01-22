from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation


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

norm = matplotlib.colors.Normalize(vmin=380, vmax=451)

title = ax.set_title('Simulation 3D')

mesh_c = polar_to_cartesian((simulation_meshes[0]))
graph = ax.scatter(mesh_c[:,:,:,0].ravel(), mesh_c[:,:,:,1].ravel(), mesh_c[:,:,:,2].ravel())
ani = matplotlib.animation.FuncAnimation(fig, update_graph,249, blit=False)

plt.show()
ani.save("simulation.gif", fps=1)