import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap, cm

def read(filename):
    f = open(filename)
    u = []
    lon = -89.875
    lat = 25.125
    ui = []
    uall = []
    for line in f:
        uu = float(line)
        if (uu!=-9999):
            ui.append(uu)
        else:
            ui.append(float('nan'))
        if (lon == 9.875):
            u.append(ui)
            ui = []
            if (lat == 74.875):
                uall.append(np.array(u))
                u = []
                lon = -89.875
                lat = 25.125
                continue
            lat += .25
            lon = -89.875
        else:
            lon += .25
    f.close()
    return uall



ufilename = "u_NATL_2005.ascii"
vfilename = "v_NATL_2005.ascii"

uall = read(ufilename)
vall = read(vfilename)

ux = uall[0]
vx = vall[0]

fig = plt.figure()
lats = np.linspace(25.125, 74.875, num=200)
lons = np.linspace(-89.875, 9.875, num=400)

m = Basemap(projection='cyl',llcrnrlat=25.125,urcrnrlat=74.875,\
            llcrnrlon=-89.875,urcrnrlon=9.875,resolution='i')
m.drawcoastlines(linewidth=.5, color='#444444')
m.fillcontinents(color='#aaaaaa',lake_color='#dddddd')

m.drawparallels(np.arange(30.,75.,20.),labels=[1,0,0,0])
m.drawmeridians(np.arange(-80.,5.,20.),labels=[0,0,0,1])
m.drawmapboundary(fill_color='#dddddd')

#xp, yp, xq, yq = m.transform_vector(ux,vx,lons,lats,40,20,returnxy=True)
xx, yy = m.transform_vector(ux,vx,lons,lats,400,200)
#im = m.imshow(np.sqrt(np.square(xx)+np.square(yy)),plt.cm.RdBu_r)
im = m.imshow(yy,plt.cm.RdBu_r)
cb = m.colorbar(im, "right", size='5%', pad='2%')
cb.set_label('m/s')
#m.quiver(xq, yq, xp, yp, scale=10)

def update_img(i):
   ux = uall[i]
   vx = vall[i]
   xx, yy = m.transform_vector(ux,vx,lons,lats,400,200)
   im.set_data(np.sqrt(np.square(xx)+np.square(yy)))
   print "generating the velocity plot for day " + str(i+1)

   plt.title("velocity magnitude (Day " + str(i+1) + ")")
   return im

ani = animation.FuncAnimation(fig,update_img,frames=365)
ani.save('velocity.mp4', dpi=300, fps=30)
#plt.show()
