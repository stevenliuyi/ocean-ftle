import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.colors import colorConverter

def read(filename):
    f = open(filename)
    u = []
    ui = []
    i = 0
    x = 0; y = 0
    for line in f:
        i += 1
        if (i==1):
            x = int(line)
        elif (i==2):
            y = int(line)
        else:
            ui.append(float(line))

        if (i>2) and ((i-2) % y == 0):
            u.append(ui)
            ui = []
    f.close()
    return (np.array(u),x,y)

# -----------------------------------------------------------------------------
# identify LCSs
def get_lcs(ftle, percent):
    (nx,ny) = ftle.shape
    ftle_max = np.nanmax(ftle)
    lcs = np.zeros((nx,ny))
    for i in range(0,nx):
        for j in range(0,ny):
            if (ftle[i][j]>percent*ftle_max):
                lcs[i][j] = 1.
            else:
                lcs[i][j] = 0.
    return lcs

# -----------------------------------------------------------------------------
def read_vel(filename):
    f = open(filename)
    u = []
    lon = -89.875
    lat = 25.125
    ui = []
    for line in f:
        uu = float(line)
        if (uu!=-9999):
            if (lon>=left_lon and lon<=right_lon):
                ui.append(uu)
        else:
            if (lon>=left_lon and lon<=right_lon):
                ui.append(float('nan'))
        if (lon == 9.875):
            if (lat>=down_lat and lat<=up_lat):
                u.append(ui)
            ui = []
            if (lat == 74.875):
                break
            lat += .25
            lon = -89.875
        else:
            lon += .25
    f.close()
    return np.array(u).transpose()

# -----------------------------------------------------------------------------
# setup data grid
def setup_grid(nx, ny, dl, ul, ll, rl):
    # grid domain
    lats = np.linspace(dl, ul, ny)
    lons = np.linspace(ll, rl, nx)

    #latitude grid size in degree
    dlat = (ul - dl) / (ny-1.)

    # WGS84 spheroid constants
    a = 6378137             # equatorial radius (m)
    c = 6356752.3142        # polar radius (m)
    e2 = 1 - c**2/a**2      # square of eccentricity

    # compute length of .25 deg of lat and 1 deg of lon
    lat_len = dlat * (np.pi*a*(1-e2)) / (180 * (1-e2*(np.sin(lats* \
            np.pi/180))**2)**1.5)
    lon_len = (np.pi*a*np.cos(lats*np.pi/180)) / (180 * (1-e2 \
            *(np.sin(lats*np.pi/180))**2)**.5)

    lats_ext = np.arange(25.125, dl, dlat)
    lat_len_ext = dlat * (np.pi*a*(1-e2)) / (180 * (1-e2*(np.sin( \
            lats_ext* np.pi/180))**2)**1.5)

    # generate grid, (40W, 25.125N) is the point of origin
    lat_sum = [sum(lat_len_ext)]
    for i in lat_len[0:-1]:
        lat_sum.append(lat_sum[-1]+i)

    x = np.outer(lons+40, lon_len)
    y = np.outer(np.ones(nx), lat_sum)
    return (x, y)

# -----------------------------------------------------------------------------
def show_traj():
    c1 = 'y'; c2 = 'c'; c3='g'
    #m.plot(-78.75, 29.625, 'o', linewidth=15, color=c3)
    #m.plot(-79, 29.875, 'o', linewidth=15, color=c3)
    #m.plot(-79.25, 30.125, 'o',linewidth=15,  color=c3)
    #m.plot(-79.5, 30.375, 'o', linewidth=15, color=c3)
    #m.plot(-79.75, 30.625, 'o', linewidth=15, color=c3)
    #m.plot(-80.0, 30.875, 'o', linewidth=15, color=c3)

    #45N/65W
    #89/36
    m.plot([-78.75, -78.6521, -78.5575, -78.4731, -78.3934, -78.3111, -78.2250, -78.1327, -78.0289, -77.9179, -77.8023, -77.6852], [29.625, 29.5528, 29.4939, 29.4622, 29.4557, 29.4692, 29.4882, 29.5049, 29.5147, 29.5055, 29.4763, 29.4449], linewidth=5, color=c2)
    #87/38
    m.plot([-79, -78.9294, -78.8497, -78.7630, -78.6741, -78.5792, -78.4865, -78.3942, -78.2948, -78.1844, -78.0620, -77.9336], [29.875, 29.8252, 29.7550, 29.6752, 29.6090, 29.5541, 29.5240, 29.5212, 29.5424, 29.5633, 29.5735, 29.5663],  linewidth=5, color=c1)
    #85/40
    m.plot([-79.25, -79.2566, -79.2722, -79.2972, -79.3295, -79.3605, -79.3771, -79.3568, -79.2665, -79.0993, -78.8312, -78.5274], [30.125, 30.1797, 30.2432, 30.3203, 30.4168, 30.5447, 30.7230, 30.9667, 31.2765, 31.6139, 31.9382, 32.1633], linewidth=5, color=c2)
    #83/42
    m.plot([-79.5, -79.5415, -79.5338, -79.4505, -79.2972, -79.0952, -78.9010, -78.7370, -78.5919, -78.4536, -78.3142, -78.1878], [30.375, 30.6112, 30.9256, 31.2982, 31.6681, 31.9684, 32.1725, 32.3035, 32.3855, 32.4510, 32.5186, 32.5849], linewidth=5, color=c1)
    #81/44
    m.plot([-79.75, -79.7402, -79.6523, -79.5316, -79.4024, -79.2750, -79.1646, -79.0790, -79.0138, -78.9649, -78.9251, -78.8866], [30.625, 31.0146, 31.3669, 31.6596, 31.8561, 32.0116, 32.1348, 32.2320, 32.3057, 32.3515, 32.3804, 32.3995], linewidth=5, color=c2)
    #79/46
    m.plot([-80.0, -79.9551, -79.8918, -79.8294, -79.7599, -79.6767, -79.5840, -79.4831, -79.3674, -79.2518, -79.1494, -79.0635], [30.875, 31.1319, 31.3184, 31.4560, 31.5623, 31.6530, 31.7470, 31.8428, 31.9427, 32.0564, 32.1510, 32.2348], linewidth=5, color=c1)
    #77/48
    #m.plot([-80.25, -80.2907, -80.3260, -80.3566, -80.3729, -80.3916, -80.4101, -80.4366, -80.4572, -80.4696, -80.5375, -80.6021], [31.125, 31.1049, 31.0985, 31.1086, 31.1207, 31.1425, 31.1912, 31.2565, 31.2524, 31.2480, 31.2504, 31.2621], linewidth=5, color=c2)

    return

o_down_lat = 25.125
o_up_lat = 74.875
o_left_lon = -89.875
o_right_lon = 9.875

#down_lat = o_down_lat; up_lat = 45
#left_lon = o_left_lon; right_lon = -65

down_lat = o_down_lat; up_lat = 55
left_lon = o_left_lon; right_lon = -30

#down_lat = o_down_lat; up_lat = 35
#left_lon = -85; right_lon = -75

#down_lat = 30; up_lat = 35
#left_lon = -80; right_lon = -75

pcmap = mpl.colors.LinearSegmentedColormap.from_list('pcmap',['white','blue'],16)
ncmap = mpl.colors.LinearSegmentedColormap.from_list('ncmap',['white','red'],16)

pcmap._init(); ncmap._init()
alphas = np.linspace(0,1,pcmap.N+3)

pcmap._lut[:,-1] = alphas
ncmap._lut[:,-1] = alphas

plot_type = "velocity"
#plot_type = "vorticity"

for t in range(6, 8, 4):
    fig = plt.figure(figsize=(14,7))

    m = Basemap(projection='cyl',llcrnrlat=down_lat,urcrnrlat=up_lat,\
            llcrnrlon=left_lon,urcrnrlon=right_lon,resolution='i')
    m.drawcoastlines(linewidth=.5, color='#444444')
    m.fillcontinents(color='#aaaaaa',lake_color='#dddddd')
    m.drawparallels(np.linspace(down_lat,up_lat,4),labels=[1,0,0,0], \
            fmt='%.1f')
    m.drawmeridians(np.linspace(left_lon,right_lon,6),labels=[0,0,0,1], \
            fmt='%.1f')
    m.drawmapboundary(fill_color='#dddddd')

    t_str = str(t).zfill(4)

    (nftle, nx, ny) = read('data/ftle_neg_' + t_str + '.txt')
    nlcs = get_lcs(nftle, .5)
    (pftle, nx, ny) = read('data/ftle_pos_' + t_str + '.txt')
    plcs = get_lcs(pftle, .5)
    
    # velocity field
    if (plot_type == 'velocity'):
        ufilename = "data/u_" + t_str + ".ascii"
        vfilename = "data/v_" + t_str + ".ascii"
        u = read_vel(ufilename)
        v = read_vel(vfilename)
        (ux, uy) = u.shape
        lons = np.linspace(left_lon, right_lon, ux)
        lats = np.linspace(down_lat, up_lat, uy)
        xx, yy = m.transform_vector(u.transpose(),v.transpose(),lons,lats,nx,ny)
        
        #xp, yp, xq, yq = m.transform_vector(u.transpose(), \
        #       v.transpose(),lons,lats,int(ux/2),int(uy/2),returnxy=True)
        
        im = m.imshow(np.sqrt(np.square(xx)+np.square(yy)),plt.cm.Blues_r)
        im.set_clim(vmin=0,vmax=1.7)
        cb = m.colorbar(im, "right", size='5%', pad='2%')
        cb.set_label('m/s')
        #m.quiver(xq, yq, xp, yp, scale=50)
    elif (plot_type == 'vorticity'):
        ofilename = "data/omega_" + t_str + ".ascii"
        omega = read_vel(ofilename)
        (ux, uy) = omega.shape
        lons = np.linspace(left_lon, right_lon, ux)
        lats = np.linspace(down_lat, up_lat, uy)
        oo = m.transform_scalar(omega.transpose(),lons,lats,nx,ny)
        im = m.imshow(oo,plt.cm.RdBu_r)
        im.set_clim(vmin=-4e-5,vmax=4.5e-5)
        cb = m.colorbar(im, "right", size='5%', pad='2%', format='%.0e')
        cb.set_label('rad/s')

    lons = np.linspace(left_lon, right_lon, nx)
    lats = np.linspace(down_lat, up_lat, ny)

    # negative LCS
    mftle = m.transform_scalar(nlcs.transpose(),lons,lats,nx,ny)
    im = m.imshow(mftle,ncmap)
    #cb = m.colorbar(im, "right", size='5%', pad='2%')
    #im.set_clim(vmin=0.000004,vmax=0.00001)
    
    # postive LCS
    mftle = m.transform_scalar(plcs.transpose(),lons,lats,nx,ny)
    im = m.imshow(mftle,pcmap)
    #cb = m.colorbar(im, "right", size='5%', pad='2%')
    #im.set_clim(vmin=0.000004,vmax=0.00001)

    # trajectory
    #show_traj()

    plt.title("pFTLE, nFTLE and " + plot_type + " [Day " + str(t) + "]")
    #fig.savefig('ftle_' + plot_type + '_' + t_str + '.png')
    print "day " + t_str + " FTLE and " + plot_type + " plot saved"
    plt.show()
    plt.close(fig)
