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
            if (lon>=m_left_lon and lon<=m_right_lon):
                ui.append(uu)
        else:
            if (lon>=m_left_lon and lon<=m_right_lon):
                ui.append(float('nan'))
        if (lon == 9.875):
            if (lat>=m_down_lat and lat<=m_up_lat):
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
# get coordinates of a point
def coord(time, i, j):
    # grid domain
    lats = np.linspace(m_down_lat, m_up_lat, my)
    lons = np.linspace(m_left_lon, m_right_lon, mx)

    xij = trajx[time][i,j]
    yij = trajy[time][i,j]

    yp=min(np.searchsorted(y[0], yij, side='right')-1, my-2)
    xm=np.zeros(mx)
    for k in range(0,mx):
        xm[k] = interp(y[0,yp], y[0,yp+1], \
                x[k,yp], x[k,yp+1], yij)
    xp=min(np.searchsorted(xm, xij, side='right')-1, mx-2)
    lon = interp(xm[xp], xm[xp+1], lons[xp], lons[xp+1], xij)
    lat = interp(y[0,yp], y[0,yp+1], lats[yp], lats[yp+1], yij)

    return (lon, lat)
    
# -----------------------------------------------------------------------------
# linear interpolation
def interp(x1, x2, y1, y2, xm):
    ym = y2 - (y2-y1)/(x2-x1) * (x2-xm)
    return ym

# -----------------------------------------------------------------------------
def show_traj(nlon, nlat):
    color1 = 'y'; color2= 'c'
    itraj = 0
    for i in range(0, mx):
        for j in range(0, my):
            valid = True
            if (i % nlon == 0) and (j % nlat == 0):
                itraj += 1
                itrajx = []; itrajy = []
                for t in range(0, 11):
                    coorx = coord(t, i, j)[0]
                    coory = coord(t, i, j)[1]
                    if (np.isnan(coorx) or np.isnan(coory)):
                        valid = False
                        break
                    itrajx.append(coorx)
                    itrajy.append(coory)
                if (valid):
                    color = color1 if (itraj % 2 == 0) else color2
                    m.plot(itrajx, itrajy, linewidth=2, color=color)
                    m.plot(itrajx[0], itrajy[0], 'o', mec='k', mfc=color, ms=4)
    return

o_down_lat = 25.125
o_up_lat = 74.875
o_left_lon = -89.875
o_right_lon = 9.875

# ftle data domain
down_lat = o_down_lat; up_lat = 55
left_lon = o_left_lon; right_lon = -30

# map domain
nx = 720; ny = 360
imin = 0; imax = 320; mx = imax - imin + 1
jmin = 0; jmax = 180; my = jmax - jmin + 1

pcmap = mpl.colors.LinearSegmentedColormap.from_list('pcmap',['white','blue'],16)
ncmap = mpl.colors.LinearSegmentedColormap.from_list('ncmap',['white','red'],16)

pcmap._init(); ncmap._init()
alphas = np.linspace(0,1,pcmap.N+3)

pcmap._lut[:,-1] = alphas
ncmap._lut[:,-1] = alphas

plot_type = "velocity"
#plot_type = "vorticity"

for t in range(6, 8, 4):
    t_str = str(t).zfill(4)

    # read trajectory files
    trajx = []; trajy = []
    for k in range(0, 11):
        (trajxk, nx, ny) = read('data/traj_x_pos_' + t_str \
                + '_' + str(k) + '.txt')
        (trajyk, nx, ny) = read('data/traj_y_pos_' + t_str \
                + '_' + str(k) + '.txt')
        trajxk = trajxk[imin:imax+1, jmin:jmax+1]
        trajyk = trajyk[imin:imax+1, jmin:jmax+1]
        trajx.append(trajxk)
        trajy.append(trajyk)
    # map domain
    lats = np.linspace(down_lat, up_lat, ny)
    lons = np.linspace(left_lon, right_lon, nx)
    m_left_lon = lons[imin]; m_right_lon = lons[imax]
    m_down_lat = lats[jmin]; m_up_lat = lats[jmax]
    # setup trajectory grid
    (x,y) = setup_grid(mx, my, m_down_lat, m_up_lat, m_left_lon, m_right_lon)

    fig = plt.figure(figsize=(14,7))

    m = Basemap(projection='cyl',llcrnrlat=m_down_lat,urcrnrlat=m_up_lat,\
            llcrnrlon=m_left_lon,urcrnrlon=m_right_lon,resolution='i')
    m.drawcoastlines(linewidth=.5, color='#444444')
    m.fillcontinents(color='#aaaaaa',lake_color='#dddddd')
    m.drawparallels(np.linspace(m_down_lat,m_up_lat,4),labels=[1,0,0,0], \
            fmt='%.1f')
    m.drawmeridians(np.linspace(m_left_lon,m_right_lon,6),labels=[0,0,0,1], \
            fmt='%.1f')
    m.drawmapboundary(fill_color='#dddddd')

    (nftle, nx, ny) = read('data/ftle_neg_' + t_str + '.txt')
    nlcs = get_lcs(nftle, .5)
    nlcs = nlcs[imin:imax+1, jmin:jmax+1]
    (pftle, nx, ny) = read('data/ftle_pos_' + t_str + '.txt')
    plcs = get_lcs(pftle, .5)
    plcs = plcs[imin:imax+1, jmin:jmax+1]
    
    # velocity field
    if (plot_type == 'velocity'):
        ufilename = "data/u_" + t_str + ".ascii"
        vfilename = "data/v_" + t_str + ".ascii"
        u = read_vel(ufilename)
        v = read_vel(vfilename)
        (ux, uy) = u.shape
        lons = np.linspace(m_left_lon, m_right_lon, ux)
        lats = np.linspace(m_down_lat, m_up_lat, uy)
        xx, yy = m.transform_vector(u.transpose(),v.transpose(),lons,lats,mx,my)
        
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
        lons = np.linspace(m_left_lon, m_right_lon, ux)
        lats = np.linspace(m_down_lat, m_up_lat, uy)
        oo = m.transform_scalar(omega.transpose(),lons,lats,mx,my)
        im = m.imshow(oo,plt.cm.RdBu_r)
        im.set_clim(vmin=-4e-5,vmax=4.5e-5)
        cb = m.colorbar(im, "right", size='5%', pad='2%', format='%.0e')
        cb.set_label('rad/s')

    lons = np.linspace(m_left_lon, m_right_lon, mx)
    lats = np.linspace(m_down_lat, m_up_lat, my)

    # negative LCS
    mftle = m.transform_scalar(nlcs.transpose(),lons,lats,mx,my)
    im = m.imshow(mftle,ncmap)
    #cb = m.colorbar(im, "right", size='5%', pad='2%')
    #im.set_clim(vmin=0.000004,vmax=0.00001)
    
    # postive LCS
    mftle = m.transform_scalar(plcs.transpose(),lons,lats,mx,my)
    im = m.imshow(mftle,pcmap)
    #cb = m.colorbar(im, "right", size='5%', pad='2%')
    #im.set_clim(vmin=0.000004,vmax=0.00001)

    # trajectory
    show_traj(8, 8)

    plt.title("pFTLE, nFTLE and " + plot_type + " [Day " + str(t) + "]")
    #fig.savefig('ftle_' + plot_type + '_' + t_str + '.png')
    print "day " + t_str + " FTLE and " + plot_type + " plot saved"
    plt.show()
    plt.close(fig)
