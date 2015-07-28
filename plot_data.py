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
            if (ftle[i,j]>percent*ftle_max):
                lcs[i,j] = 1.
            else:
                lcs[i,j] = 0.
    return (lcs, ftle_max)

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

# -----------------------------------------------------------------------------
# read trajectory files
def read_traj():
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
    return (trajx, trajy)

# -----------------------------------------------------------------------------
# identify saddle points
def saddle():
    return

# -----------------------------------------------------------------------------
# add a current time saddle point candidate
def add_saddle(x, y):
    global saddles
    dist = 4    # range of the neighbors
    for i in range(0, len(saddles)):
        neighbors = []
        for p in range(x-dist, x+dist+1):
            for q in range(y-dist, y+dist+1):
                neighbors.append((p,q))
        for k in neighbors:
            if k in saddles[i]:
                saddles[i].append((x,y))
                return
    saddles.append([(x,y)])
    return

o_down_lat = 25.125
o_up_lat = 74.875
o_left_lon = -89.875
o_right_lon = 9.875

# ftle data domain
down_lat = o_down_lat; up_lat = 55
left_lon = o_left_lon; right_lon = -30

# map domain
nx = 1200; ny = 600
imin = 0; imax = 600; mx = imax - imin + 1
jmin = 0; jmax = 300; my = jmax - jmin + 1

# LCS threshold
thres = .5

pcmap = mpl.colors.LinearSegmentedColormap.from_list('pcmap',['white','blue'],16)
ncmap = mpl.colors.LinearSegmentedColormap.from_list('ncmap',['white','red'],16)

pcmap._init(); ncmap._init()
alphas = np.linspace(0,1,pcmap.N+3)

pcmap._lut[:,-1] = alphas
ncmap._lut[:,-1] = alphas

plot_type = "velocity"
#plot_type = "vorticity"

saddle_series = []

for t in range(2874, 2875, 1):
    t_str = str(t).zfill(4)

    # read trajectory files
    # (trajx, trajy) = read_traj()
    # (x, y) = (trajx[0], trajy[0])

    # map domain
    lats = np.linspace(down_lat, up_lat, ny)
    lons = np.linspace(left_lon, right_lon, nx)
    m_left_lon = lons[imin]; m_right_lon = lons[imax]
    m_down_lat = lats[jmin]; m_up_lat = lats[jmax]

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

    (nftle, nx, ny) = read('ftle_neg_' + t_str + '.txt')
    (nlcs, nmax) = get_lcs(nftle, thres)
    nlcs = nlcs[imin:imax+1, jmin:jmax+1]
    (pftle, nx, ny) = read('ftle_pos_' + t_str + '.txt')
    (plcs, pmax) = get_lcs(pftle, thres)
    plcs = plcs[imin:imax+1, jmin:jmax+1]
    
    # velocity field
    if (plot_type == 'velocity'):
        ufilename = "u_" + t_str + ".ascii"
        vfilename = "v_" + t_str + ".ascii"
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
    mnftle = m.transform_scalar(nlcs.transpose(),lons,lats,mx,my)
    im = m.imshow(mnftle,ncmap)
    #cb = m.colorbar(im, "right", size='5%', pad='2%')
    #im.set_clim(vmin=0.000004,vmax=0.00001)
    
    # postive LCS
    mpftle = m.transform_scalar(plcs.transpose(),lons,lats,mx,my)
    im = m.imshow(mpftle,pcmap)
    #cb = m.colorbar(im, "right", size='5%', pad='2%')
    #im.set_clim(vmin=0.000004,vmax=0.00001)

    # trajectory
    # show_traj(8, 8)

    # identify saddle points
    saddles = []
    for i in range(0, len(lons)):
        for j in range(0, len(lats)):
            if (nftle[i,j] > thres*nmax) and (pftle[i,j] > thres*pmax):
                add_saddle(i,j)

    # the first element identify if the series is active (Ture)
    # or inactive (False)
    for i in range(0, len(saddle_series)):
        saddle_series[i][0] = False

    for i in range(0, len(saddles)):
        maxftle = 0; xmax = 0; ymax = 0
        for j in range(0, len(saddles[i])):
            x = saddles[i][j][0]
            y = saddles[i][j][1]
            multiply = nftle[x,y] * pftle[x,y]
            if (multiply > maxftle):
                maxftle = multiply; xmax = x; ymax = y
        # determine whether create a new saddle serires 
        # or add the saddle point to a existed one
        dist = 6
        neighbors = []
        for p in range(xmax-dist, xmax+dist+1):
            for q in range(ymax-dist, ymax+dist+1):
                neighbors.append((p,q))
        exist = False
        for l in range(0, len(saddle_series)):
            for k in neighbors:
                if (k == saddle_series[l][-1]):
                    saddle_series[l].append((xmax, ymax))
                    saddle_series[l][0] = True
                    exist = True
                    break
            if (exist): break
        if (not exist): saddle_series.append([True, (xmax, ymax)])
            
    # remove inactive series
    saddle_series_all = list(saddle_series)
    saddle_series = []
    for i in range(0, len(saddle_series_all)):
        if (saddle_series_all[i][0]):
            saddle_series.append(saddle_series_all[i])
        
    # plot saddle point series
    for i in range(0, len(saddle_series)):
        saddle_x = []
        saddle_y = []
        for j in range(1, len(saddle_series[i])):
            saddle_x.append(lons[saddle_series[i][j][0]])
            saddle_y.append(lats[saddle_series[i][j][1]])
        m.plot(saddle_x[-1], saddle_y[-1], 'o', mfc='y', ms=5)
        m.plot(saddle_x, saddle_y, '-', color='y', linewidth=2)

    plt.title("pFTLE, nFTLE and " + plot_type + " [Day " + str(t) + "]")
    fig.savefig('ftle_' + plot_type + '_' + t_str + '.png')
    print "day " + t_str + " FTLE and " + plot_type + " plot saved"
    plt.show()
    plt.close(fig)
