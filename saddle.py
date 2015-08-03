import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import time

# -----------------------------------------------------------------------------
# read data
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

folder = '/Volumes/TOSHIBA EXT/Ocean/'


o_down_lat = 25.125
o_up_lat = 74.875
o_left_lon = -89.875
o_right_lon = 9.875

# ftle data domain
down_lat = o_down_lat; up_lat = 55
left_lon = o_left_lon; right_lon = -30

# map domain
nx = 1200; ny = 600
imin = 0; imax = 1199; mx = imax - imin + 1
jmin = 0; jmax = 599; my = jmax - jmin + 1

# LCS threshold
thres = .5

# map domain
lats = np.linspace(down_lat, up_lat, ny)
lons = np.linspace(left_lon, right_lon, nx)
m_left_lon = lons[imin]; m_right_lon = lons[imax]
m_down_lat = lats[jmin]; m_up_lat = lats[jmax]

read_data = True
save_data = True

if read_data:
    saddle_series_archive = list(np.load('saddle_series_archive.npy'))
    saddle_series = list(np.load('saddle_series.npy'))
else:
    saddle_series_archive = []
    saddle_series = []

start_time = time.time()
for t in range(2000, 2500, 1):
    t_str = str(t).zfill(4)

    (nftle, nx, ny) = read(folder + 'FTLE/ftle_neg_' + t_str + '.txt')
    (nlcs, nmax) = get_lcs(nftle, thres)
    nlcs = nlcs[imin:imax+1, jmin:jmax+1]
    (pftle, nx, ny) = read(folder + 'FTLE/ftle_pos_' + t_str + '.txt')
    (plcs, pmax) = get_lcs(pftle, thres)
    plcs = plcs[imin:imax+1, jmin:jmax+1]

    # identify saddle points
    saddles = []
    for i in range(0, len(lons)):
        for j in range(0, len(lats)):
            if (nftle[i,j] > thres*nmax) and (pftle[i,j] > thres*pmax):
                add_saddle(i,j)

    # the first element identify if the series is active (True)
    # or inactive (False)
    # set default value to be inactive
    for i in range(0, len(saddle_series)):
        saddle_series[i][0] = False

    for i in range(0, len(saddles)):
        # find the repsentative sadde point in a neighborhood
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
            
    # archive inactive series
    saddle_series_all = list(saddle_series)
    saddle_series = []
    for i in range(0, len(saddle_series_all)):
        if (saddle_series_all[i][0]):
            saddle_series.append(saddle_series_all[i])
        else:
            if len(saddle_series_all[i]) > 15:
                saddle_series_archive.append(saddle_series_all[i])
    
    print 'day ' + t_str + ' calculated (time: ' + str(time.time() - start_time) + ')'
    if save_data:
        np.save('saddle_series.npy', np.array(saddle_series))
        np.save('saddle_series_archive.npy', np.array(saddle_series_archive))

# plot saddle point series
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

for i in range(0, len(saddle_series_archive)):
    if len(saddle_series_archive[i]) <= 60: continue
    saddle_x = []
    saddle_y = []
    for j in range(1, len(saddle_series_archive[i])):
        saddle_x.append(lons[saddle_series_archive[i][j][0]])
        saddle_y.append(lats[saddle_series_archive[i][j][1]])
    m.plot(saddle_x[0], saddle_y[0], 'o', mfc='r', ms=4)
    m.plot(saddle_x, saddle_y, '-', color='r', linewidth=1)

plt.title("saddle points")
fig.savefig('saddle.png')
#plt.show()
plt.close(fig)

