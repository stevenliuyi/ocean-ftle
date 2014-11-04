import numpy as np

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
def read_vel(filename):
    f = open(filename)
    u = []
    lon = -89.875
    lat = 25.125
    ui = []
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
                break
            lat += .25
            lon = -89.875
        else:
            lon += .25
    f.close()
    return np.array(u).transpose()


(nx, ny) = (400, 200)

o_down_lat = 25.125         # original domain
o_up_lat = 74.875
o_left_lon = -89.875
o_right_lon = 9.875
dlat = (o_up_lat - o_down_lat) / (ny-1)
dlon = (o_right_lon - o_left_lon) / (nx-1)

# grid data
(x ,y) = setup_grid(nx, ny, o_down_lat, o_up_lat, \
    o_left_lon, o_right_lon)

for t in range(8, 156, 4):
    # velocity field
    t_str = str(t).zfill(3)
    ufilename = "data/u_" + t_str + ".ascii"
    vfilename = "data/v_" + t_str + ".ascii"
    u = read_vel(ufilename)
    v = read_vel(vfilename)

    filename = "data/omega_" + t_str + ".ascii"
    f = open(filename, 'w')
    
    # calculate dudy, dvdx
    for j in range(0, ny):
        for i in range(0, nx):
            if (i == 0):
                dvdx = (v[i+1][j] - v[i][j]) / (x[i+1][j] - x[i][j])
            elif (i == nx-1):
                dvdx = (v[i][j] - v[i-1][j]) / (x[i][j] - x[i-1][j])
            else:
                dvdx = (v[i+1][j] - v[i-1][j]) / (x[i+1][j] - x[i-1][j])

            if (j == 0):
                dudy = (u[i][j+1] - u[i][j]) / (y[i][j+1] - y[i][j])
            elif (j == ny-1):
                dudy = (u[i][j] - u[i][j-1]) / (y[i][j] - y[i][j-1])
            else:
                dudy = (u[i][j+1] - u[i][j-1]) / (y[i][j+1] - y[i][j-1])

            omega = dvdx - dudy
            f.write(str(omega) + '\n')

            

