import sys
import time
import numpy as np

# -----------------------------------------------------------------------------
# arguments: t_start, direction, scale, delta, inttime, method

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
# read velocity from data file
def read_vel(filename):
    f = open(filename)
    u = []
    lon = -89.875
    lat = 25.125
    ui = []
    for line in f:
        uu = float(line)
        if (uu != -9999):
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

# -----------------------------------------------------------------------------
# read u, v data at given time
def read_data(t):
    ufile = 'u_' + str(t).zfill(4) + '.ascii'
    vfile = 'v_' + str(t).zfill(4) + '.ascii'
    return (read_vel(ufile), read_vel(vfile))

# -----------------------------------------------------------------------------
# linear interpolation
def interp(x1, x2, y1, y2, xm):
    ym = y2 - (y2-y1)/(x2-x1) * (x2-xm)
    return ym

# -----------------------------------------------------------------------------
# interpolate velocity in space
def interp_vel(velt, xij, yij):
    (nx, ny) = x_data.shape
    yp=min(np.searchsorted(y_data[0], yij)-1, ny-2)
    xm=np.zeros(nx)
    for i in range(0,nx):
        xm[i] = interp(y_data[0,yp], y_data[0,yp+1], \
                x_data[i,yp], x_data[i,yp+1], yij)
    xp=min(np.searchsorted(xm, xij)-1, nx-2)

    u_left = interp(x_data[xp,yp], x_data[xp,yp+1], \
            velt[xp][yp][0], velt[xp][yp+1][0], xm[xp])
    u_right = interp(x_data[xp+1,yp], x_data[xp+1,yp+1], \
            velt[xp+1][yp][0], velt[xp+1][yp+1][0], xm[xp+1])
    um = interp(xm[xp], xm[xp+1], u_left, u_right, xij)

    v_left = interp(x_data[xp,yp], x_data[xp,yp+1], \
            velt[xp][yp][1], velt[xp][yp+1][1], xm[xp])
    v_right = interp(x_data[xp+1,yp], x_data[xp+1,yp+1], \
            velt[xp+1][yp][1], velt[xp+1][yp+1][1], xm[xp+1])
    vm = interp(xm[xp], xm[xp+1], v_left, v_right, xij)

    return (um,vm)

# -----------------------------------------------------------------------------
# set velocity in the trajectory grid
def set_vel(tt):
    # interpolate velocity in time
    (nx, ny) = u.shape
    velt = np.zeros((nx,ny,2))
    for i in range(0,nx):
        for j in range(0,ny):
            velt[i][j][0] = interp(spd*t_jump, \
                    spd*(t_jump+direction), u1[i][j], u[i][j], tt*delta*direction)
            velt[i][j][1] = interp(spd*t_jump, \
                    spd*(t_jump+direction), v1[i][j], v[i][j], tt*delta*direction)

    # interpolate velocity in space
    #(ox, oy) = traj_x.shape
    #vel = np.zeros((ox,oy,2))
    #for i in range(0,ox):
    #    for j in range(0,oy):
    #        (vel[i][j][0], vel[i][j][1]) = interp_vel(velt, \
    #                x_data, y_data, traj_x[i][j], traj_y[i][j])
    return velt

# -----------------------------------------------------------------------------
# update trajectory
def update_traj(method):
    global traj_x, traj_y
    for i in range(0,ox):
        for j in range(0,oy):
            (vx, vy) = interp_vel(velt, traj_x[i][j], traj_y[i][j])
            if (method == 1):
                traj_x[i][j] += vx*delta*direction
                traj_y[i][j] += vy*delta*direction
            elif (method == 2):
                k1x = vx*delta*direction
                k1y = vy*delta*direction
                (vx1, vy1) = interp_vel(velt, traj_x[i][j]+k1x/2, \
                        traj_y[i][j]+k1y/2)
                (vx2, vy2) = interp_vel(velt2, traj_x[i][j]+k1x/2, \
                        traj_y[i][j]+k1y/2)
                k2x = (vx1+vx2)/2*delta*direction
                k2y = (vy1+vy2)/2*delta*direction
                traj_x[i][j] += k2x
                traj_y[i][j] += k2y
            elif (method == 3):
                k1x = vx*delta*direction
                k1y = vy*delta*direction
                (vx, vy) = interp_vel(velt2, traj_x[i][j]+k1x/2, \
                        traj_y[i][j]+k1y/2)
                k2x = vx*delta*direction
                k2y = vy*delta*direction
                (vx, vy) = interp_vel(velt2, traj_x[i][j]+k2x/2, \
                        traj_y[i][j]+k2y/2)
                k3x = vx*delta*direction
                k3y = vy*delta*direction
                (vx, vy) = interp_vel(velt2, traj_x[i][j]+k3x, \
                        traj_y[i][j]+k3y)
                k4x = vx*delta*direction
                k4y = vy*delta*direction
                traj_x[i][j] += k1x/6 + k2x/3 + k3x/3 + k4x/6
                traj_y[i][j] += k1y/6 + k2y/3 + k3y/3 + k4y/6

    return

# -----------------------------------------------------------------------------
# calculate FTLE field
def calc_ftle(num):
    global ftle
    for i in range(0,ox):
        for j in range(0,oy):
            # index 0:left, 1:right, 2:down, 3:up
            xt = np.zeros(4); yt = np.zeros(4)
            xo = np.zeros(4); yo = np.zeros(4)

            # central differencing except end points
            if (i==0):
                xt[0] = traj_x[i][j]; xt[1] = traj_x[i+1][j]
                yt[0] = traj_y[i][j]; yt[1] = traj_y[i+1][j]
                xo[0] = x[i][j]; xo[1] = x[i+1][j]
            elif (i==ox-1):
                xt[0] = traj_x[i-1][j]; xt[1] = traj_x[i][j]
                yt[0] = traj_y[i-1][j]; yt[1] = traj_y[i][j]
                xo[0] = x[i-1][j]; xo[1] = x[i][j]
            else:
                xt[0] = traj_x[i-1][j]; xt[1] = traj_x[i+1][j]
                yt[0] = traj_y[i-1][j]; yt[1] = traj_y[i+1][j]
                xo[0] = x[i-1][j]; xo[1] = x[i+1][j]

            if (j==0):
                xt[2] = traj_x[i][j]; xt[3] = traj_x[i][j+1]
                yt[2] = traj_y[i][j]; yt[3] = traj_y[i][j+1]
                yo[2] = y[i][j]; yo[3] = y[i][j+1]
            elif (j==oy-1):
                xt[2] = traj_x[i][j-1]; xt[3] = traj_x[i][j]
                yt[2] = traj_y[i][j-1]; yt[3] = traj_y[i][j]
                yo[2] = y[i][j-1]; yo[3] = y[i][j]
            else:
                xt[2] = traj_x[i][j-1]; xt[3] = traj_x[i][j+1]
                yt[2] = traj_y[i][j-1]; yt[3] = traj_y[i][j+1]
                yo[2] = y[i][j-1]; yo[3] = y[i][j+1]
    
            lambdas = eigs(xt, yt, xo, yo)
            if (lambdas=='nan'):
                ftle[i][j] = float('nan')
            else:
                ftle[i][j] = .5*np.log(max(lambdas))/(num*delta)
    return

# -----------------------------------------------------------------------------
# calculate eigenvalues of [dx/dx0]^T[dx/dx0]
def eigs(xt, yt, xo, yo):
    ftlemat = np.zeros((2,2))
    ftlemat[0][0] = (xt[1]-xt[0])/(xo[1]-xo[0])
    ftlemat[1][0] = (yt[1]-yt[0])/(xo[1]-xo[0])
    ftlemat[0][1] = (xt[3]-xt[2])/(yo[3]-yo[2])
    ftlemat[1][1] = (yt[3]-yt[2])/(yo[3]-yo[2])
    if (True in np.isnan(ftlemat)): return 'nan'
    ftlemat = np.dot(ftlemat.transpose(), ftlemat)
    w, v = np.linalg.eig(ftlemat)

    return w

# -----------------------------------------------------------------------------
# identify LCSs
def get_lcs(percent):
    global lcs
    ftle_max = np.nanmax(ftle)
    for i in range(0,ox):
        for j in range(0,oy):
            if (ftle[i][j]>percent*ftle_max):
                lcs[i][j] = 1.
            else:
                lcs[i][j] = 0.
    return 

# -----------------------------------------------------------------------------
# write data to files
def write_data(data, filename):
    (ix, iy) = data.shape
    f = open(filename, 'w')
    f.write(str(ix) + '\n')
    f.write(str(iy) + '\n')
    for i in range(0,ix):
        for j in range(0,iy):
            f.write(str(data[i][j]) + '\n')
    f.close()
    return
# -----------------------------------------------------------------------------
# get longitude and latitude form indeces
def get_coor(trajx, trajy):
    yp = min(np.searchsorted(y[0], trajy)-1, oy-2)
    xm=np.zeros(ox)
    for i in range(0,ox):
        xm[i] = interp(y[0,yp], y[0,yp+1], \
                x[i,yp], x[i,yp+1], trajy)
    xp = min(np.searchsorted(xm, trajx)-1, ox-2)
    ii = interp(xm[xp], xm[xp+1], xp, xp+1, trajx)
    jj = interp(y[0,yp],y[0,yp+1], yp, yp+1, trajy)

    return (left_lon+ii*(right_lon-left_lon)/(ox-1.), \
            down_lat+jj*(up_lat-down_lat)/(oy-1.))

# -----------------------------------------------------------------------------
# grid info
(nx, ny) = (400, 200)       # number of grids in one direction
o_down_lat = 25.125         # original domain
o_up_lat = 74.875
o_left_lon = -89.875
o_right_lon = 9.875
dlat = (o_up_lat - o_down_lat) / (ny-1)
dlon = (o_right_lon - o_left_lon) / (nx-1)
scale = (float(sys.argv[3]), float(sys.argv[3]))

# part of domain need to be computed
down_lat = o_down_lat; up_lat = o_up_lat
left_lon = o_left_lon; right_lon = o_right_lon

#down_lat = o_down_lat; up_lat = 45
#left_lon = o_left_lon; right_lon = -65

down_lat = o_down_lat; up_lat = 55
left_lon = o_left_lon; right_lon = -30

#down_lat = o_down_lat; up_lat = 35
#left_lon = -85; right_lon = -75

#down_lat = 30; up_lat = 35
#left_lon = -80; right_lon = -75

# integration info
spd = 86400                 # seconds per day
inttime = spd*(float(sys.argv[5]))
delta = spd*float(sys.argv[4])
direction = int(sys.argv[2])
method = int(sys.argv[6])

# -----------------------------------------------------------------------------
# setup data grid
(x_data,y_data) = setup_grid(nx, ny, o_down_lat, o_up_lat, \
        o_left_lon, o_right_lon)

# setup trajectory grid
ox = int(int((right_lon-left_lon)/dlon+1)*scale[0])
oy = int(int((up_lat-down_lat)/dlat+1)*scale[1])
(x,y) = setup_grid(ox, oy, down_lat, up_lat, left_lon, right_lon)

t_start = int(sys.argv[1])

t_jump = 0
(u1,v1) = read_data(t_start+t_jump)
(u,v) = read_data(t_start+t_jump+direction)

traj_x = np.copy(x)
traj_y = np.copy(y)

# initialize FTLE field
ftle = np.zeros((ox,oy))

#show_traj(5,5,'r')
start_time = time.time()
#print get_coor(traj_x[77][48], traj_y[77][48])

dir_str = 'pos' if direction==1 else 'neg'
# write the trajectory data at the initial moment
write_data(x, 'traj_x_' + dir_str + '_' + str(t_start).zfill(4) + '_0.txt')
write_data(y, 'traj_y_' + dir_str + '_' + str(t_start).zfill(4) + '_0.txt')

# start FTLE integration
for t in range(0, int(inttime/delta)):
    print 'time: ' + str(t_start) + ' (' + dir_str + \
          '), integration time: ' + str(t+1) + \
          '/' + str(int(inttime/delta))
    
    if (t == 0):
        if (abs(t*delta) > abs(spd*(t_jump+direction))):
            t_jump += direction
            (u1,v1) = (u,v)
            (u,v) = read_data(t_start+t_jump+direction)
        # interpolate velocity in time and space
        velt = set_vel(t)
    else:
        velt =velt2

    if (abs((t+1)*delta) > abs(spd*(t_jump+direction))):
        t_jump += direction
        (u1,v1) = (u,v)
        (u,v) = read_data(t_start+t_jump+direction)
    velt2 = set_vel(t+1)

    # use updated velocity to compute trajectory 
    update_traj(method)
    #show_traj(5,5,'b')
    #print get_coor(traj_x[77][48], traj_y[77][48])
    write_data(traj_x, 'traj_x_' + dir_str + '_' + str(t_start).zfill(4) \
            + '_' + str(t+1) + '.txt')
    write_data(traj_y, 'traj_y_' + dir_str + '_' + str(t_start).zfill(4) \
            + '_' + str(t+1) + '.txt')

    print 'calculation time: ' + str(time.time()-start_time)

# cacluate FTLE field
calc_ftle(t)

# write the FTLE data and the trajectory data
write_data(ftle, 'ftle_' + dir_str + '_' + \
           str(t_start).zfill(4) + '.txt')
print dir_str + ' FTLE field at t=' + str(t_start) + \
      ' has been calculated'

