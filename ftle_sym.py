import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# calculate velocity formula
def sym_vec():
    x, y = sp.symbols('x y', real=True)

    phi0  = sc * lam
    k     = 2 * sp.pi / L
    yc    = A * sp.sin(k * x)
    alpha = sp.atan(A * k * sp.cos(k * x))

    phi   = phi0 * (1 - sp.tanh((y - yc) / (lam / sp.cos(alpha)))) + \
            cx * y
    phi   = -phi
    u     = sp.lambdify((x, y),  sp.diff(phi, y))
    v     = sp.lambdify((x, y), -sp.diff(phi, x))
    return (u, v)

# -----------------------------------------------------------------------------
# calculate velocity for a given point
def calc_vec(xi, yi):
    return (u(xi, yi), v(xi, yi))

# -----------------------------------------------------------------------------
# update trajectory (4-order Runge-Kutta)
def update_traj():
    global traj_x, traj_y
    for i in range(0, nx):
        for j in range(0, ny):
            ui, vi = calc_vec(traj_x[i,j], traj_y[i,j])
            (k1x, k1y) = (ui*delta*direction, vi*delta*direction)
            ui, vi = calc_vec(traj_x[i,j] + k1x/2, \
                              traj_y[i,j] + k1y/2)
            (k2x, k2y) = (ui*delta*direction, vi*delta*direction)
            ui, vi = calc_vec(traj_x[i,j] + k2x/2, \
                              traj_y[i,j] + k2y/2)
            (k3x, k3y) = (ui*delta*direction, vi*delta*direction)
            ui, vi = calc_vec(traj_x[i,j] + k3x, \
                              traj_y[i,j] + k3y)
            (k4x, k4y) = (ui*delta*direction, vi*delta*direction)
            traj_x[i,j] += k1x/6 + k2x/3 + k3x/3 + k4x/6
            traj_y[i,j] += k1y/6 + k2y/3 + k3y/3 + k4y/6
    return

# -----------------------------------------------------------------------------
# calculate FTLE field
def calc_ftle(num):
    global ftle
    for i in range(0,nx):
        for j in range(0,ny):
            # index 0:left, 1:right, 2:down, 3:up
            xt = np.zeros(4); yt = np.zeros(4)
            xo = np.zeros(4); yo = np.zeros(4)

            # central differencing except end points
            if (i==0):
                xt[0] = traj_x[i,j]; xt[1] = traj_x[i+1,j]
                yt[0] = traj_y[i,j]; yt[1] = traj_y[i+1,j]
                xo[0] = x[i,j]; xo[1] = x[i+1,j]
            elif (i==nx-1):
                xt[0] = traj_x[i-1,j]; xt[1] = traj_x[i,j]
                yt[0] = traj_y[i-1,j]; yt[1] = traj_y[i,j]
                xo[0] = x[i-1,j]; xo[1] = x[i,j]
            else:
                xt[0] = traj_x[i-1,j]; xt[1] = traj_x[i+1,j]
                yt[0] = traj_y[i-1,j]; yt[1] = traj_y[i+1,j]
                xo[0] = x[i-1,j]; xo[1] = x[i+1,j]

            if (j==0):
                xt[2] = traj_x[i,j]; xt[3] = traj_x[i,j+1]
                yt[2] = traj_y[i,j]; yt[3] = traj_y[i,j+1]
                yo[2] = y[i,j]; yo[3] = y[i,j+1]
            elif (j==ny-1):
                xt[2] = traj_x[i,j-1]; xt[3] = traj_x[i,j]
                yt[2] = traj_y[i,j-1]; yt[3] = traj_y[i,j]
                yo[2] = y[i,j-1]; yo[3] = y[i,j]
            else:
                xt[2] = traj_x[i,j-1]; xt[3] = traj_x[i,j+1]
                yt[2] = traj_y[i,j-1]; yt[3] = traj_y[i,j+1]
                yo[2] = y[i,j-1]; yo[3] = y[i,j+1]
    
            lambdas = eigs(xt, yt, xo, yo)
            if (lambdas=='nan'):
                ftle[i,j] = float('nan')
            else:
                ftle[i,j] = .5*np.log(max(lambdas))/(num*delta)
    return

# -----------------------------------------------------------------------------
# calculate eigenvalues of [dx/dx0]^T[dx/dx0]
def eigs(xt, yt, xo, yo):
    ftlemat = np.zeros((2,2))
    ftlemat[0,0] = (xt[1]-xt[0])/(xo[1]-xo[0])
    ftlemat[1,0] = (yt[1]-yt[0])/(xo[1]-xo[0])
    ftlemat[0,1] = (xt[3]-xt[2])/(yo[3]-yo[2])
    ftlemat[1,1] = (yt[3]-yt[2])/(yo[3]-yo[2])
    if (True in np.isnan(ftlemat)): return 'nan'
    ftlemat = np.dot(ftlemat.transpose(), ftlemat)
    w, v = np.linalg.eig(ftlemat)

    return w

# -----------------------------------------------------------------------------
# show trajectory
def show_traj():
    for i in range(0, nx):
        for j in range(0, ny):
            plt.plot([x[i,j], traj_x[i,j]], \
                     [y[i,j], traj_y[i,j]], color='k')
    return

# -----------------------------------------------------------------------------
# main program

# wave parameters
A   = 50
L   = 400
sc  = 100
cx  = 20
lam = 40

# integration info
inttime = 8.0
delta = 0.5
percent = 0.7

nx = 200; ny = 100
xx = np.linspace(0,    700, num=nx)
yy = np.linspace(-250, 250, num=ny)

x, y = np.meshgrid(xx, yy)
x    = np.transpose(x)
y    = np.transpose(y)

u, v = sym_vec()
vec  = np.zeros((nx, ny))
for i in range(0, nx):
    for j in range(0, ny):
        ui, vi = calc_vec(x[i,j], y[i,j])
        vec[i,j] = np.sqrt(ui**2 + vi**2)

# vec_contour = plt.contourf(x, y, vec, 20, cmap=plt.cm.Reds)
# plt.colorbar(vec_contour)

# -----------------------------------------------------------------------------
# pFTLE

traj_x = np.copy(x)
traj_y = np.copy(y)

# start positive FTLE integration
direction = 1
for t in range(0, int(inttime/delta)):
    update_traj()
    #show_traj()

# cacluate positive FTLE field
ftle = np.zeros((nx,ny))
calc_ftle(t)
plt.contourf(x, y, ftle, 100, cmap=plt.cm.Blues)
plt.show()

lcs = np.zeros((nx, ny))
for i in range(0, nx):
    for j in range(0, ny):
        ftle_max = np.nanmax(ftle)
        if ftle[i,j] > percent * ftle_max: lcs[i,j] = 1.

plt.contourf(x, y, lcs, 3, cmap=plt.cm.Blues)
plt.show()

# -----------------------------------------------------------------------------
# nFTLE

traj_x = np.copy(x)
traj_y = np.copy(y)

# start negative FTLE integration
direction = -1
for t in range(0, int(inttime/delta)):
    update_traj()
    #show_traj()

# cacluate negative FTLE field
ftle = np.zeros((nx,ny))
calc_ftle(t)
plt.contourf(x, y, ftle, 100, cmap=plt.cm.Reds)
plt.show()

lcs = np.zeros((nx, ny))
for i in range(0, nx):
    for j in range(0, ny):
        ftle_max = np.nanmax(ftle)
        if ftle[i,j] > percent * ftle_max: lcs[i,j] = 1.

plt.contourf(x, y, lcs, 3, cmap=plt.cm.Reds)
plt.show()
