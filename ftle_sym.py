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
inttime = 1.
delta = 0.2
direction = 1

nx = 50; ny = 25
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

traj_x = np.copy(x)
traj_y = np.copy(y)

# start FTLE integration
for t in range(0, int(inttime/delta)):
    update_traj()
    show_traj()


fig  = plt.figure()
plt.contourf(x, y, vec, 10, cmap=plt.cm.Reds)
plt.show()
