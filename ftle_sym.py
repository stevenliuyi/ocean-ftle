import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

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

def calc_vec(xi, yi):
    return (u(xi, yi), v(xi, yi))


# wave parameters
A   = 50
L   = 400
sc  = 100
cx  = 20
lam = 40

nx = 200; ny = 100
x = np.linspace(0,    700, num=nx)
y = np.linspace(-250, 250, num=ny)
vec = np.zeros((nx, ny))

u, v = sym_vec()
for i in range(0, nx):
    for j in range(0, ny):
        ui, vi = calc_vec(x[i], y[j])
        vec[i,j] = np.sqrt(ui**2 + vi**2)

X, Y = np.meshgrid(x, y)
fig  = plt.figure()
plt.contourf(np.transpose(X), np.transpose(Y), vec, 10, cmap=plt.cm.Reds)
plt.show()
