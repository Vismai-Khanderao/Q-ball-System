import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.gridspec import GridSpec
from scipy.integrate import odeint


def qball_maker(w, pos, v, data):
    v_x, v_y = v
    speed = (v_x**2+v_y**2)**0.5
    gamma = 1/(1-speed**2)**0.5

    x, y = np.meshgrid(np.linspace(-domain, domain, n), np.linspace(-domain, domain, n))
    x0, y0 = pos

    # rotating so that qball is travelling in the x direction for lorentz transformation
    theta = -np.arctan2(v_y,v_x)

    x_rotated = x * np.cos(theta) - y * np.sin(theta)
    y_rotated = x * np.sin(theta) + y * np.cos(theta)
    x0_rotated = x0 * np.cos(theta) - y0 * np.sin(theta)
    y0_rotated = x0 * np.sin(theta) + y0 * np.cos(theta)

    gamma_x = (x_rotated-x0_rotated) * gamma
    boosted_coords = (gamma_x**2 + (y_rotated-y0_rotated)**2)**0.5    

    sigma_gamma = np.array([[lag(i, data) for i in row] for row in boosted_coords])

    phi_01_boost = sigma_gamma*np.cos(w*gamma*speed*(x_rotated-x0_rotated))
    phi_02_boost = sigma_gamma*np.sin(w*gamma*speed*(x_rotated-x0_rotated))

    dsdx = np.gradient(sigma_gamma, np.linspace(-domain, domain, n), axis=1)
    dsdy = np.gradient(sigma_gamma, np.linspace(-domain, domain, n), axis=0)
    dsdu = (dsdx*v_x/speed + dsdy*v_y/speed)/gamma

    pi_01_boost = dsdu*np.cos(w*gamma*speed*(x_rotated-x0_rotated))*gamma*speed - w*phi_02_boost*gamma
    pi_02_boost = dsdu*np.sin(w*gamma*speed*(x_rotated-x0_rotated))*gamma*speed + w*phi_01_boost*gamma

    return phi_01_boost, phi_02_boost, pi_01_boost, pi_02_boost

def lag(x, data):
    i = (np.abs(data[:,0] - x)).argmin()     

    low_i = max(0, i-2)
    high_i = min(len(data[:,0]), low_i+4)

    if high_i==len(data[:,0]):
        low_i = len(data[:,0])-4

    xdata = [data[:,0][n] for n in range(low_i, high_i)]
    ydata = [data[:,1][n] for n in range(low_i, high_i)]

    data_points = len(xdata)  # 4 in this case

    # main lagrange interpolation function
    y = 0
    for i in range(data_points):
        p = 1        
        for j in range(data_points):
            if i != j:
                p = p * (x - xdata[j])/(xdata[i] - xdata[j])
        y = y + p * ydata[i]   

    return y

def ode_solver(s0, bounds):
    w_high, w_low = bounds
    w = (w_high+w_low)/2

    r0 = [s0,0]
    
    r = np.linspace(0, 100, int(1e3))
    y_x = odeint(qball_ode, r0, r, (w,))

    while abs(w_high-w_low)>1e-15:
        inf = 0  

        for y in y_x[:,0]:                                   
            if abs(y) > s0:   
                w_high = w
                inf = 1
                break

        if inf == 0:
            w_low = w

        w = (w_high+w_low)/2                                    
        y_x = odeint(qball_ode, r0, r, (w,)) 
    
    print(w)

    return exp_fit(w, r, y_x)

def qball_ode(y, r, w):
    s = y[0]
    u = y[1]

    if r == 0:
        dydr = [u, (1/2)*((w**2 - 1)*s0-s0**2+B*s0**3)]
    else:
        dydr = [u, -(1/r)*u-s*(w**2 - 1)-A*s**2 + B*s**3]   

    return dydr

def exp_fit(w, r, y_x):
    start = 0
    diff = abs(y_x[:,0][0]-y_x[:,0][1])
    fit_threshold =  1e-5

    for i in range(len(y_x[:,0])):
        new_diff = abs(y_x[:,0][i]-y_x[:,0][i+1])
        if abs(y_x[:,0][i])<=fit_threshold and not start:
            start = i
        if start and new_diff>diff:
            end = i-1                                       
            break
        diff = new_diff

    x1 = r[start]
    x2 = r[end]
    y1 = y_x[:,0][start]
    y2 = y_x[:,0][end]

    b = (np.log(y2)-np.log(y1))/(x1-x2)
    a = y1/(np.exp(-b*x1))

    popt = [a,b]

    # making and saving the piecewise function for the s0
    piecewise_y = np.append(y_x[:,0][:start], exp_func(r[start:], *popt))

    return [w, np.transpose([r, piecewise_y])]

def exp_func(x, a, b):
    return a * np.exp(-b * x)


params_list = np.loadtxt('params.csv', delimiter=' ', skiprows=1)
if params_list.ndim == 1:
    params_list = np.array([params_list])
A, B, domain, n = np.genfromtxt('params.csv', delimiter=' ', skip_footer=len(params_list))
domain = int(domain)
n = int(n)

sigma0_data = {}
phi01_system = np.array([])
phi02_system = np.array([])
pi01_system = np.array([])
pi02_system = np.array([])

for params in params_list:
    s0, high_bound, low_bound, pos_x, pos_y, v_x, v_y = params
    pos = [pos_x, pos_y]
    v = [v_x, v_y]

    if s0 not in sigma0_data:
        bounds = [high_bound, low_bound]
        sigma0_data[s0] = ode_solver(s0, bounds)

    w, data = sigma0_data[s0]

    phi01_curr_qball, phi02_curr_qball, pi01_curr_qball, pi02_curr_qball = qball_maker(w, pos, v, data)
    if phi01_system.size == 0:               # first qball
        phi01_system = phi01_curr_qball
        phi02_system = phi02_curr_qball
        pi01_system = pi01_curr_qball
        pi02_system = pi02_curr_qball
    else:
        phi01_system = phi01_system + phi01_curr_qball
        phi02_system = phi02_system + phi02_curr_qball
        pi01_system = pi01_system + pi01_curr_qball
        pi02_system = pi02_system + pi02_curr_qball

x, y = np.meshgrid(np.linspace(-domain, domain, n), np.linspace(-domain, domain, n))

# plotting stuff
fig = plt.figure(figsize = (10, 7))
gs = GridSpec(2, 12)

ax1 = fig.add_subplot(gs[0,0:6], projection ="3d")
surf = ax1.plot_surface(x, y, phi01_system, cmap=cm.magma)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel(r'$\phi^0_1$')
ax1.set_title(r'$\phi^0_1=\sigma(\gamma x,y) \cos(\omega\gamma vx)$')

ax2 = fig.add_subplot(gs[0,6:12], projection ="3d")
surf2 = ax2.plot_surface(x, y, phi02_system, cmap=cm.magma)
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel(r'$\phi^0_2$')
ax2.set_title(r'$\phi^0_2=\sigma(\gamma x,y) \sin(\omega\gamma vx)$')

ax3 = fig.add_subplot(gs[1,0:6], projection ="3d")
surf3 = ax3.plot_surface(x, y, pi01_system, cmap=cm.viridis)
ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax3.set_zlabel(r'$\pi^0_1$')
ax3.set_title(r'$\pi^0_1=v\gamma \cos(\omega\gamma vx)\frac{\partial\sigma(\gamma x,y)}{\partial x}-\omega\gamma\phi_2$')

ax4 = fig.add_subplot(gs[1,6:12], projection ="3d")
surf4 = ax4.plot_surface(x, y, pi02_system, cmap=cm.viridis)
ax4.set_xlabel('x')
ax4.set_ylabel('y')
ax4.set_zlabel(r'$\pi^0_2$')
ax4.set_title(r'$\pi^0_2=v\gamma \sin(\omega\gamma vx)\frac{\partial\sigma(\gamma x,y)}{\partial x}+\omega\gamma\phi_1$')

fig.colorbar(surf, ax=ax1, shrink=0.5, aspect=5)
fig.colorbar(surf2, ax=ax2, shrink=0.5, aspect=5)
fig.colorbar(surf3, ax=ax3, shrink=0.5, aspect=5)
fig.colorbar(surf4, ax=ax4, shrink=0.5, aspect=5)

# outputting stuff
xys = [x.flatten(), y.flatten(), phi01_system.flatten(), phi02_system.flatten(), pi01_system.flatten(), pi02_system.flatten()]
np.savetxt('initdata_new.csv', np.transpose(xys), delimiter=' ')

fig.tight_layout()
plt.show()



