import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['mathtext.fontset'] = 'stix'  # Use STIX fonts for math symbols
plt.rcParams['text.usetex'] = True # enable LaTeX rendering

def ypoints(a, b, Ny):
    ypts = np.zeros(Ny)
    c = 0.5 * (b + a)
    r = 0.5 * (b - a)
    piN = np.pi / (Ny - 1)
    for i in range(Ny):
        ypts[i] = c + r * np.cos(piN * i)
    return ypts

def smooth(x,y,z,scale):
    from scipy.interpolate import griddata
    X, Y = np.meshgrid(x, y)
    x_fine = np.linspace(x[0], x[-1], scale*len(x))
    y_fine = np.linspace(y[0], y[-1], scale*len(y))
    X_fine, Y_fine = np.meshgrid(x_fine, y_fine)
    Z_fine = griddata((X.flatten(), Y.flatten()), z.flatten(), (X_fine, Y_fine), method='cubic')
    return x_fine,y_fine,Z_fine

def scale(x,y,z,Nx=-1,Ny=-1):
    from scipy.interpolate import griddata
    X, Y = np.meshgrid(x, y)
    x_fine = x
    if Nx != len(x) and Nx > 0:
        x_fine = np.linspace(x[0], x[-1], Nx)
    y_fine = y
    if Ny != len(y) and Ny > 0:
        y_fine = np.linspace(y[0], y[-1], Ny)
    X_fine, Y_fine = np.meshgrid(x_fine, y_fine)
    Z_fine = griddata((X.flatten(), Y.flatten()), z.flatten(), (X_fine, Y_fine), method='cubic')
    return x_fine,y_fine,Z_fine

def x_averaged_scalar(file_path,scale=1):
    import netCDF4 as nc
    ds = nc.Dataset(file_path)
    scalar = ds.variables['Component_0'][:]
    scalar=np.concatenate((scalar, scalar[0,:, :][np.newaxis,:,:]), axis=0) # merge [z,y,x] and [0,y,x]
    y = ds.variables['Y'][:]
    z = ds.variables['Z'][:]
    z = np.concatenate((z, (z[-1]+z[1]-z[0])[np.newaxis]), axis=0)
    averaged_scalar = np.mean(scalar, axis=-1)
    return smooth(y,z,averaged_scalar,scale)



def x_averaged_density(file_path_temp,file_path_salt,Rrho=2,scale=1):
    [y,z,averaged_temp] = x_averaged_scalar(file_path_temp,scale)
    [y,z,averaged_salt] = x_averaged_scalar(file_path_salt,scale)
    averaged_density = (Rrho*averaged_salt-averaged_temp)/(Rrho-1)
    return y,z,averaged_density

def xz_averaged_salar(file_path):
    import netCDF4 as nc
    df = nc.Dataset(file_path)
    salar = df.variables['Component_0'][:]
    averaged_scalar = np.mean(salar, axis=(0,2))
    y = df.variables['Y'][:]
    return y,averaged_scalar

def xz_averaged_density(file_path_temp,file_path_salt,Rrho):
    [y,averaged_temp] = xz_averaged_salar(file_path_temp)
    [y,averaged_salt] = xz_averaged_salar(file_path_salt)
    averaged_density = (Rrho*averaged_salt-averaged_temp)/(Rrho-1)
    return y,averaged_density

def readdata(path,index):
    import pandas as pd
    df = pd.read_csv(path, sep='\s+')
    array = df.values
    mu = array[:,0]
    vals = array[:,index]
    return mu,vals

def readdata_po(path,index):
    import pandas as pd
    minvals = []
    maxvals = []
    [mu_values, directories] = parse_file(path+'MuD.asc') # get dir of each search
    for mu, directory in zip(mu_values, directories):
        try:
            df = pd.read_csv(path+directory+'energy.asc', sep='\s+')
            array = df.values
            vals = array[:,index]
            minvals.append(np.min(vals))
            maxvals.append(np.max(vals))
        except FileNotFoundError:
            print("Error: The file was not found.")
            minvals.append(minvals[-1])
            maxvals.append(maxvals[-1])
    return mu_values,minvals,maxvals

def x_averaged_vec(file_path,scale=1):
    import netCDF4 as nc
    ds = nc.Dataset(file_path)
    vec_x = ds.variables['Velocity_X'][:]
    vec_y = ds.variables['Velocity_Y'][:]
    vec_z = ds.variables['Velocity_Z'][:]
    y = ds.variables['Y'][:]
    z = ds.variables['Z'][:]
    vec_x=np.concatenate((vec_x, vec_x[0,:, :][np.newaxis,:,:]), axis=0)
    vec_y=np.concatenate((vec_y, vec_y[0,:, :][np.newaxis,:,:]), axis=0)
    vec_z=np.concatenate((vec_z, vec_z[0,:, :][np.newaxis,:,:]), axis=0)
    z = np.concatenate((z, (z[-1]+z[1]-z[0])[np.newaxis]), axis=0)
    averaged_vec_x = np.mean(vec_x, axis=-1)
    averaged_vec_y = np.mean(vec_y, axis=-1)
    averaged_vec_z = np.mean(vec_z, axis=-1)
    yy,zz,averaged_vec_x=smooth(y,z,averaged_vec_x,scale)
    yy,zz,averaged_vec_y=smooth(y,z,averaged_vec_y,scale)
    yy,zz,averaged_vec_z=smooth(y,z,averaged_vec_z,scale)
    return yy,zz,averaged_vec_x,averaged_vec_y,averaged_vec_z

def z_averaged_scalar(file_path,scale=1):
    import netCDF4 as nc
    ds = nc.Dataset(file_path)
    scalar = ds.variables['Component_0'][:] # [z,y,x]
    scalar=np.concatenate((scalar, scalar[:,:,0][:,:,np.newaxis]), axis=2) # merge [z,y,x] and [z,y,x0]
    y = ds.variables['Y'][:]
    x = ds.variables['X'][:]
    x = np.concatenate((x, (x[-1]+x[1]-x[0])[np.newaxis]), axis=0)
    averaged_scalar = np.mean(scalar, axis=0)
    return smooth(y,x,averaged_scalar.T,scale)

def z_averaged_vec(file_path,scale=1):
    import netCDF4 as nc
    ds = nc.Dataset(file_path)
    vec_x = ds.variables['Velocity_X'][:] # [z,y,x]
    vec_y = ds.variables['Velocity_Y'][:]
    vec_z = ds.variables['Velocity_Z'][:]
    y = ds.variables['Y'][:]
    x = ds.variables['X'][:]
    vec_x=np.concatenate((vec_x, vec_x[:,:,0][:,:,np.newaxis]), axis=2)
    vec_y=np.concatenate((vec_y, vec_y[:,:,0][:,:,np.newaxis]), axis=2)
    vec_z=np.concatenate((vec_z, vec_z[:,:,0][:,:,np.newaxis]), axis=2)
    x = np.concatenate((x, (x[-1]+x[1]-x[0])[np.newaxis]), axis=0)
    averaged_vec_x = np.mean(vec_x, axis=0)
    averaged_vec_y = np.mean(vec_y, axis=0)
    averaged_vec_z = np.mean(vec_z, axis=0)
    yy,xx,averaged_vec_x=smooth(y,x,averaged_vec_x.T,scale)
    yy,xx,averaged_vec_y=smooth(y,x,averaged_vec_y.T,scale)
    yy,xx,averaged_vec_z=smooth(y,x,averaged_vec_z.T,scale)
    return yy,xx,averaged_vec_x,averaged_vec_y,averaged_vec_z


def xz_averaged_vec(file_path,scale=1):
    import netCDF4 as nc
    ds = nc.Dataset(file_path)
    vec_x = ds.variables['Velocity_X'][:] # [z,y,x]
    vec_y = ds.variables['Velocity_Y'][:]
    vec_z = ds.variables['Velocity_Z'][:]
    y = ds.variables['Y'][:]
    averaged_vec_x = np.mean(vec_x, axis=(0,2))
    averaged_vec_y = np.mean(vec_y, axis=(0,2))
    averaged_vec_z = np.mean(vec_z, axis=(0,2))
    return y,averaged_vec_x,averaged_vec_y,averaged_vec_z

def plot_scalar(x,y,scalar,vmin=-1,vmax=1,colorbar=False):
    X, Y = np.meshgrid(x, y)
    cax = plt.pcolormesh(X,Y,scalar.T, shading='gouraud', cmap='bwr',vmin=vmin,vmax=vmax, zorder=1)
    plt.xlabel(r'$x$', fontsize=16)
    plt.ylabel(r'$z$', fontsize=16)
    plt.xticks([0,np.pi,2*np.pi], [r'0', r'$\pi$', r'$2\pi$'],fontsize=16)
    plt.yticks([-0.5, 0, 0.5],[r'$-0.5$',r'0',r'0.5'],fontsize=16)
    plt.xlim(0, 2*np.pi)
    plt.ylim(-0.5, 0.5)
    if colorbar:
        ticks=[vmin, 0, vmax]
        cbar = plt.colorbar(cax, ticks=ticks, shrink=1.0, aspect=10)
        # cbar.ax.set_aspect(1)
        cbar.ax.set_yticklabels([r'{:.2f}'.format(vmin), r'0', r'{:.2f}'.format(vmax)])

def plot_streamline(x,y,vec_x,vec_y):
    X, Y= np.meshgrid(x,y)
    sort_indices = np.argsort(Y[:,0]) 
    Y_sorted = Y[sort_indices, :] 
    X_sorted = X[sort_indices, :] 
    vec_x_sorted = vec_x[:, sort_indices] 
    vec_y_sorted = vec_y[:, sort_indices] 
    speed = np.sqrt(vec_x_sorted**2 + vec_y_sorted**2)
    # lw = 1.2 * speed / speed.max()
    lw = 1.0 * speed / speed.max()
    # contour_levels = [0.1, 0.2, 0.5, 1.0]*speed.max()
    plt.streamplot(X_sorted, Y_sorted, vec_x_sorted.T, vec_y_sorted.T, 
                   density = 0.7, color ='k', linewidth = lw.T, arrowstyle ='->',arrowsize = 0.8, broken_streamlines=True) 
    # plt.streamplot(X_sorted, Y_sorted, vec_x_sorted.T, vec_y_sorted.T, 
    #                density = 1, linewidth = lw.T,color ='k', broken_streamlines=False) 
    # contour = plt.contour(X, Y, vec_x_sorted.T, levels = contour_levels, colors='k',linewidths = 0.6)

def plot_contour(x,y,vec_x):
    # X, Y= np.meshgrid(x,y)
    # sort_indices = np.argsort(Y[:,0]) 
    # vec_x_sorted = vec_x[:, sort_indices] 
    # contour_levels = np.linspace(vec_x_sorted.min(), vec_x_sorted.max(), 10)
    # contour = plt.contour(X, Y, vec_x_sorted.T, levels = contour_levels, colors='k',linewidths = 0.6)
    X, Y= np.meshgrid(x,y)
    # sort_indices = np.argsort(Y[:,0]) 
    # vec_x_sorted = vec_x[:, sort_indices] 
    contour_levels = np.linspace(vec_x.min(), vec_x.max(), 10)
    contour = plt.contour(X, Y, vec_x.T, levels = contour_levels, colors='k',linewidths = 0.6)


def plot_contour_vorticity(x,y,vec_x,vec_y):
    X, Y= np.meshgrid(x,y)
    sort_indices = np.argsort(Y[:,0]) 
    Y_sorted = Y[sort_indices, :] 
    X_sorted = X[sort_indices, :] 
    vec_x_sorted = vec_x[:, sort_indices] 
    vec_y_sorted = vec_y[:, sort_indices] 
    dVdx = np.gradient(vec_y_sorted, axis=1)  # ∂v/∂x
    dUdy = np.gradient(vec_x_sorted, axis=0)  # ∂u/∂y
    vorticity = dVdx - dUdy  # ω = ∂v/∂x - ∂u/∂y
    contour_levels = np.linspace(vorticity.min(), vorticity.max(), 10)
    contour = plt.contour(X, Y, vorticity.T, levels = contour_levels, colors='k',linewidths = 0.6)

def plot_contour_streamfunc(x,y,vec_x,vec_y):
    X, Y= np.meshgrid(x,y)
    sort_indices = np.argsort(Y[:,0]) 
    Y_sorted = Y[sort_indices, :] 
    X_sorted = X[sort_indices, :] 
    vec_x_sorted = vec_x[:, sort_indices] 
    vec_y_sorted = vec_y[:, sort_indices] 
    from scipy import integrate
    int_x = integrate.cumtrapz(vec_y_sorted, X, axis=1, initial=0)[0]  # Integrate v along x
    int_y = integrate.cumtrapz(vec_x_sorted, Y, axis=0, initial=0)      # Integrate u along y
    streamfunction = -int_x + int_y
    contour_levels = np.linspace(streamfunction.min(), streamfunction.max(), 10)
    contour = plt.contour(X, Y, streamfunction.T, levels = contour_levels, colors='k',linewidths = 0.6)

def standard_deviation(filepath):
    [y,T] = xz_averaged_salar(filepath)
    return y,abs(T)
def conductive_flux(filepath,smooth=1,delta=1,H=1):
    import scipy.interpolate as interp
    [y,T] = xz_averaged_salar(filepath)
    interp_func = interp.BarycentricInterpolator(y, T)
    y_smooth = np.linspace(min(y), max(y), len(y)*smooth)
    T_smooth = interp_func(y_smooth)
    print(np.shape(y_smooth),np.shape(T_smooth))
    return y_smooth,-np.gradient(T_smooth-y_smooth, y_smooth)/(delta/H)
def convective_flux(tpath,upath,kappa,smooth=1,delta=1,H=1):
    [y,_,T] = z_averaged_scalar(tpath)
    [_,_,_,w,_] = z_averaged_vec(upath)
    wT = w*T
    print(np.shape(wT))
    wT_m = np.mean(wT, axis=(0))
    return y,wT_m/(kappa*delta/H)
def standard_deviation(filepath):
    [y,T] = xz_averaged_salar(filepath)
    return y,abs(T)
def thickness(path):
    y,dT = standard_deviation(path)
    from scipy.signal import find_peaks
    peaks, _ = find_peaks(dT[:])
    thickness = abs(y[0]-y[peaks[0]])
    # plt.plot(y[peaks],dT[peaks], "ro", label="Peaks")
    # plt.plot(y,dT)
    return thickness
def thickness_ratio_eq(tpath,spath):
    ht = thickness(tpath)
    hc = thickness(spath)
    return ht/hc
def flux_ratio(tpath,spath,Lambda=2,tau=0.01):
    Nu = Nusselt(tpath)
    Sh = Nusselt(spath)
    return Sh/Nu*Lambda*tau

def flux_ratio_po(tpath,spath,upath,kappaT,kappaS,Lambda=2,tau=0.01):
    [y,Fv_t]=convective_flux(tpath,upath,kappaT)
    [y,Fd_t]=conductive_flux(tpath)
    Nu = Fv_t+Fd_t
    # plt.plot(y,Nu)
    [y,Fv_c]=convective_flux(spath,upath,kappaS)
    [y,Fd_c]=conductive_flux(spath)
    Sh = Fv_c+Fd_c
    # plt.plot(y,Sh)
    return y,Sh/Nu*Lambda*tau

def Nusselt(path):
    [y,scalar_h] = xz_averaged_salar(path)
    totscalar_h = scalar_h-y
    nusselt = abs((totscalar_h[1]-totscalar_h[0])/(y[1]-y[0]))
    return nusselt
def get_Nu(dir):
    [Ra, Dirs] = parse_file(dir+"MuD.asc")
    Ra = np.array(Ra, dtype=float)
    Nu = np.zeros(len(Ra))
    for i in range(len(Dirs)):
        Nu[i]=Nusselt(dir+Dirs[i]+'tbest.nc')
    return Ra,Nu
def get_Sh(dir):
    [Ra, Dirs] = parse_file(dir+"MuD.asc")
    Ra = np.array(Ra, dtype=float)
    Sh = np.zeros(len(Ra))
    for i in range(len(Dirs)):
        Sh[i]=Nusselt(dir+Dirs[i]+'sbest.nc')
    return Ra,Sh

def get_bifurcation_data_po(foldername, indexE, stableResidual=1e-6):
    mu, dirs = parse_file(foldername+"/MuD.asc")
    isunstable = []
    for index in range(len(mu)):
        # print(foldername+"/"+dirs[index]+"Residu.asc")
        real_list, imag_list = load_eigenvals(foldername+"/"+dirs[index]+"lambda.asc",foldername+"/"+dirs[index]+"Residu.asc",stableResidual)
        if len(real_list) > 0:
            max_idx = np.argmax(real_list)
            max_real = real_list[max_idx]
            print(dirs[index],mu[index],max_real,imag_list[max_idx])
            if max_real>stableResidual:
                isunstable.append(True)
            else:
                isunstable.append(False)
        else:
            # print("The array is empty.")
            isunstable.append(True)
    return mu, isunstable

def get_bifurcation_data(foldername, indexE, stableResidual=1e-6):
    mu, dirs = parse_file(foldername+"/MuD.asc")
    mu, refes=readdata(foldername+"/MuE.asc",indexE)
    if len(dirs) != len(refes):
        print("Sizes of MuE = "+str(len(refes))+" and MuD = "+str(len(dirs))+" are not consistent!")
        return
    isunstable = []
    for index in range(len(mu)):
        # print(foldername+"/"+dirs[index]+"Residu.asc")
        real_list, imag_list = load_eigenvals(foldername+"/"+dirs[index]+"lambda.asc",foldername+"/"+dirs[index]+"Residu.asc",stableResidual)
        if len(real_list) > 0:
            # max_real = max(real_list)
            # print(mu[index],max_real)
            max_idx = np.argmax(real_list)
            # print(max_idx)
            # print(real_list)
            max_real = real_list[max_idx]
            print(dirs[index],mu[index],max_real,imag_list[max_idx])
            if max_real>0:
                isunstable.append(True)
            else:
                isunstable.append(False)
        else:
            # print("The array is empty.")
            isunstable.append(True)
    return mu, refes, isunstable
    
def separate_into_segments(mu,refes,isunstable):
    # Initialize lists for segments
    segments = []
    current_segment = {'mu': [], 'refes': [], 'isunstable': isunstable[0]}
    # Separate the data into segments based on the condition
    for i in range(len(mu)):
        if isunstable[i] == current_segment['isunstable']:
            current_segment['mu'].append(mu[i])
            current_segment['refes'].append(refes[i])
        else:
            segments.append(current_segment)
            current_segment = {'mu': [mu[i-1], mu[i]], 'refes': [refes[i-1], refes[i]], 'isunstable': isunstable[i]}
    segments.append(current_segment)
    return segments

def psd(t,x0,name):
    # Parameters
    distance = 10.0 # time period for computing (from time start of signal)
    delta_t = t[1] - t[0] # time step
    pos_start = int(distance // delta_t) # integer position at time to start computing PSD

    # modify signal and reset initial point
    t = t[pos_start:] - t[pos_start]
    x0 = x0[pos_start:]
    # L = int(np.floor(distance / delta_t)) + 1
    L = len(x0) + 1

    # define frequency vector
    f_low  = 1./t[len(t)-1] # lowest frequency for this range
    fs = 1./delta_t # sampling frequency
    f = np.arange(1, L) * f_low

    # Calculate FFT
    FFT0 = np.fft.fft(x0)

    # double sided spectrum
    P1_0 = np.abs(FFT0) / L

    # Single-sided spectrum
    P2_0 = P1_0[:L//2]
    P2_0[1:len(P2_0)-2] = 2 * P2_0[1:len(P2_0)-2]
    f2_0 = f[:L//2]

    # Find maximum PSD value and corresponding frequency
    # max_PSD = np.max(P2_0)
    # max_freq_index = np.argmax(P2_0)
    # max_freq = f2_0[max_freq_index]

    # # Normalize max PSD
    # k = 10000000 / max_PSD
    # while abs(k) > 100:
    #     k /= 10
    # k = abs(k)
    # max_PSD *= k

    # # Find minimum PSD value and corresponding frequency
    # min_PSD = np.min(P2_0)
    # min_freq_index = np.argmin(P2_0)
    # min_freq = f2_0[min_freq_index]

    # # Normalize min PSD
    # h = 100 / min_PSD
    # while abs(h) > 1:
    #     h /= 10
    # min_PSD *= h

    # # Find max and min frequencies
    # max_x = np.max(f2_0)
    # min_x = np.min(f2_0)

    # # Print results (optional)
    # print(f"Max PSD: {max_PSD:.6f}")
    # print(f"Min PSD: {min_PSD:.6f} at frequency {min_freq:.2f} Hz")
    # print(f"Max frequency: {max_x:.2f} Hz")
    # print(f"Min frequency: {min_x:.2f} Hz")

    # Plot the PSD figure
    plt.figure(figsize=(4, 3))
    # plt.gca().set_aspect(ASPECT*0.5)
    plt.loglog(f2_0[2:len(f2_0)-2], P2_0[2:len(P2_0)-2], label="E1")
    plt.xlabel("St")
    plt.ylabel("PSD")
    plt.grid(True, which="both", ls="-")
    # plt.legend()

    # save figure
    plt.savefig(name+"_PSD_St.png", dpi=300, bbox_inches="tight")
    plt.show()

def psd2(t, x0, name, L_ref=1.0, U_ref=1.0):
    # Parameters
    distance = 10.0  # time period for computing
    delta_t = t[1] - t[0]
    pos_start = int(distance // delta_t)

    # Modify signal and reset time
    t = t[pos_start:] - t[pos_start]
    x0 = x0[pos_start:]
    L = len(x0)

    # Frequency vector
    fs = 1.0 / delta_t  # sampling frequency
    f = np.fft.fftfreq(L, d=delta_t)
    f2_0 = f[:L//2]

    # FFT
    FFT0 = np.fft.fft(x0)
    P1_0 = (np.abs(FFT0) ** 2) / L  # power spectrum
    P2_0 = P1_0[:L//2]
    P2_0[1:-1] = 2 * P2_0[1:-1]  # single-sided

    # Convert to Strouhal number
    St = f2_0 * L_ref / U_ref

    # Find maximum PSD and corresponding St
    max_idx = np.argmax(P2_0)
    max_PSD = P2_0[max_idx]
    max_St = St[max_idx]

    # print(f"Max PSD = {max_PSD:.6e} at St = {max_St:.6f}")

    # Plot PSD
    plt.figure(figsize=(4, 3))
    plt.loglog(St, P2_0, label="PSD")
    plt.xlabel("Frequency", fontsize=14)
    plt.ylabel("PSD", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim(1e-3, 0.5)
    plt.ylim(bottom=1e-4, top=1e4)
    plt.grid(True, which="both", ls="-")

def mergeImagesToVideo(path,videoname,frames=10):
    import os
    comd = "ffmpeg -r "+str(frames)+" -i "+path+" -c:v libx264 -pix_fmt yuv420p -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" "+videoname
    print(comd)
    os.system(comd)

def getpath(path):
    import glob
    from natsort import natsorted
    file_paths = glob.glob(path)
    return natsorted(file_paths)

def loaddata_channelflow_ddc(path,dt=1):
    import numpy as np
    import glob
    import pandas as pd
    file_list = getpath(path)
    nstep = len(file_list)
    
    t = np.linspace(0, nstep, nstep) * dt

    # load data
    T_mean = []
    S_mean = []
    T_gradient = []
    S_gradient = []
    T_flux = []
    S_flux = []
    y = []
    for file in file_list:
        df = pd.read_csv(file, skiprows=1)
        array = df.values
        y = array[:,0]
        
        T_mean.append(array[:,1])
        S_mean.append(array[:,2])
        T_gradient.append(array[:,3])
        S_gradient.append(array[:,4])
        T_flux.append(array[:,5])
        S_flux.append(array[:,6])
    T_mean = np.array(T_mean)
    S_mean = np.array(S_mean)
    T_gradient = np.array(T_gradient)
    S_gradient = np.array(S_gradient)
    T_flux = np.array(T_flux)
    S_flux = np.array(S_flux)

    return t,y,T_mean,S_mean,T_gradient,S_gradient,T_flux,S_flux

def plot_meandata(t,y,tm,sm,tgrad,sgrad,tflux,sflux,Rrho,startT,endT,name='data'):
    import matplotlib.pyplot as plt
    X, Y = np.meshgrid(t, y)
    fig = plt.figure(figsize=(14, 4))

    ax1 = plt.subplot(2, 3, 1)
    plt.pcolormesh(X,Y,tm.T, cmap='bwr',shading='gouraud')
    plt.ylabel(r'$z$',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
    plt.title(r'$\langle T\rangle_h$')

    ax2 = plt.subplot(2, 3, 2)
    plt.pcolormesh(X,Y,tgrad.T, cmap='bwr',shading='gouraud')
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=False)
    plt.title(r'$\langle\partial_y T\rangle_h$')

    ax3 = plt.subplot(2, 3, 3)
    plt.pcolormesh(X,Y,tflux.T, cmap='bwr',shading='gouraud')
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=False)
    plt.title(r'$\langle vT\rangle_h$')
    
    ax4 = plt.subplot(2, 3, 4)
    plt.pcolormesh(X,Y,sm.T, cmap='bwr',shading='gouraud')
    plt.xlabel(r'$t$',fontsize=14)
    plt.ylabel(r'$z$',fontsize=14)
    plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.title(r'$\langle S\rangle_h$')
    
    ax5 = plt.subplot(2, 3, 5, adjustable='box')
    plt.pcolormesh(X,Y,sgrad.T, cmap='bwr',shading='gouraud')
    plt.xlabel(r'$t$',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=False)
    plt.title(r'$\langle \partial_y S\rangle_h$')
    
    ax6 = plt.subplot(2, 3, 6)
    plt.pcolormesh(X,Y,sflux.T, cmap='bwr',shading='gouraud')
    plt.xlabel(r'$t$',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=False)
    plt.title(r'$\langle vS\rangle_h$')

    plt.savefig(name+'.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

    mean_T = np.mean(tm[startT:endT,:], axis=-2)
    mean_S = np.mean(sm[startT:endT,:], axis=-2)
    mean_Dens = (Rrho*mean_S-mean_T)/(Rrho-1)
    mean_T_grad = np.mean(tgrad[startT:endT,:], axis=-2)
    mean_S_grad = np.mean(sgrad[startT:endT,:], axis=-2)
    mean_T_flux = np.mean(tflux[startT:endT,:], axis=-2)
    mean_S_flux = np.mean(sflux[startT:endT,:], axis=-2)

    data = np.array([y, mean_T, mean_S, mean_Dens, mean_T_grad, mean_S_grad, mean_T_flux, mean_S_flux])
    np.savetxt(name+'_meanProfiles.csv', data.T, delimiter=', ')

    plt.figure(figsize=(10, 2))
    ax1 = plt.subplot(1, 7, 1)
    plt.plot(mean_T,y)
    plt.ylabel(r'$z$',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
    plt.title(r'$\overline{T}$')
    ax2 = plt.subplot(1, 7, 2)
    plt.plot(mean_S,y)
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=False)
    plt.title(r'$\overline{S}$')
    ax3 = plt.subplot(1, 7, 3)
    plt.plot(mean_Dens,y)
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=False)
    plt.title(r'$\overline{\rho}$')
    ax4 = plt.subplot(1, 7, 4)
    plt.plot(mean_T_grad,y)
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=False)
    plt.title(r'$\overline{\partial_y T}$')
    ax5 = plt.subplot(1, 7, 5)
    plt.plot(mean_S_grad,y)
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=False)
    plt.title(r'$\overline{\partial_y S}$')
    ax6 = plt.subplot(1, 7, 6)
    plt.plot(mean_T_flux,y)
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=False)
    plt.title(r'$\overline{v T}$')
    ax7 = plt.subplot(1, 7, 7)
    plt.plot(mean_S_flux,y)
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=False)
    plt.title(r'$\overline{v S}$')

    plt.savefig(name+'_meanProfiles.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

def plot_meanprofiles(t,y,tm,sm,startT,endT,name='data',Rrho=2):
    import matplotlib.pyplot as plt

    mean_T = np.mean(tm[startT:endT,:], axis=-2)
    mean_S = np.mean(sm[startT:endT,:], axis=-2)
    mean_Dens = (Rrho*mean_S-mean_T)/(Rrho-1)

    plt.figure(figsize=(3, 3))
    plt.plot(mean_T,y,label=r'$\overline{T}$')
    plt.plot(mean_S,y,label=r'$\overline{S}$')
    plt.plot(mean_Dens,y,label=r'$\overline{\rho}$')
    plt.xlabel(r'$\overline{T},\overline{S},\overline{\rho}$',fontsize=14)
    plt.ylabel(r'$z$',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks([-1, 0, 1],['-1','0','1'],fontsize=14)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
    # plt.legend()

    plt.savefig(name+'_meanProfiles.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

def load_eigenvals(lambda_path,residual_path,residual=1e-6):
    import pandas as pd
    try:
        lambda_list = pd.read_csv(lambda_path,skiprows=0,sep=' ') 
        residual_list = pd.read_csv(residual_path,skiprows=0,sep=' ')  
        merged_array = pd.concat([lambda_list, residual_list], axis=1)
        condition = merged_array.iloc[:,2] < residual
        filtered_eigenvals = merged_array[condition]
        return np.array(filtered_eigenvals.iloc[:,0]),np.array(filtered_eigenvals.iloc[:,1])
    except FileNotFoundError:
        print("Error: The file was not found. " + lambda_path + " and "+ residual_path)
        return [],[]
    

def check_instability(real,imag):
    real = np.array(real)
    imag = np.array(imag)
    stable_index = real < 1e-6
    unstable_index = real > 0
    unstable_saddle_index = (real>0) & (imag==0)
    unstable_spiral_index = (real>0) & (imag!=0)
    return stable_index,unstable_index,unstable_saddle_index,unstable_spiral_index

def plot_eigenvals(real,imag,idensity=True):
    if idensity:
        [stable,unstable,saddle,spiral]=check_instability(real,imag)
        plt.scatter(real[stable], imag[stable],color = 'k')
        plt.scatter(real[saddle], imag[saddle],color = 'tab:orange')
        plt.scatter(real[spiral], imag[spiral], marker='s',color = 'tab:green')
    else:
        plt.scatter(real, imag, marker='*', color = 'k')
    plt.xlabel(r'$\mu$',fontsize=14)
    plt.ylabel(r'$\omega$',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
def plot_stability_bifurcation(segments,label=None,color='k',unstableline="dashed"):
    if len(segments)>1:
        # this bifurcation has both stable and unstable bifurcations
        # add label to first stable bifurcation
        addedLabel = False
        for segment in segments:
            if segment['isunstable']:
                plt.plot(segment['mu'], segment['refes'],label=None, linestyle=unstableline, color=color)
            else:
                if addedLabel:
                    plt.plot(segment['mu'], segment['refes'],label=None, linestyle='-', color=color)
                else:
                    plt.plot(segment['mu'], segment['refes'],label=label, linestyle='-', color=color)
                    addedLabel = True
    else:
        # this bifurcation has only one stable or unstable on the whole bifurcation
        # add label to anywhere 
        for segment in segments:
            if segment['isunstable']:
                plt.plot(segment['mu'], segment['refes'],label=label, linestyle=unstableline, color=color)
            else:
                plt.plot(segment['mu'], segment['refes'],label=label, linestyle='-', color=color)
    



def parse_file(filename):
    mu_values = []
    directories = []
    with open(filename, 'r') as file:
        next(file)  # Skip the header row
        for line in file:
            if line.strip():  # Skip empty lines
                columns = line.split()
                mu_values.append(float(columns[0]))  # Extract the first value for Ra
                directories.append(columns[-1].lstrip('./'))  # Extract the last value for directory
    return mu_values, directories

def get_eigen_count(folder):
    mu_values, directories = parse_file(folder+'MuD.asc')
    saddle_count = np.zeros(len(mu_values))
    spiral_count = np.zeros(len(mu_values))
    for index in range(len(directories)):
        [real,imag] = load_eigenvals(folder+directories[index]+'Lambda.asc',folder+directories[index]+'Residu.asc',1e-6)
        [stable,unstable,saddle,spiral]=check_instability(real,imag)
        saddle_count[index] = sum(saddle)
        spiral_count[index] = sum(spiral)
        print(f'mu = {mu_values[index]} in folder \'{directories[index]}\' have {saddle_count[index]} real eigen and {spiral_count[index]} cc eigen')
    return mu_values,saddle_count,spiral_count

def plot_ddc(folder,index,streamline=True):
    for id in index:
        if (id<0):
            filename = folder+'initial-'+str(np.abs(id)-1)+'/'
        else:
            filename = folder+'search-'+str(id)+'/'
        print(filename)
        plt.figure(figsize=(6, 1.5))
        [y,z,scalar] = z_averaged_scalar(filename+'tbest.nc')
        plot_scalar(z,y,scalar,np.min(scalar),np.max(scalar))
        print(np.min(scalar),np.max(scalar))
        if streamline:
            [y,x,vec_x,vec_y,vec_z] = z_averaged_vec(filename+'ubest.nc')
            plot_streamline(x,y,vec_x,vec_y)
        plt.xlabel(r'$x$', fontsize=14)
        plt.ylabel(r'$z$', fontsize=14)
        plt.xticks([0,np.pi,2*np.pi], ['0', r'$\pi$', r'$2\pi$'],fontsize=14)
        plt.yticks([-0.5, 0, 0.5],['-0.5','0','0.5'],fontsize=14)
        plt.xlim(0, 2*np.pi)
        plt.ylim(-0.5, 0.5)
        plt.savefig(filename+'temp.png', dpi=300, bbox_inches='tight')
        plt.show()
        plt.close()
    
        plt.figure(figsize=(6, 1.5))
        [y,z,scalar] = z_averaged_scalar(filename+'sbest.nc')
        plot_scalar(z,y,scalar,np.min(scalar),np.max(scalar))
        print(np.min(scalar),np.max(scalar))
        if streamline:
            [y,x,vec_x,vec_y,vec_z] = z_averaged_vec(filename+'ubest.nc')
            plot_streamline(x,y,vec_x,vec_y)
        plt.xlabel(r'$x$', fontsize=14)
        plt.ylabel(r'$z$', fontsize=14)
        plt.xticks([0,np.pi,2*np.pi], ['0', r'$\pi$', r'$2\pi$'],fontsize=14)
        plt.yticks([-0.5, 0, 0.5],['-0.5','0','0.5'],fontsize=14)
        plt.xlim(0, 2*np.pi)
        plt.ylim(-0.5, 0.5)
        plt.savefig(filename+'salt.png', dpi=300, bbox_inches='tight')
        plt.show()
        plt.close()