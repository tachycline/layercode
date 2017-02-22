import numpy as np
from scipy.integrate import ode
from matplotlib import pyplot as plt
from matplotlib import animation


class flow(object):
    def __init__(self, nx=64, ny=64, lx=1.0e6, ly=2.0e6,
                 beta=6.0e-10, kappa=10.0, w=4.0e-6, forcing_mode=(6,9),
                 method='dopri5', dt=2*3600):
        x = np.linspace(-np.pi, np.pi, nx, endpoint=False)
        y = np.linspace(-np.pi, np.pi, ny, endpoint=False)
        self.xx, self.yy = np.meshgrid(x,y, indexing='ij')
        
        kx = np.fft.fftfreq(nx, lx/nx)
        ky = np.fft.rfftfreq(ny, ly/ny)
        self.kkx, self.kky = np.meshgrid(kx,ky, indexing='ij')
        
        self.ksq = self.kkx**2 + self.kky**2
        self.ksq[0,0] += 1.0e-15
        self.kmax = np.min(np.max(kx), np.max(ky))
        self.kxmax = np.max(kx)
        self.kymax = np.max(ky)
        
        self.nx = nx
        self.ny = ny
        self.lx = lx
        self.ly = ly
        self.beta = beta
        self.kappa = kappa
        self.w = w
        
        self.forcing_mode = forcing_mode

        self.psihat = np.zeros_like(self.ksq, dtype=complex)

        # should we do random noise instead for initial condition?
#        self.psihat[2,4] = nx*ny*4
#        self.psihat[4,2] = nx*ny
#        self.psi = np.fft.irfft2(self.psihat)
        
#        self.qhat = -self.psihat*self.ksq
#        self.q = np.fft.irfft2(self.qhat)

        self.qhat = self.forcing(0.0)
        self.q = np.fft.irfft2(self.qhat)
        self.psihat = self.get_psihat_from_qhat(self.qhat)
        self.psi = np.fft.irfft2(self.psihat)


        self.integrator = ode(self.rhs).set_integrator(method)
        self.results = [ self.qhat ]
        self.dt = dt
        self.t = 0.0

        
    def integrate(self, tf):
        
        if tf < self.t + self.dt:
            print("won't integrate backward in time from %f to %f" %(self.t,tf))
            return
        
        y0 = self.munge(self.qhat)
        self.integrator.set_initial_value(y0, self.t)
        while self.integrator.successful() and self.integrator.t < tf:
            self.integrator.integrate(self.integrator.t + self.dt)
            self.results.append(self.unmunge(self.integrator.y))
            
        self.qhat = self.unmunge(self.integrator.y)
        self.psihat = self.get_psihat_from_qhat(self.qhat)
        self.psi = np.fft.irfft2(self.psihat)
        self.q = np.fft.irfft2(self.qhat)
        self.t = self.integrator.t
        
        return
            
    def plot_psi(self):
        return plt.contour(self.xx, self.yy, self.psi)
        
    def plot_q(self):
        return plt.contour(self.xx, self.yy, self.q)
        
        
    def get_psihat_from_qhat(self, qhat):
        """What it says on the tin.        
        """
        
        psihat = qhat/self.ksq
        return psihat
        
    def waveterm(self, psihat):
        """Compute the beta wave term.
        
        Assume that we start and end in Fourier space.
        """
        return self.beta*psihat*self.kkx*(0.0+1.0j)
    
    def dissipation(self, qhat):
        """Dissipation term, all in Fourier space."""
        
        return self.kappa*qhat*self.ksq
    
    def forcing(self, t):
        """Forcing term goes here.
        
        This ought to be random phases into a k-space anulus, but I'll use something simple for now.
        """
        fzero = 1.0e-4
        thick = 1.0e3
        famp = self.w * fzero/thick

        phases = np.random.uniform(-np.pi, np.pi, size=self.ksq.shape)
        
        forcing = np.zeros_like(self.ksq, dtype=complex)
        forcing[np.abs(self.ksq - self.kmax**2/9) < 3e-11] = famp
        forcing *= (np.cos(phases) + np.sin(phases)*(0.0+1.0j))
        
        return forcing
    
    def nlterm(self, qhat, psihat):
        """Compute the jacobian determinant."""
        
        psihat_x = psihat*self.kkx*(0.0+1.0j)
        psihat_y = psihat*self.kky*(0.0+1.0j)
        qhat_x = qhat*self.kkx*(0.0+1.0j)
        qhat_y = qhat*self.kky*(0.0+1.0j)
        
        psi_x = np.fft.irfft2(psihat_x)
        psi_y = np.fft.irfft2(psihat_y)
        q_x = np.fft.irfft2(qhat_x)
        q_y = np.fft.irfft2(qhat_y)
        
        jac = psi_x*q_y - psi_y*q_x
        
        jachat = np.fft.rfft2(jac)
        
        # dealias
        jachat[self.kkx > 2/3 * self.kxmax] = 0.0
        jachat[self.kky > 2/3 * self.kymax] = 0.0
        
        return jachat
    
    def rhs(self, arg1, arg2):
        """The time derivative, ready for the integrator."""
        
        # scipy.ode and scipy.odeint use opposite call signatures
        # so we have to figure out which of the arguments is a float
        # and which is an array
        
        if type(arg1) == type(0.0):
            t = arg1
            q_reshaped = arg2
        else:
            t = arg2
            q_reshaped = arg1
            
        qhat = self.unmunge(q_reshaped)
        
        psihat = self.get_psihat_from_qhat(qhat)
        nlterm = self.nlterm(qhat, psihat)
        waveterm = self.waveterm(psihat)
        dissipation = self.dissipation(qhat)
        forcing = self.forcing(t)
        
        return self.munge(forcing - dissipation + waveterm + nlterm)
    
    def calc_energy(self, psihat):
        uhat = psihat*self.kkx*(0.0+1.0j)
        vhat = psihat*self.kky*(0.0+1.0j)
        
        u = np.fft.irfft2(uhat)
        v = np.fft.irfft2(vhat)
        
        efield = u**2 + v**2
        return np.sum(efield)
    
    def calc_enstrophy(self, qhat):
        q = np.fft.irfft2(qhat)
        return np.sum(q**2)
    
    def munge(self, qhat):
        """format a complex k-space field for odeint"""
        
        r = qhat.real
        i = qhat.imag
        z = np.array([r,i])
        return z.reshape(-1)
    
    def unmunge(self, munged):
        """Return the 1d real sequence to its 2d complex state"""
        
        z = munged.reshape((2,self.nx,int(self.ny/2+1)))
        r = z[0]
        i = z[1]
        return r + (0+1.0j)*i
    
    
    def animate_results(self, filename, stride=1):
        fig = plt.figure(figsize=(10,10))
        ax = plt.axes()

        plt.xlabel(r'x')
        plt.ylabel(r'y')

        def animate(i):
            qhat = self.results[int(i*stride)]
            psihat = self.get_psihat_from_qhat(qhat)
            z = np.fft.irfft2(psihat)
            ax.clear()
            cont = plt.contour(self.xx, self.yy, z)
            return cont

        anim = animation.FuncAnimation(fig, animate, frames=len(self.results)//stride, blit=False)
        mywriter = animation.FFMpegWriter()
        anim.save(filename, bitrate=10000)
