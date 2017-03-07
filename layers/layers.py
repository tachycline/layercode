"""Layers of flow."""

import numpy as np
from scipy.integrate import ode
from sqlalchemy import Column, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
import h5py

BASE = declarative_base()

DEFAULTS = {'nx': 64,
            'ny': 64,
            'lx': 1.0e6,
            'ly': 1.0e6,
            'beta': 6.0e-10,
            'kappa': 10.0,
            'w': 4.0e-6,
            'method': 'dopri5',
            'initial_amp': 1.2e5,
            'dt': 2*3600.0,
            'seed': None,
           }

class FlowORM(BASE):
    """The ORM class corresponding to the flow class defined below."""
    __tablename__ = 'flows'

# should probably include a date and a version
    
    id = Column(Integer, primary_key=True)
    nx = Column(Integer)
    ny = Column(Integer)
    lx = Column(Float)
    ly = Column(Float)
    beta = Column(Float)
    kappa = Column(Float)
    w = Column(Float)
    method = Column(String(10))
    dt = Column(Float)
    initial_amp = Column(Float)
    seed = Column(Integer)
    resultsfile = Column(String(50))

    def __repr__(self):
        return "<Flow {}>".format(self.id)

class Flow(object):
    """A two dimensional geophysical fluid flow."""
    def __init__(self, **kwargs):
        """Set up the appropriate parameters and initial conditions."""

        # set defaults and parse args
        for key in DEFAULTS:
            setattr(self, key, kwargs.get(key, DEFAULTS[key]))

        # if we didn't get a seed, choose one
        if self.seed is None:
            self.seed = np.random.randint(999999999)

        # initialize the random number generator
        np.random.seed(self.seed)

        # x space dimensions
        x = np.linspace(-np.pi, np.pi, self.nx, endpoint=False)
        y = np.linspace(-np.pi, np.pi, self.ny, endpoint=False)
        self.xx, self.yy = np.meshgrid(x, y, indexing='ij')

        # k space dimensions
        kx = np.fft.fftfreq(self.nx, self.lx/self.nx)
        ky = np.fft.rfftfreq(self.ny, self.ly/self.ny)
        self.kkx, self.kky = np.meshgrid(kx, ky, indexing='ij')
        self.ksq = self.kkx**2 + self.kky**2
        self.ksq[0, 0] += 1.0e-15
        self.kmax = np.min(np.max(kx), np.max(ky))
        self.kxmax = np.max(kx)
        self.kymax = np.max(ky)

        # initial state
        self.psihat = np.zeros_like(self.ksq, dtype=complex)
        self.qhat = self.forcing(0.0) * self.initial_amp
        self.q = np.fft.irfft2(self.qhat)
        self.psihat = self.get_psihat_from_qhat(self.qhat)
        self.psi = np.fft.irfft2(self.psihat)

        # integrator method and initial conditions
        self.integrator = ode(self.rhs).set_integrator(self.method)
        self.results = [self.qhat]
        self.t = 0.0

        # filename for results
        self.resultsfile, self.h5out = self.setup_hdf5()
    

    def setup_hdf5(self):
        """Set up the hdf5 archive."""
        fname =  "flow-{}-{}-{}.h5".format(self.w, self.kappa, self.seed)
        h5out = h5py.File(self.resultsfile, "w")
        h5out.create_dataset('qhat', (None, self.qhat.shape[0], self.qhat.shape[1]), dtype='c16')
        h5out.create_dataset('psihat', (None, self.qhat.shape[0], self.qhat.shape[1]), dtype='c16')
        h5out.create_dataset('q', (None, self.q.shape[0], self.q.shape[1]))
        h5out.create_dataset('psi', (None, self.q.shape[0], self.q.shape[1]))
        h5out.create_dataset('times', (None))

        h5out.create_dataset('energy', (None))
        h5out.create_dataset('enstrophy', (None))
        
        return fname, h5out

    def integrate(self, tf):
        """Run the integrator from the current time to tf.

        Integrator parameters are defined in __init__
        """

        # make sure we're walking forward in time
        if tf < self.t + self.dt:
            print("won't integrate backward in time from %f to %f" %(self.t,
                                                                     tf))
            return

        y0 = self.munge(self.qhat)
        self.integrator.set_initial_value(y0, self.t)

        # do the integration
        while self.integrator.successful() and self.integrator.t < tf:
            self.integrator.integrate(self.integrator.t + self.dt)
            self.results.append(self.unmunge(self.integrator.y))

        # update object state
        self.qhat = self.unmunge(self.integrator.y)
        self.psihat = self.get_psihat_from_qhat(self.qhat)
        self.psi = np.fft.irfft2(self.psihat)
        self.q = np.fft.irfft2(self.qhat)
        self.t = self.integrator.t

        # serialize

        return


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

        This provides random phases in a k-space annulus. Should be more
        configurable.
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

        if isinstance(arg1, type(0.0)):
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
        """Calculate the domain-integrated energy."""

        uhat = psihat*self.kkx*(0.0+1.0j)
        vhat = psihat*self.kky*(0.0+1.0j)

        u = np.fft.irfft2(uhat)
        v = np.fft.irfft2(vhat)

        efield = u**2 + v**2
        return np.sum(efield)

    def calc_enstrophy(self, qhat):
        """Calculate the domain-integrated enstrophy."""

        q = np.fft.irfft2(qhat)
        return np.sum(q**2)

    def munge(self, qhat):
        """format a complex k-space field for odeint"""

        r = qhat.real
        i = qhat.imag
        z = np.array([r, i])
        return z.reshape(-1)

    def unmunge(self, munged):
        """Return the 1d real sequence to its 2d complex state"""

        z = munged.reshape((2, self.nx, int(self.ny/2+1)))
        r = z[0]
        i = z[1]
        return r + (0+1.0j)*i

    def serialize(self):
        """Create an object to go into the database."""

        metadata = FlowORM(nx=self.nx, ny=self.ny, lx=self.lx, ly=self.ly,
                           beta=self.beta, kappa=self.kappa, w=self.w,
                           method=self.method, dt=self.dt,
                           initial_amp=self.initial_amp,
                           seed=self.seed)

        return metadata
