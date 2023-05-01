"""Run Helmholtz simulations."""

from helmholtz import helmsolve, wasserstein_rms, get_spectra
from geometry import basic_room, speaker_radial_boundary
import numpy as np
from numpy import pi
from scipy.optimize import basinhopping, brute
from scipy.optimize import Bounds


def f_range(f1, f2, N=10):
    f1, f2 = np.log10([f1, f2])
    return np.logspace(f1, f2, N)


def speaker_optimization_problem(geofunc=basic_room, f=[20, 40], N=100):
    def objective_function(U):
        """_summary_

        Parameters
        ----------
        U : arraylike
            For a single speaker U = [x, y, theta].

        Returns
        -------
        Wrms : float
            Wasserstein distance RMS.
        """
        
        *pos, theta = U
        
        # find solution
        ph = helmsolve(basic_room, pos, theta, f)
        
        # convert to frequency spectra
        spectra = get_spectra(ph)
        
        # use random sample
        np.random.shuffle(spectra)
        
        print("Solved Helmholtz.")
        
        # return w_rms    
        return wasserstein_rms(spectra[:N, :])
    
    return objective_function


def rectangular_bounds(x, y, theta=[-pi/4, pi/4]):
    U = np.array([x, y, theta])

    R = speaker_radial_boundary()
    U[:2, :] += np.array([[R, -R]])
    
    lb, ub = U[:, 0], U[:, 1]
    return U.tolist(), Bounds(lb, ub)


def basinhop_optimization(geofunc=basic_room, U0=[2, 5, 0], f=[20, 40],
                         I_max=1, N_sample=100,
                         x_bound=[0,5], y_bound=[0, 5]):
    func = speaker_optimization_problem(geofunc, f, N_sample)
    _, bounds = rectangular_bounds(x_bound, y_bound)
    min_kw = {"bounds": bounds}
    sol = basinhopping(func, U0, niter=I_max, minimizer_kwargs=min_kw)
    return sol


def brute_optimization(geofunc=basic_room, f=[20, 40],
                       N_splits=3, N_sample=50,
                       x_bound=[0,5], y_bound=[0, 5], 
                       plot=False, warp=0):
    func = speaker_optimization_problem(geofunc, f, N_sample)
    range, _ = rectangular_bounds(x_bound, y_bound)
    sol = brute(func, range, Ns=N_splits, finish=None)
    
    if plot:
        helmsolve(pos=sol[:2], angle=sol[2], f=f, plot=True, warp=warp)
    
    return sol
