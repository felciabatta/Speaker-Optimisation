"""Run Helmholtz simulations."""

from helmholtz import helmsolve, wasserstein_rms, get_spectra
from geometry import basic_room
import numpy as np
from numpy import pi


def speaker_optimization_problem(geofunc=basic_room, f=[20], N=100):
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
        
        # return w_rms    
        return wasserstein_rms(spectra[:N, :])
    
    return objective_function
