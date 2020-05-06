import numpy as np
import scipy.interpolate
from astropy.table import Table

class SNlightcurve:
    """

    # example

    # create the lightcurve object
    LC = SNlightcurve()
    # make an array of times to use when evaluating
    # this doesn't have to be the same as the input array of times
    # but it shouldn't be before the beginning or after the end
    # when including the shift & stretch
    t = np.arange(-10,80,0.5) + 10

    # you can also use the LC object in a fitting function
    # e.g.,
    # popt, pcov = scipy.optimize.curve_fit(LC.compute_lightcurve, t, m, sigma=m_err, p0 = [m0_guess, t0_guess, s_guess], absolute_sigma = True)
    # then follow the curve_fit documentation to interpret the output
    # to construct chi^2:
    # chisq = (((LC.compute_lightcurve(t, m0, t0, s) - m)/m_err)**2).sum()
    # dof = len(m) - 3    

    plt.clf()
    plt.plot(t,LC.compute_lightcurve(t, 15, 10, 1.1),label='stretch=1.1')
    plt.plot(t,LC.compute_lightcurve(t, 15, 10, 1.0),label='stretch=1.0')
    plt.plot(t,LC.compute_lightcurve(t, 15, 10, 0.9),label='stretch=0.9')
    plt.plot(10, 15, 'ro',label='max light')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.xlabel('Time (days)')
    plt.ylabel('Magnitude')

    """


    file = 'SNIa_lc_template.dat'

    def __init__(self, band='B'):
        
        # read in the data (once)
        self.data = Table.read(self.file,
                               format='ascii.commented_header')

        # make sure that the band we picked
        # is one of the available ones
        assert band.upper() in self.data.colnames[1:]

        # ignore the first point - it appears bad
        self.t = self.data['epoch'][1:]
        self.band = band
        self.mag = self.data[self.band][1:]
        
        # make an interpolating function to use for
        # going smoothly between the available points
        self.mag_interp = scipy.interpolate.interp1d(self.t,
                                                     self.mag,
                                                     )

        
    def compute_lightcurve(self,
                           t,
                           m0,
                           t0,
                           s):
        """
        lightcurve = compute_lightcurve(self,t,m0,t0,s)

        return the lightcurve evaluated at times t (in days)
        shifted to t0 (in days)
        shifted to magnitude m0
        and stretched by a factor s
        """

        return m0 + self.mag_interp((t - t0)/s)
                        
        

