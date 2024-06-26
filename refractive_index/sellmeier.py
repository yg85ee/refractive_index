# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 14:56:49 2024

@author: ygkim
"""

import numpy as np

class BK7:
    """
    BK7 is a class to calculate the refractive index of BK7 glass at a given wavelength.
    It includes methods to calculate group velocity dispersion and other properties.
    
    Methods:
        __call__(wl): Calculate the refractive index at wavelength wl.
        ng(wl): Calculate the group velocity at wavelength wl.
        GVD(wl): Calculate the group veolcity dispersion at wavelength wl in fs^2/mm.
        D(wl): Calculate the dispersion parameter D at wavelength wl in ps/(nm km).
    """
    def __call__(self,wl):
        """
        Calculate the refractive index of BK7 glass using the Sellmeier equation.

        Parameters:
            wl (float): Wavelength in meters.

        Returns:
            float: Refractive index of BK7 at the given wavelength.
        """
        x= wl*1e6
        n=(1+1.03961212/(1-0.00600069867/x**2)+0.231792344/(1-0.0200179144/x**2)+1.01046945/(1-103.560653/x**2))**.5
        return n
    
    def ng(self,wl):
        """
        Calculate the group refractive index of BK7 glass using the Sellmeier equation.

        Parameters:
            wl (float): Wavelength in meters.

        Returns:
            float: group index of BK7 at the given wavelength.
        """
        disp=(self.__call__(wl+0.00001*wl)-self.__call__(wl-0.00001*wl))/(0.00002*wl)
        n0=self.__call__(wl)
        return n0-wl*disp
    
    def GVD(self,wl):
        """
        Calculate the group veolcity dispersion of BK7 glass using the Sellmeier equation in fs^2/mm.

        Parameters:
            wl (float): Wavelength in meters.

        Returns:
            float: group velocity dispersion of BK7 at the given wavelength in fs^2/mm.
        """
        c0 = 299792458.0
        tmp=wl**3/(2*np.pi*c0**2)*((self.__call__(wl+0.00002*wl)-self.__call__(wl))-(self.__call__(wl)-self.__call__(wl-0.00002*wl)))/(0.00002*wl)**2
        return tmp*1e30/1e3

    def D(self,wl):
        """
        Calculate the dispersion parameter of BK7 glass using the Sellmeier equation in ps/(nm km).

        Parameters:
            wl (float): Wavelength in meters.

        Returns:
            float: the dispersion parameter of BK7 at the given wavelength in ps/(nm km).
        """
        c0 = 299792458.0
        tmp=wl**3/(2*np.pi*c0**2)*((self.__call__(wl+0.00002*wl)-self.__call__(wl))-(self.__call__(wl)-self.__call__(wl-0.00002*wl)))/(0.00002*wl)**2
        tmp*=-2*np.pi*c0/wl**2
        return tmp*1e6

if __name__ == '__main__':
    bk7 = BK7()
    print(bk7(0.5876e-6))
    print(bk7.ng(0.5876e-6))
    print(bk7.GVD(0.5876e-6))
    print(bk7.D(0.5876e-6))
    
    