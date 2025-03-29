import os
import re
from math import gcd
from ase.io import vasp
from math import pi, sqrt,exp
from numpy import cross, linalg,log, expm1
import numpy as np
from scipy.integrate import quad



class ThermalConductivityCalculator:
    def __init__(self,dimension="3D"):
        self.k = 1.380649e-23  # Boltzmann constant, J/K
        self.hbar = 1.0545718e-34  # Reduced Planck's constant, J·s
        self.Na = 6.02214076e23
        self.AMU_to_grams = 1.66053906660e-27
        self.dimension = dimension




    def tau_U(self, omega, T, gamma, Theta_D, M, v_s):
        # Constants
        k_B = self.k
        hbar = self.hbar
        b = 3  # Constant for Umklapp processes
        
        # Calculate A_U
        A_U = (gamma ** 2 * k_B * Theta_D) / (M * v_s ** 2)
        
        # Umklapp scattering rate
        exp_term = np.exp(-Theta_D / (b * T))
        inv_tau_U = A_U * omega ** 2 * T * exp_term
        return 1 / inv_tau_U

    def tau_N(self, omega, T, gamma, M, v_s):
        # Constants
        k_B = self.k
        hbar = self.hbar
        
        # Calculate A_N
        A_N = (k_B * gamma ** 2) / (hbar * M * v_s ** 2)
        
        # Normal scattering rate
        inv_tau_N = A_N * omega ** 2 * T ** 3
        return 1 / inv_tau_N



    def callaway_integrand(self, x, t_ph):
        try:
            return ( x ** 4 * exp(x)) /t_ph/ expm1(x)**2 
        except OverflowError as e:
            print(f"OverflowError in callaway_integrand: {e}")
            return np.nan

    def slack_integrand(self, x, t_ph):
        return x  * t_ph  if  self.dimension =="2D" else  x ** 2  * t_ph 


         

    def cahill_integrand(self, x):
        """
        Integrand function for Cahill thermal conductivity.

        Args:
            x: (float) (hbar * omega) / (k * T)

        Returns:
            (float) Value of the integrand at x
        """
        try:
            if self.dimension == "2D":
                return (x ** 2 * exp(x)) / expm1(x) ** 2
            else:
                return (x ** 3 * exp(x)) / expm1(x) ** 2
        except OverflowError as e:
            print(f"OverflowError in cahill_integrand: {e}")
            return np.nan
#return (x **2 * exp(x)) / (expm1(x) ** 2) if self.dimension =="2D" else (x ** 3 * exp(x)) / (expm1(x) ** 2)            
            

    def cahill_integrand_summation(self, v_i, T, theta):
        """
        Calculate the summation term for the Cahill thermal conductivity integrand model.

        Args:
            v_i: (float) sound velocity for the acoustic mode i (in m/s)
            T: (float) absolute temperature (in K)
            theta: (float) Debye temperature (in K)

        Returns:
            (float) Summation term for one acoustic mode i
        """
        integrand = lambda x: self.cahill_integrand(x)
        try:
            integral_result, error = quad(integrand, 0, theta / T)
            return v_i * (T / theta) ** 2 * integral_result
        except Exception as e:
            print(f"Error in cahill_integrand_summation: {e}")
            return np.nan


       
    def slack_low_temp(self, v_m, theta, T, t_c):
        """
        Calculate Slack thermal conductivity.

        Args:
            v_m (float): Sound velocity.
            theta (float): Debye temperature.
            T (float): Temperature
            t_c (float): Scattering rate.

        Returns:
            float: Slack thermal conductivity (in SI units, i.e. W/(m-K))
        """
        try:
            integrand = lambda x: self.callaway_integrand(x, t_c)
            integral_result, error = quad(integrand, 0, theta / T)
            return (self.k / (2 * pi ** 2 * v_m)) * (self.k * T / self.hbar) ** 3 * integral_result
        except Exception as e:
            print(f"Error in slack_low_temp: {e}")
            return np.nan
            

    def slack_high_temp(self, v_m, theta, T=300, t_c=1E+12):
        """
        Calculate Slack thermal conductivity.

        Args:
            v_m (float): Sound velocity.
            theta (float): Debye temperature.
            T (float): Temperature
            t_c (float): Scattering rate.

        Returns:
            float: Slack thermal conductivity (in SI units, i.e. W/(m-K))
        """
    
        integrand = lambda x: self.slack_integrand(x, t_c)
        integral_result, error = quad(integrand, 0, theta / T)
        return (self.k / (2 * pi ** 2 * v_m)) * (self.k * T / self.hbar) ** 3 * integral_result


    def cahill_thermalconduct(self, velocities, T, theta, n, V):
        """
        Calculate Cahill thermal conductivity using the integrand model.

        Args:
            velocities: (list of float) sound velocities for each acoustic mode (in m/s)
            T: (float) absolute temperature (in K)
            theta: (float) Debye temperature (in K)
            Ref: 10.1038/nature13184

        Returns:
            (float) Cahill thermal conductivity (in W/(m*K))
        """
        try:
            cahill_sum = cahill_sum = self.cahill_integrand_summation(velocities, T, theta)
            return (pi / 6) ** (1. / 3.) * self.k * ((n/V))  * cahill_sum  if self.dimension == "2D" else  (pi / 6) ** (1. / 3.) * self.k * ((n /V)) ** (2. / 3.)  * cahill_sum
        except Exception as e:
            print(f"Error in cahill_thermalconduct: {e}")
            return np.nan
            

    def callaway_thermal_conductivity(self, v_s, Theta_D, T, gamma, M):
        """
        Calculate thermal conductivity using the Callaway model for both 2D and 3D materials.
        
        Args:
            v_s (float): Average sound velocity (m/s).
            Theta_D (float): Debye temperature (K).
            T (float): Temperature (K).
            gamma (float): Grüneisen parameter.
            M (float): Average atomic mass (amu).

        Returns:
            float: Thermal conductivity (W/(m·K)).
        """
        omega_D = self.k * Theta_D / self.hbar
        k_B = self.k
        hbar = self.hbar
        M = M * 1.660539e-27  # Convert atomic mass units to kg

        def tau_U(omega):
            return self.tau_U(omega, T, gamma, Theta_D, M, v_s)

        def tau_N(omega):
            return self.tau_N(omega, T, gamma, M, v_s)

        def tau_eff(omega):
            inv_tau_eff = 1 / tau_U(omega) + (omega / omega_D) ** 2 / tau_N(omega)
            return 1 / inv_tau_eff

        def integrand(x):
            omega = x * k_B * T / hbar
            tau = tau_eff(omega)
            if self.dimension == "2D":
                numerator = x ** 3 * np.exp(x)
                denominator = (np.expm1(x)) ** 2
                return tau * numerator / denominator
            else:  # Default to 3D
                numerator = x ** 4 * np.exp(x)
                denominator = (np.expm1(x)) ** 2
                return tau * numerator / denominator

        # Integration limits
        x_max = Theta_D / T

        integral_result, error = quad(integrand, 0, x_max)

        if self.dimension == "2D":
            pre_factor = (k_B / (2 * np.pi * v_s)) * (k_B * T / hbar) ** 2
        else:  # Default to 3D
            pre_factor = (k_B / (2 * np.pi ** 2 * v_s)) * (k_B * T / hbar) ** 3

        kappa = pre_factor * integral_result
        return kappa



            
                            
    def debye_model(self, M, E, rho):
        """
        Calculate Debye thermal conductivity based on the dimensionality of the material.

        Args:
            M (float): Molecular mass (in atomic mass units).
            E (float): Young's modulus (in Pascals).
            rho (float): Density (in kg/m^3 for 3D or kg/m^2 for 2D).
            dim (str): Dimensionality of the material ('2D' or '3D', or '1D').

        Returns:
            float: Debye thermal conductivity (in SI units, i.e., W/(m·K)).
        """
        return 2.489e-11 * (M ** (-1/3)) * (E ** 0.5) * (rho ** (-1/2)) if self.dimension == "2D" else 2.489e-11 * (M ** (-1/3)) * (E ** 0.5) * (rho ** (-1/6))

    def constant_volume_hc(self, t_d, n,temp):
        """
        Calculate the molar heat capacity at constant volume.
        
        Args:
            N (int): Number of atoms in the system = n*Na.
            t_d (float): Debye temperature.
            temp (float): Current temperature.
        
        Returns:
            float: Molar heat capacity at constant volume in J/(mol·K).
        """
        
        def integrand(x):
            exm1 = expm1(x)
            if exm1 == 0:
                return 0
            else:
                return (x**3 * exp(x)) / exm1**2 if self.dimension == "2D" else  (x**4 * exp(x)) / exm1**2

        t_ratio = temp / t_d
        upper_limit = min(t_d / temp, 10)
        
        #integrand = lambda x: (x**4 * exp(x) / (exp(x) - 1)**2)
        if self.dimension == "2D":
            integral_result, _ = quad(integrand, 0, upper_limit)
            t_ratio = t_ratio**2
            const = 3
        else:
            integral_result, _ = quad(integrand, 0, upper_limit)
            t_ratio = t_ratio**3
            const = 9
        
        c_v = const * self.Na*n*self.k * t_ratio * integral_result #quad(integrand, 0, 1/t_ratio)[0]
        return c_v


    def debye_entropy(self, t_d, n, temp):
         
        # Debye integral for entropy => \[ S(T) = 3Nk\left(-\ln\left[1 - e^{-\frac{\Theta_D}{T}}\right] + 4\frac{T^3}{\Theta_D^3}\int_{0}^{\frac{\Theta_D}{T}}\frac{x^3}{e^x - 1}dx\right) \]
        # J/(mol·K)
        def integrand(x):
            return (x**2) / expm1(x) if self.dimension == "2D" else (x**3) / expm1(x)

                    
        if self.dimension == "2D":
            integral_result, _ = quad(integrand, 0, t_d/temp)
            localterm = 2 * (temp/t_d)**3
            const = 1.0
        else:
            integral_result, _ = quad(integrand, 0, t_d/temp)
            localterm = 4 * (temp/t_d)**3
            const = 3.0

        # Calculating entropy
        log_arg = 1. - exp(-t_d/temp)
        if log_arg <= 0:
            log_arg = 1e-10  # A small positive number to avoid log(0) or log(negative)
        s = const * n * self.Na * self.k * ( localterm * integral_result - log(log_arg) )
        return s

    def entropy(self,Theta_D, n, T):
        def cv_over_t(T_prime):
            if T_prime <= 0:
                return 0
            else:
                return self.constant_volume_hc(Theta_D, n,T_prime)/ T_prime   
        
        # Integrate C_v/T' from a small positive value near 0 to T to avoid division by zero
        S, _ = quad(cv_over_t, 1e-3, T)#, limit=100) 
        return S
            
    
    
    def slack_simple_model(self, M, theta, v_a, gamma, n, T): 
        """
        Calculate the simple Slack thermal conductivity

        References
            - DOI: 10.1007/0-387-25100-6_2 (Title: "High Lattice Thermal Conductivity Solids")
            - DOI: 10.1002/adfm.201600718 (Title: "Minimum Thermal Conductivity in Weak Topological Insulators with
            Bismuth-Based Stack Structure")

        Args:
            A_0: (float) constant, depends on Gruneisen parameter
            M: (float) average atomic mass
            theta: (float) Debye temperature (K)
            v_a: (float) (v_a)**3 is the volume per atom (Angstroms)
            gamma: (float) Gruneisen parameter
            n: (int) number of atoms in primitive cell
            T: (float) absolute temperature (K)

        Returns: (float) Slack thermal conductivity (in SI units, i.e. W/(m.K))

        """
        if self.dimension == '2D':
            powr = 1.
            scale_factor = 1E10
        else:
            powr = 2./3.
            scale_factor = 1.
        A_0 = 2.43*1.0E-8/(1.0 -0.514/gamma +0.228/gamma**2 ) # See Phys. Rev. 137, A128 (1965)   
        return (A_0 * M * theta ** 3 * v_a) / (gamma * n ** (powr) * T)*scale_factor


    def cahill_thermal_conductivity(self, n, V,v_l, *v_ts):
        """
        Calculate Cahill thermal conductivity.
        Args:
            n: (int) number of atoms in the unit cell
            V: (float) unit cell volume (in SI units, i.e., m^3)
            v_l: (float) longitudinal sound velocity (in SI units, i.e., m/s)
            *v_ts: (float) variable number of transverse sound velocities (in SI units, i.e., m/s)
                   Accepts either 1 or 2 velocities for 2D or 3D materials, respectively.

        Returns: (float) Cahill thermal conductivity (in W/(m·K))
        """
        # Sum the longitudinal and all transverse sound velocities
        #total_velocity = v_l + sum(v_ts)
        if self.dimension == "1D":
            total_velocity = v_l  # Only longitudinal velocity is considered for 1D
        else:
            total_velocity = v_l + sum(v_ts) 
            
        if self.dimension =="2D":
            factor = (n / V)
        else:
            factor = (n / V) ** (2. / 3.)
        thermal_conductivity = 1./2. * ((pi / 6) ** (1. / 3.)) * self.k * factor * total_velocity
        # k/2.48 * (rho_2D*Na/M/1E-3) *(V_ma)*1E3
        return thermal_conductivity
        

    def clarke_model(self, n, E, rho, m):
        """
        Calculate Clarke thermal conductivity for both 3D and 2D materials.

        Args:
            n: (int) number of atoms in primitive cell
            E: (float) Young's modulus (in SI units, i.e., Kgm(s)^(-2))
            rho: (float) density (in SI units, i.e., Kg/m^3 for 3D or kg/m^2 for 2D)
            m: (float) total mass per unit cell (in SI units, i.e., Kg)
            dim: (str) dimensionality of the material ('3D' or '2D')

        Returns:
            (float) Clarke thermal conductivity (in SI units, i.e., W/(m·K))
        """
        if self.dimension == '2D':
            return 0.87 * self.k/m * (np.abs(E) ** (1. / 2.)) *np.sign(E)* (rho ** (1. / 2.))
        else:
            return 0.87 * self.k /(m ** (2. / 3.)) * (np.abs(E) ** (1. / 2.))*np.sign(E) * (rho ** (1. / 6.))
            


    def internalenergy(self, t_d, n,temp):
        """
        Calculate the Internal energy.
        
        Args:
            N (int): Number of atoms in the system = n*Na.
            t_d (float): Debye temperature.
            temp (float): Current temperature.
        
        Returns:
            float: Internal energy in J.
        """
        
        def integrand(x):
            exm1 = expm1(x)
            if exm1 == 0:
                return 0
            else:
                return x**3 / exm1


        def integrand_2D(x):
            exm1 = expm1(x)
            if exm1 == 0:
                return 0
            else:
                return x**2 / exm1
                
            
        t_ratio = temp / t_d
        upper_limit = min(t_d / temp, 10)
        
 
        if self.dimension == "2D":
            integral_result, _ = quad(integrand_2D, 0, upper_limit)
            t_ratio = t_ratio**2
            const = 1
        else:
            integral_result, _ = quad(integrand, 0, upper_limit)
            t_ratio = t_ratio**3
            const = 9
        
        U = const * self.Na*n*self.k * t_ratio * integral_result
        return U



    def enthalpy(self, t_d, n,temp):
      return (self.internalenergy(t_d, n,temp) - temp*self.entropy(t_d, n, temp))*n/self.Na/1.602E-19

    
    def meanfreepath2(self, M, rho, theta, n, T, vm, gamma):
        """
        Calculate the mean free path of phonons in a material at a given temperature.

        Parameters:
        - M: Mass of the material in AMU.
        - rho: Density of the material in kilograms per cubic meter (kg/m^3) or (kg/m^2) for 2D.
        - theta: Debye temperature of the material in Kelvin (K).
        - n: (int) number of atoms in primitive cell
        - T: Temperature at which the mean free path is calculated, in Kelvin (K).
        - vm: Velocity of sound in the material in meters per second (m/s).
        - gamma: (float) Gruneisen parameter

        Returns:
        - meanpath: The mean free path of phonons in the material at temperature T, in meters (m).
        """

        molar_mass_kg = M * self.AMU_to_grams * self.Na # Ensure M is in kg/mol

        # Calculate volumetric heat capacity
        cv_volume = self.constant_volume_hc(theta, n, T) * rho / molar_mass_kg  # J/(m^3·K)
        #print("cv_volume ", cv_volume)

        # Calculate mean free path
        meanpath = 3 * self.slack_simple_model(M/n, theta, vm, gamma, n, T) / (vm * cv_volume)
        #slack_low_temp(self, vm, theta, T=300, t_c=1E+12)
        return meanpath


    def meanfreepath(self, M, rho, theta, n, T, vm, t_c):
        """
        Calculate the mean free path of phonons in a material at a given temperature.

        Parameters:
        - M: Mass of the material in AMU.
        - rho: Density of the material in kilograms per cubic meter (kg/m^3) or (kg/m^2) for 2D.
        - theta: Debye temperature of the material in Kelvin (K).
        - n: (int) number of atoms in primitive cell
        - T: Temperature at which the mean free path is calculated, in Kelvin (K).
        - vm: Velocity of sound in the material in meters per second (m/s).
        - gamma: (float) Gruneisen parameter
        - t_c : lifetime 

        Returns:
        - meanpath: The mean free path of phonons in the material at temperature T, in meters (m).
        """

        molar_mass_kg = M * self.AMU_to_grams * self.Na # Ensure M is in kg/mol

        # Calculate volumetric heat capacity
        cv_volume = self.constant_volume_hc(theta, n, T) * rho / molar_mass_kg  # J/(m^3·K)
        #print("cv_volume ", cv_volume)

        # Calculate mean free path
        meanpath = 3 * self.slack_low_temp(vm, theta, T, t_c) / (vm * cv_volume)
        #
        return meanpath
        
                                
    def compute_thermal_conductivity_over_temperature_range(self, v_m, theta, M, v_a, gamma, n, vol, start, interval, end, rho, t_c=1E+12,thickness=None):
        temperatures = np.arange(start, end + interval, interval)
        thermal_conductivities = []
        thermal_conductivities_slack_simple = []
        thermal_cal_cahill = []
        #callaway_thermal = []
        heat_capacity = []
        entropies =  []
        debye_entropies = []
        free_energies = []
        mpf = []

        for T in temperatures:
            conductivity = self.slack_low_temp(v_m, theta, T, t_c)
            conductivity_slack_simple = self.slack_simple_model(M/n, theta, v_a, gamma, n, T)
            conductivity_cahill = self.cahill_thermalconduct(v_m, T, theta,n,vol)
            
            conductivity_callaway = self.callaway_thermal_conductivity(v_m, theta, T, gamma, M)
            #print(conductivity_callaway)
            
            capacity = self.constant_volume_hc(theta, n,T)
            entropy = self.entropy(theta, n, T)
            debye_entropy = self.debye_entropy(theta, n, T)
            free_energy = self.enthalpy(theta,n,T)
                       
            meanpath = self.meanfreepath(M, rho, theta, n, T, v_m, t_c) * (thickness if thickness is not None else 1) * 1E+6
            heat_capacity.append(capacity)
            thermal_conductivities.append(conductivity)
            thermal_conductivities_slack_simple.append(conductivity_slack_simple)
            #callaway_thermal.append(conductivity_callaway)
            entropies.append(entropy)
            debye_entropies.append(debye_entropy)
            free_energies.append(free_energy)
            thermal_cal_cahill.append(conductivity_cahill)
            mpf.append(meanpath)
            


        return temperatures, thermal_conductivities, thermal_conductivities_slack_simple, thermal_cal_cahill, heat_capacity, entropies, debye_entropies,free_energies,mpf

