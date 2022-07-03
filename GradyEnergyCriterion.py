import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math


class GradyEnergyCriterion:
    def __init__(self, metal_data, strain_rate):
        self.R = metal_data[0]
        self.B0 = metal_data[1]*10**9  # Bulk modulus at standard state, Pa
        self.C0 = metal_data[2]
        self.Kmin = metal_data[3]*10**6  # fracture toughness N/m^(3/2)
        self.Kmax = metal_data[4]*10**6
        self.Kavg = (self.Kmin + self.Kmax)/2
        self.gamma = self.Kavg**2/(2*self.R*self.C0**2)  # Fracture surface energy per unit area

        self.Ymin = metal_data[5]*10**9  # flow stress
        self.Ymax = metal_data[6]*10**9
        self.Yavg = (self.Ymin + self.Ymax)/2

        self.Ec = metal_data[7]  # void fracture
        self.E_rate = strain_rate  # Estimated
        # print(self.R, self.B0,self.C0,self.Kmin, self.Kmax, self.Ymin, self.Ymax, self.Ec )

        Ec_rate = GradyEnergyCriterion.transition_strain(self)
        ductile_min = GradyEnergyCriterion.ductile_spall(self, self.Ymin)
        ductile_max = GradyEnergyCriterion.ductile_spall(self, self.Ymax)
        brittle_min = GradyEnergyCriterion.brittle_spall(self, self.Kmin)
        brittle_max = GradyEnergyCriterion.brittle_spall(self, self.Kmax)
        data = np.array([ductile_min, ductile_max, brittle_min, brittle_max])
        print(Ec_rate)
        print(data)
        print(gamma)

        GradyEnergyCriterion.plot_energy_unit_volume(self)

    def transition_strain(self):
        Ec_rate_min = math.sqrt((8*self.B0**2*(self.Ymin*self.Ec)**3)/(9*self.R*self.Kmax**4))
        Ec_rate_max = math.sqrt((8*self.B0**2*(self.Ymax*self.Ec)**3)/(9*self.R*self.Kmin**4))
        return np.array([Ec_rate_min, Ec_rate_max])

    def ductile_spall(self, Y):  # Flow stress dominated
        Ps = math.sqrt(2*self.R*(self.C0**2)*Y*self.Ec)  # Pa
        ts = math.sqrt(2*Y*self.Ec/(self.R*(self.C0**2)*(self.E_rate**2)))  # s
        s = 2*self.C0*ts  # m
        return np.array([Ps, ts, s])

    def brittle_spall(self, K):  # Fracture toughness dominated
        Ps = (3*self.R*self.C0*(K**2)*self.E_rate)**(1/3)  # Pa
        ts = (1/self.C0)*(math.sqrt(3)*K/(self.R*self.C0*self.E_rate))**(2/3)  # s
        s = 2*self.C0*ts  # m
        return np.array([Ps, ts, s])

    def plot_energy_unit_volume(self):
        t = np.linspace(0.01*10**-6, 10**-6, 101)
        fracture_E = 3*self.gamma/(self.C0*t)
        elastic_E = self.R*(self.C0*self.E_rate*t)**2/2

        plt.gca.set_xlim([t[0], t[100]])
        font = {'family': 'calibri', 'color': 'black', 'size': 15}
        plt.xlabel('time or correlation horizon, s', fontdict=font)
        plt.ylabel('Energy, J/m^3', fontdict=font)
        plt.plot(t, fracture_E, 'r', label='Fracture Surface Energy')
        plt.plot(t, elastic_E, 'b', label='Elastic strain Energy')
        plt.legend()
        plt.show()

        t0 = 0.1*10**-6
        t_intercept = fsolve(GradyEnergyCriterion.equilibrium_fragmentation, t0, args=(self.Kavg, self.R, self.C0, self.E_rate), xtol=1e-6)
        E_intercept = 3*gamma/(self.C0*t_intercept[0])
        print(t_intercept, E_intercept)

    def equilibrium_fragmentation(t, Kavg, R, C0, E_rate):
        gamma = Kavg**2/(2*R*C0**2)  # Fracture surface energy per unit area
        fracture_E = 3*gamma/(C0*t)
        elastic_E = R*(C0*E_rate*t)**2/2

        return fracture_E - elastic_E
