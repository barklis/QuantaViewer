#Virtual spectrometer module - B. Klis, 2022

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

class Spectrometry:
    def __init__(self, ch_num, raw=None, ene_data=None, eff_data=None, spect=None, norm_spect=None):
        #Channels number for particular detector
        self.channels = ch_num
        #Raw data
        if raw is None:
            self.raw_data = np.zeros(ch_num, dtype=int)
        else:
            self.raw_data = raw
        #Energy calibration data
        if ene_data is None:
            self.energy_data = np.zeros(ch_num, dtype=float)
        else:
            self.energy_data = ene_data
        #Efficiency calibration data
        if eff_data is None:
            self.efficiency_data = np.zeros(ch_num, dtype=float)
        else:
            self.efficiency_data = eff_data
        #Gamma spectrum data
        if spect is None:
            self.spectra = np.zeros((ch_num, 2), dtype=float)
        else:
            self.spectra = spect
        #Normalized spectrum data
        if norm_spect is None:
            self.normalized_spectra = np.zeros((ch_num, 2), dtype=float)
        else:
            self.normalized_spectra = norm_spect

    def __add__(self, other):
        ch_num = self.channels
        r_data = np.zeros(ch_num)
        eff_dat = np.zeros(ch_num)
        spec = np.zeros((ch_num, 2))
        norm_spect = np.zeros((ch_num, 2))
        for i in range(ch_num):
            r_data[i] = self.raw_data[i] + other.raw_data[i]
            eff_dat[i] = self.efficiency_data[i]
            spec[i][0] = self.spectra[i][0]
            spec[i][1] = self.spectra[i][1] + other.spectra[i][1]
            norm_spect[i][0] = self.normalized_spectra[i][0]
            norm_spect[i][1] = self.normalized_spectra[i][1] + other.spectra[i][1]
        return Spectrometry(ch_num, r_data, eff_dat, spec, norm_spect)

    def __sub__(self, other):
        ch_num = self.channels
        r_data = np.zeros(ch_num)
        eff_dat = np.zeros(ch_num)
        spec = np.zeros((ch_num, 2))
        norm_spect = np.zeros((ch_num, 2))
        for i in range(ch_num):
            r_data[i] = abs(self.raw_data[i] - other.raw_data[i])
            eff_dat[i] = self.efficiency_data[i]
            spec[i][0] = self.spectra[i][0]
            spec[i][1] = abs(self.spectra[i][1] - other.spectra[i][1])
            norm_spect[i][0] = self.normalized_spectra[i][0]
            norm_spect[i][1] = abs(self.normalized_spectra[i][1] - other.spectra[i][1])
        return Spectrometry(ch_num, r_data, eff_dat, spec, norm_spect)

    def load_data_from_file(self, file_name):
        with open(file_name, 'r') as f:
            channel = 0
            for line in f:
                self.raw_data[channel] = int(line)
                channel += 1
        for i in range(self.channels):
            self.spectra[i,0] = float(i)
            self.spectra[i,1] = float(self.raw_data[i])

    def load_data_from_tka_file(self, file_name):
        with open(file_name, 'r') as f:
            channel = 0
            head = 0
            for line in f:
                if head == 0:
                    self.online_time = line
                    head += 1
                    continue
                elif head == 1:
                    self.measurement_time = line
                    head += 1
                else:
                    self.raw_data[channel] = int(line)
                    channel += 1
        for i in range(self.channels):
            self.spectra[i,0] = float(i)
            self.spectra[i,1] = float(self.raw_data[i])

    def calibrate_energy(self, function=None):
        if function == None:
            for i in range(self.channels):
                self.spectra[i,0] = self.energy_data[i]
            pass
        for i in range(self.channels):
            self.spectra[i,0] = function(float(i))
            self.energy_data[i] = function(float(i))

    def calibrate_efficiency(self, function):
        for i in range(self.channels):
            energy = self.spectra[i,0]
            self.efficiency_data[i] = function(energy)
        
    def calculate_normalized_spectra(self):
        for i in range(self.channels):
            self.normalized_spectra[i,0] = self.spectra[i,0]
            self.normalized_spectra[i,1] = ( 1 / self.efficiency_data[i] ) * self.spectra[i,1]

    def suming_method(self, ranges):
        num_ranges = int(len(ranges)/2)
        sum_ranges = np.zeros((num_ranges, 2))
        for i in range(num_ranges):
            sum_ranges[i,0] = (ranges[2*i] + ranges[2*i+1]) / 2
            for n in range(self.channels):
                energy = self.spectra[n,0]
                counts = self.spectra[n,1]
                if ranges[2*i] < energy < ranges[2*i+1]: 
                    sum_ranges[i,1] += counts
        return sum_ranges

    def print_raw_data(self):
        for i in range(self.channels):
            print('Ch:' + str(i) + ' Counts: ' + str(self.raw_data[i]))

    def print_spectrum_data(self):
        for i in range(self.channels):
            print('Energy:' + str(self.spectra[i][0]) + ' Counts: ' + str(self.spectra[i][1]))
    
    def print_normalized_spectrum(self):
        for i in range(self.channels):
            print('Energy:' + str(self.normalized_spectra[i][0]) + ' Counts: ' + str(self.normalized_spectra[i][1]))

    def plot_raw_data(self):
        plt.plot(range(self.channels), self.raw_data[:], color='black')
        plt.yscale('log')
        plt.xlabel('Channel  (number)', fontweight='bold')
        plt.ylabel('Number of counts (counts)', fontweight='bold')
        #plt.fill_between(range(self.channels), 0, self.raw_data[:], color='blue')
        plt.show()

    def plot_spectrum(self):
        plt.plot(self.spectra[:,0], self.spectra[:,1], color='black')
        plt.yscale('log')
        plt.xlabel('Photon energy  (keV)', fontweight='bold')
        plt.ylabel('Number of counts (counts)', fontweight='bold')
        #plt.fill_between(range(self.channels), 0, self.raw_data[:], color='blue')
        plt.show()

    def plot_normalized_spectrum(self):
        plt.plot(self.normalized_spectra[:,0], self.normalized_spectra[:,1], color='black')
        plt.yscale('log')
        plt.xlabel('Photon energy (keV)', fontweight='bold')
        plt.ylabel('Number of counts (counts)', fontweight='bold')
        #plt.fill_between(self.norm_spectity[:,0], 0, self.norm_spectity[:,1], color='green')
        plt.show()

    def plot_energy_calibration(self):
        plt.plot(range(self.channels), self.energy_data[:], color='black')
        #plt.ylim( 0.0 , 2500.0)
        plt.xlabel('Channel  (number)', fontweight='bold')
        plt.ylabel('Photon energy (keV)', fontweight='bold')
        plt.show()

    def plot_efficiency_calibration(self):
        plt.plot(self.spectra[:,0], self.efficiency_data[:], color='black')
        #plt.ylim( 0.0 , 0.5)
        plt.xlabel('Photon energy (keV)', fontweight='bold')
        plt.ylabel('Absolute full-energy peak efficiency (a.u.)', fontweight='bold')
        plt.show()

#Default function for energy calibration
def energy_function(a:float=0, b:float=0, c:float=0, d:float=0):
    def e_func(channel):
        return a + b * channel + c * channel **2 + d * channel **3
    return e_func

#Default function for efficiency calibration
def efficency_function(e_low:list, e_high:list, limit=500.0):
    def eff_func(energy):
        if energy < limit:
            ln_eff = e_low[0] + e_low[1] * np.log(energy) + e_low[2] * (np.log(energy) ** 2) + e_low[3] * (np.log(energy) ** 3) + e_low[4] * (np.log(energy) ** 4) + e_low[5] * (np.log(energy) ** 5)
        else:
            ln_eff = e_high[0] + e_high[1] * np.log(energy) + e_high[2] * (np.log(energy) ** 2) + e_high[3] * (np.log(energy) ** 3) + e_high[4] * (np.log(energy) ** 4) + e_high[5] * (np.log(energy) ** 5)
        return np.exp(ln_eff)
    return eff_func

#EXAMPLE usage
def main():
    print('Virtual spectrometer module - example of usage, based on HPGe detector data')

    #Energy calibration data - linear function
    ene_function = energy_function(-1.249e-1, 3.284e-1)

    #Efficiency calibration data - fitted logarythmic polynomials
    eff_low = np.array([0, 0, 0, 0, 0, 0])
    eff_high = np.array([-9.29e1, 6.832e1, -2.066e1, 3.167e0, -2.491e-1, 8.023e-3])
    limit = 0 # Split point for dual calibration
    eff_function = efficency_function(eff_low, eff_high, limit)

    #Create virtual device
    spectrum = Spectrometry(8192) #Define channels number

    #Load data
    spectrum.load_data_from_tka_file('example_gamma_spec.TKA')

    #Calibrate data
    spectrum.calibrate_energy(ene_function)
    spectrum.calibrate_efficiency(eff_function)
    spectrum.calculate_normalized_spectra()

    #Printing data
    spectrum.print_raw_data()
    spectrum.print_spectrum_data()
    spectrum.print_normalized_spectrum()

    #Plot data
    spectrum.plot_raw_data()
    spectrum.plot_spectrum()
    spectrum.plot_normalized_spectrum()

    #Plot calibration
    spectrum.plot_energy_calibration()
    spectrum.plot_efficiency_calibration()

if __name__ == '__main__': main()