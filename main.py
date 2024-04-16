import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from astropy import modeling

HEL1_PATH = r'.\KvantMätningar\helspektrum_mörkt_1mmEntrance.xlsx'
HEL2_PATH = r'.\KvantMätningar\helspektrum_mörkt_2_1mmEntrance.xlsx'
HEL3_PATH = r'.\KvantMätningar\helspektrum_mörkt_3_1mmEntrance.xlsx'

TOPP3_1_PATH = r'.\KvantMätningar\topp3_1mm_mörkt_1.xlsx'
TOPP3_2_PATH = r'.\KvantMätningar\topp3_1mm_mörkt_2.xlsx'
TOPP3_3_PATH = r'.\KvantMätningar\topp3_1mm_mörkt_3.xlsx'

TOPP4_1_PATH = r'.\KvantMätningar\topp4_1.xlsx'
TOPP4_2_PATH = r'.\KvantMätningar\topp4_2.xlsx'
TOPP4_3_PATH = r'.\KvantMätningar\topp4_3.xlsx'

B1, A1 = (608,616), 1000
B2, A2 = (616,624), 600
B3, A3 = (624,632), 600
B4, A4 = (632,641), 500
B5, A5 = (640,650), 500
B6, A6 = (650,658), 500
B7, A7 = (658,668), 550
B8, A8 = (666,676), 550
B9, A9 = (676,686), 500
B10, A10 = (684,696), 500
B11, A11 = (696,704), 500
B12, A12 = (704,714), 500
B13, A13 = (726,734), 530
B14, A14 = (736,744), 520
Bs, As = [B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,], [A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14]

class Measurement:
    def __init__(self, path):
        data_dict = pd.read_excel(path).to_dict(orient='list')
        keys = list(data_dict.keys())
        self.data = {'Wavelength': np.array(data_dict[keys[0]])/1e3, 'Current': np.array(data_dict[keys[1]])}

    def table(self):
        print(pd.DataFrame.from_dict(self.data, orient='columns'))
    
    def __call__(self):
        plt.plot(self.data['Wavelength'], self.data['Current'], color='black')
        plt.xlabel('Wavelength')
        plt.ylabel('Current')
        plt.show()

    def fit_gauss(self, bounds):
        plt.plot(self.data['Wavelength'][bounds[0]:bounds[1]], self.data['Current'][bounds[0]:bounds[1]], color='black')
        def gauss(x, a, b, c):
            return a*np.exp((x-b)**2/(2*c**2))
        params = curve_fit(gauss, self.data['Wavelength'][bounds[0]:bounds[1]], self.data['Current'][bounds[0]:bounds[1]])[0]
        x = np.linspace(bounds[0], bounds[1], 1000)
        plt.plot(x, gauss(x, *tuple(params)))
        plt.show()

class Measurement_group:
    def __init__(self, *paths, offset=0):
        self.Measurements = []
        for path in paths:
            self.Measurements.append(Measurement(path))
        self.condition = True
        for i in range(len(self.Measurements)):
            if not self.Measurements[0].data['Wavelength'] in self.Measurements[i].data['Wavelength']:
                self.condition = False
        if self.condition:
            self.averaged_data = {'Wavelength': self.Measurements[0].data['Wavelength'], 'Current': sum(np.array([self.Measurements[i].data['Current'] for i in range(len(self.Measurements))]))/len(self.Measurements)+offset}
        self.find_peaks(Bs, As)

    def show_data(self):
        fig, axs = plt.subplots(1, len(self.Measurements))
        for i in range(len(self.Measurements)):
            axs[i].plot(self.Measurements[i].data['Wavelength'], self.Measurements[i].data['Current'], color='black')
            axs[i].set_xlabel('Wavelength')
            axs[i].set_ylabel('Current')
        plt.show()
    
    def __call__(self):
        plt.plot(self.averaged_data['Wavelength'], self.averaged_data['Current'], color='black')
        plt.show()

    def compare(self, condition=None):
        if not self.condition:
            print('Wavelengths do not match')
            return
        if condition==None:
            fig, axs = plt.subplots(1, len(self.Measurements)+1)
            axs[0].plot(self.averaged_data['Wavelength'], self.averaged_data['Current'], color='black')
            axs[0].set_xlabel('Wavelength')
            axs[0].set_ylabel('Current')
            axs[0].set_title('Averaged')
            for i in range(len(self.Measurements)):
                axs[i+1].plot(self.Measurements[i].data['Wavelength'], self.Measurements[i].data['Current'], color='black')
                axs[i+1].set_xlabel('Wavelength')
                axs[i+1].set_ylabel('Current')
    
        else:
            fig, axs = plt.subplots(1, len(self.Measurements))
            for i in range(len(self.Measurements)):
                axs[i].plot(self.averaged_data['Wavelength'], self.averaged_data['Current'], color='red')
                axs[i].plot(self.Measurements[i].data['Wavelength'], self.Measurements[i].data['Current'], '--k')
                axs[i].set_xlabel('Wavelength')
                axs[i].set_ylabel('Current')
        plt.xlabel('Wavelength')
        plt.ylabel('Current')
        plt.show()
    
    def fit_gauss(self, Bounds, **guess):
        bounds = self.get_closest_index(Bounds)
        fitter = modeling.fitting.LevMarLSQFitter()
        model = modeling.models.Gaussian1D(**guess, mean=sum(Bounds)/2, stddev=2)
        fitted_model = fitter(model, self.averaged_data['Wavelength'][bounds[0]:bounds[1]], self.averaged_data['Current'][bounds[0]:bounds[1]])
        # print(f'{fitted_model.amplitude.value:.3e}, {fitted_model.mean.value:.3e}, {fitted_model.stddev.value:.3e}')
        return fitted_model.amplitude.value, fitted_model.mean.value, fitted_model.stddev.value
    
    def get_closest_index(self, Bounds):
        wvls = self.averaged_data['Wavelength']
        index_b = []
        for Bound in Bounds:
            for i in range(len(wvls)):
                if Bound-wvls[i] >=0 and Bound-wvls[i+1] <= 0:
                    index_b.append(i)
        return index_b
    
    def plot_peaks(self, params_list=None):
        if params_list==None:
            params_list=self.params_list
        wvls = self.averaged_data['Wavelength']
        amps = self.averaged_data['Current']
        x = np.linspace(min(wvls), max(wvls), 10000)
        plt.plot(wvls, amps, '-k')
        for params in params_list:
            plt.plot(x, modeling.models.Gaussian1D(amplitude=params[0], mean=params[1], stddev=params[2])(x), '--r')
        plt.show()

    def find_peaks(self, Bs, As):
        self.params_list = []
        for i in range(len(Bs)):
            self.params_list.append(self.fit_gauss(Bounds=Bs[i], amplitude=As[i]))

HEL = Measurement_group(HEL1_PATH, HEL2_PATH, HEL3_PATH, offset=0)

HEL.plot_peaks()
