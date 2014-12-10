"""
An optics calculator for circular accelerator machines such as synchrotrons and storage rings.

De 2014 Xaratustrah

"""

import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


class Component:
    def __init__(self, typ, name, length, strength, param=0):
        self.length = length
        self.strength = strength
        self.param = param
        self.name = name
        self.typ = typ.lower()

        if self.typ == 'qf':
            self.mat = self.m_qf(self.length, self.strength)

        if self.typ == 'qd':
            self.mat = self.m_qd(self.length, self.strength)

        if self.typ == 'b':
            self.mat = self.m_b(self.length, self.strength)

        if self.typ == 'eb':
            self.mat = self.m_eb(self.strength, self.param)

        if self.typ == 'd':
            self.mat = self.m_d(self.length)

    def m_d(self, s):
        return np.array([[1, s, 0, 0, 0],
                         [0, 1, 0, 0, 0],
                         [0, 0, 1, s, 0],
                         [0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 1]])

    def m_qf(self, s, k):
        if k >= 0:
            return np.zeros((5, 5))
        sq_abs_k = np.sqrt(np.abs(k))
        Omega = sq_abs_k * s
        return np.array([[np.cos(Omega), np.sin(Omega) / sq_abs_k, 0, 0, 0],
                         [-sq_abs_k * np.sin(Omega), np.cos(Omega), 0, 0, 0],
                         [0, 0, np.cosh(Omega), np.sinh(Omega) / sq_abs_k, 0],
                         [0, 0, sq_abs_k * np.sinh(Omega), np.cosh(Omega), 0],
                         [0, 0, 0, 0, 1]])

    def m_qd(self, s, k):
        if k <= 0:
            return np.zeros((5, 5))
        sq_abs_k = np.sqrt(np.abs(k))
        Omega = sq_abs_k * s
        return np.array([[np.cosh(Omega), np.sinh(Omega) / sq_abs_k, 0, 0, 0],
                         [sq_abs_k * np.sinh(Omega), np.cosh(Omega), 0, 0, 0],
                         [0, 0, np.cos(Omega), np.sin(Omega) / sq_abs_k, 0],
                         [0, 0, -sq_abs_k * np.sin(Omega), np.cos(Omega), 0],
                         [0, 0, 0, 0, 1]])

    def m_b(self, s, R):
        if R <= 0:
            return np.zeros((5, 5))
        chi = s / R
        return np.array([[np.cos(chi), R * np.sin(chi), 0, 0, R * (1 - np.cos(chi))],
                         [-np.sin(chi) / R, np.cos(chi), 0, 0, np.sin(chi)],
                         [0, 0, 1, s, 0],
                         [0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 1]])

    def m_eb(self, R, psi):
        if R <= 0:
            return np.zeros((5, 5))
        psi = np.deg2rad(psi)
        return np.array([[1, 0, 0, 0, 0],
                         [np.tan(psi / R), 1, 0, 0, 0],
                         [0, 0, 1, 0, 0],
                         [0, 0, -np.tan(psi) / R, 1, 0],
                         [0, 0, 0, 0, 1]])


class Cell:
    def __init__(self, filename):

        self.blau = '#3090C7'
        self.gruen = '#507642'
        self.rot = '#D23641'

        self.filename = filename
        with open(self.filename) as f:
            self.title = f.readline().replace('#', '')
        lat = np.genfromtxt(filename, dtype=None)
        self.component_list = []
        self.length_list = []
        self.component_number = np.size(lat)
        self.length_list_sum = [0]

        for i in range(np.size(lat)):
            if np.size(lat) == 1:
                c = Component(lat[()][0].decode("utf-8"), lat[()][1], lat[()][2], lat[()][3], lat[()][4])
            else:
                c = Component(lat[i][0].decode("utf-8"), lat[i][1], lat[i][2], lat[i][3], lat[i][4])
            self.component_list.append(c)
            self.length_list.append(c.length)
            self.length_list_sum.append(self.length_list_sum[i] + c.length)

    def beta_a_x(self, m):
        return np.sqrt(-(m[0][1] * m[1][1]) / (m[1][0] * m[0][0]))

    def beta_b_x(self, m):
        return -1 / self.beta_a_x(m) * m[0][1] / m[1][0]

    def beta_a_y(self, m):
        return np.sqrt(-(m[2][3] * m[3][3]) / (m[3][2] * m[3][3]))

    def beta_b_y(self, m):
        return -1 / self.beta_a_y(m) * m[2][3] / m[3][2]

    def d_a_x(self, m):
        return -m[1][4] / m[1][0]

    def d_a_y(self, m):
        return -m[0][0] * m[1][4] / m[1][0] + m[0][4]

    def total_matrix(self):
        m = np.identity(5)
        # for i in range(self.component_number):
        for i in range(len(self.component_list)):
            m = np.dot(m, self.component_list[i].mat)
        return m

    def plot_cell(self):
        bx, by = self.evolve_beta()
        plt.xlabel('s [m]')
        plt.plot(self.length_list_sum, bx, label='beta_x', color=self.gruen, linewidth=2)
        ax1 = plt.gca()
        ax1.set_ylabel('beta_x [m]', color=self.gruen)

        line = plt.gca().get_lines()[0]
        maximum = int(line.get_ydata().max()) - 1
        cnt = 0
        global rect  # makes python happy
        for i in range(self.component_number):
            if self.component_list[i].typ == 'qf':
                rect = mpatches.Rectangle((cnt, maximum), self.length_list[i], 0.5, ec=self.gruen, fill=False,
                                          linewidth=2,
                                          joinstyle='round')
            if self.component_list[i].typ == 'qd':
                rect = mpatches.Rectangle((cnt, maximum - 0.5), self.length_list[i], 0.5, ec=self.rot, fill=False,
                                          linewidth=2,
                                          joinstyle='round')
            if self.component_list[i].typ == 'b':
                rect = mpatches.Rectangle((cnt, maximum - 0.5), self.length_list[i], 1, ec=self.blau, fill=False,
                                          linewidth=2,
                                          joinstyle='round')
            if self.component_list[i].typ == 'd':
                rect = mpatches.Rectangle((cnt, maximum), self.length_list[i], 0.1, ec="black", color='black', )
            cnt = cnt + self.length_list[i]
            plt.gca().add_patch(rect)

        ax2 = ax1.twinx()
        ax2.plot(self.length_list_sum, by, label='beta_y', color=self.rot, linewidth=2)
        ax2.set_ylabel('beta_y [m]', color=self.rot)
        plt.grid(True)
        plt.title(self.title)

        plt.show()
        plt.savefig(os.path.splitext(self.filename)[0] + '.pdf')

    def evolve_beta(self):
        Bx0 = self.beta_a_x(self.total_matrix())
        By0 = self.beta_a_y(self.total_matrix())
        B0 = np.array([[Bx0, 0, 0, 0, 0],
                       [0, 1 / Bx0, 0, 0, 0],
                       [0, 0, By0, 0, 0],
                       [0, 0, 0, 1 / By0, 0],
                       [0, 0, 0, 0, 1]
        ])
        B = B0
        beta_x_list, beta_y_list = [], []
        for i in range(len(self.component_list)):
            M = self.component_list[i].mat
            B = np.dot(np.dot(M, B), M.T)
            beta_x_list.append(B[0][0])
            beta_y_list.append(B[3][3])
        beta_x_list.append(self.beta_b_x(self.total_matrix()))
        beta_y_list.append(self.beta_b_y(self.total_matrix()))
        return beta_x_list, beta_y_list
