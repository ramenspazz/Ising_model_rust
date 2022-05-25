import matplotlib.pyplot as plt
# import numpy as np

My_data: list = list()

Ey_data: list = list()

Xy_data: list = list()

Cy_data: list = list()

beta_data: list = list()

try:
    print("Opening file...")
    fname = "beta.dat"
    with open(file=fname) as data_file:
        for i, cur_ln in enumerate(data_file):
            beta_data.append(float(cur_ln))
except Exception:
    print("an error occured!")

try:
    print("Opening file...")
    fname = "mag_data.dat"
    with open(file=fname) as data_file:
        for i, cur_ln in enumerate(data_file):
            My_data.append(float(cur_ln))
except Exception:
    print("an error occured!")

try:
    print("Opening file...")
    fname = "energy_data.dat"
    with open(file=fname) as data_file:
        for i, cur_ln in enumerate(data_file):
            Ey_data.append(float(cur_ln))
except Exception:
    print("an error occured!")

try:
    print("Opening file...")
    fname = "susceptibility_data.dat"
    with open(file=fname) as data_file:
        for i, cur_ln in enumerate(data_file):
            Xy_data.append(float(cur_ln))
except Exception:
    print("an error occured!")

try:
    print("Opening file...")
    fname = "heatcap_data.dat"
    with open(file=fname) as data_file:
        for i, cur_ln in enumerate(data_file):
            Cy_data.append(float(cur_ln))
except Exception:
    print("an error occured!")

fig, axs = plt.subplots(2, 2)
axs[0, 0].semilogx(beta_data, My_data)
axs[0, 0].set_title('Magnitization v. Beta')
axs[0, 1].semilogx(beta_data, Ey_data, 'tab:orange')
axs[0, 1].set_title('Dimentionless Energy v. Beta')
axs[1, 0].semilogx(beta_data, Xy_data, 'tab:green')
axs[1, 0].set_title('Susceptibility v. Beta')
axs[1, 1].semilogx(beta_data, Cy_data, 'tab:red')
axs[1, 1].set_title('Heat Capacity v. Beta')

plt.show()
