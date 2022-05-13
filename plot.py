import matplotlib.pyplot as plt
# import numpy as np

Mx_data: list = list()
My_data: list = list()

Ex_data: list = list()
Ey_data: list = list()

Xx_data: list = list()
Xy_data: list = list()

Cx_data: list = list()
Cy_data: list = list()

try:
    print("Opening file...")
    fname = "mag_data.dat"
    with open(file=fname) as data_file:
        for i, cur_ln in enumerate(data_file):
            Mx_data.append(i)
            My_data.append(float(cur_ln))
except Exception:
    print("an error occured!")

try:
    print("Opening file...")
    fname = "energy_data.dat"
    with open(file=fname) as data_file:
        for i, cur_ln in enumerate(data_file):
            Ex_data.append(i)
            Ey_data.append(float(cur_ln))
except Exception:
    print("an error occured!")

try:
    print("Opening file...")
    fname = "susceptibility_data.dat"
    with open(file=fname) as data_file:
        for i, cur_ln in enumerate(data_file):
            Xx_data.append(i)
            Xy_data.append(float(cur_ln))
except Exception:
    print("an error occured!")

try:
    print("Opening file...")
    fname = "heatcap_data.dat"
    with open(file=fname) as data_file:
        for i, cur_ln in enumerate(data_file):
            Cx_data.append(i)
            Cy_data.append(float(cur_ln))
except Exception:
    print("an error occured!")

fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(Mx_data, My_data)
axs[0, 0].set_title('Magnitization v. Iteration')
axs[0, 1].plot(Ex_data, Ey_data, 'tab:orange')
axs[0, 1].set_title('Dimentionless Energy v. Iteration')
axs[1, 0].plot(Xx_data, Xy_data, 'tab:green')
axs[1, 0].set_title('Susceptibility v. Iteration')
axs[1, 1].plot(Cx_data, Cy_data, 'tab:red')
axs[1, 1].set_title('Heat Capacity v. Iteration')


# plt.scatter(Mx_data, My_data)
plt.show()
