import numpy as np
import matplotlib.pyplot as plt

num_particles = 5000
deltaT = 0.005
data_file = "run_log"
comp_data_file = "dataFVM/2ml.dat"

with open(data_file,"r") as f:

    lines = f.readlines()

stats_line = []

for line in lines:
    line_split = line.split()

    if line_split[0] == '[captureStatistics]' and len(line_split) == 4:
        stats_line.append(int(line_split[3].split('=')[1]))

e_t = np.zeros(shape=np.asarray(stats_line).shape)

for i, line in enumerate(stats_line):
    if i == 0:
        summer = 0
    else:
        summer+= line
        e_t[i] = ((line - stats_line[i-1])/(num_particles*deltaT))

t = np.arange(0, deltaT*len(e_t)-0.5*deltaT, deltaT)
#t = np.arange(0, deltaT*len(e_t), deltaT)
print(t)
print(e_t)

with open(comp_data_file, "r") as f:

    lines = f.readlines()

C_t_comp = []
delta_Ts = []

for line in lines:
    C_t_comp.append(float(line.split()[1]))
    delta_Ts.append(float(line.split()[2]))

C_t_comp = np.asarray(C_t_comp)
e_t_comp = C_t_comp / np.trapz(C_t_comp, delta_Ts)

# e_t_comp = np.zeros(shape=C_t_comp.shape)

# for i, line in enumerate(C_t_comp):
#     if i!= 0:
#         e_t_comp[i] = ((line - C_t_comp[i-1])/(delta_Ts[i] - delta_Ts[i-1]))

print(f"Area under the E curve  -- > {np.sum(e_t*deltaT)}")
print(f"Mean Residence Time LBM -- > {np.trapz(t*e_t, t, deltaT)}")
print(f"Mean Residence Time FVM -- > {np.trapz(delta_Ts*e_t_comp, delta_Ts)}")
print(f"Space time (V/Q)  -- > {1000000*0.02*(np.pi*(0.005**2)) / 1.985}")

plt.plot(t, e_t, label="lbm")
plt.plot(delta_Ts, e_t_comp, label="fvm")
plt.xlabel("Time (s)")
plt.ylabel("E(t)")
plt.legend()
plt.savefig("test.png")   
