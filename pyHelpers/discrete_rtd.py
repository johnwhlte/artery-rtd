import numpy as np
import matplotlib.pyplot as plt

num_particles = 30000
deltaT = 0.001
data_file = "success_3"

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

t = np.arange(0, deltaT*len(e_t), deltaT)

print(f"Area under the E curve  -- > {np.sum(e_t*deltaT)}")
print(f"Mean Residence Time  -- > {np.trapz(t*e_t, t, deltaT)}")
print(f"Space time (V/Q)  -- > {1000000*0.02*(np.pi*(0.005**2)) / 1.985}")

plt.plot(t, e_t)
#plt.xlim((0,2))
plt.savefig("test.png")   
