import math
import os.path
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

types = ['float', 'float', 'float', 'float', 'float', 'float', 'float']
domainTypes = ['float', 'float']

# Set up the codec for the video file
Writer = animation.writers['ffmpeg']
writer = Writer(fps=25, metadata=dict(artist='Lucas H'), bitrate=1800)
fig = plt.figure(figsize=(7, 7))

# this counts the number of frames
path = "./cmake-build-debug/datas/"
num_files = len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])
print(num_files)
totalFrames = num_files - 1
energy = []
potentialEnergy = []
kineticEnergy = []
kineticEnergyRot = []
for i in range(0, totalFrames):
    fileName = path + "grain" + str(i) + ".txt"
    data = np.genfromtxt(fileName,
                         delimiter=',',
                         dtype=types,
                         names=["ID", 'x', 'y', 'vx', 'vy', 'theta', 'radius'])
    eTotal = 0
    potentialEnergyTot = 0
    kineticEnergyTot = 0
    kineticEnergyRotTot = 0
    if data.size > 1:
        for j in range(0, data.size):
            mass = math.pi * pow(data[j]['radius'],2)
            potentialEnergyTot += data[j]['y'] * mass * 9.81
            eTotal += potentialEnergyTot
            kineticEnergyTot += .5 * mass * (data[j]['vx'] * data[j]['vx'] + data[j]['vy'] * data[j]['vy'])
            eTotal += kineticEnergyTot
            kineticEnergyRotTot += .5 * .5 * mass * pow(data[j]['radius'], 2) * pow(data[j]['theta'], 2)
            eTotal += kineticEnergyRotTot
    else:
        mass = math.pi * data['radius'] * data['radius']
        potentialEnergyTot += data['y'] * mass * 9.81
        eTotal += potentialEnergyTot
        kineticEnergyTot += .5 * mass * (data['vx'] * data['vx'] + data['vy'] * data['vy'])
        eTotal += kineticEnergyTot
        kineticEnergyRotTot += .5 * .5 * mass * pow(data['radius'], 2) * pow(data['theta'], 2)
        eTotal += kineticEnergyRotTot
    energy.insert(i, eTotal)
    potentialEnergy.insert(i, potentialEnergyTot)
    kineticEnergy.insert(i, kineticEnergyRotTot + kineticEnergyTot)
    # kineticEnergyRot.insert(i, kineticEnergyRotTot + kineticEnergyTot)

plt.xlabel('x')
plt.ylabel('y')
plt.plot(energy)
plt.plot(potentialEnergy)
plt.plot(kineticEnergy)
plt.legend(['energy', 'potentialEnergy', 'totalKineticEnergy'])
# plt.plot(kineticEnergyRot)


plt.savefig('exports/energy.png')
