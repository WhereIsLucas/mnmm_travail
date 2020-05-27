import math
import os.path
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

types = ['float', 'float', 'float', 'float', 'float', 'float', 'float']
typesBarrel = ['float', 'float', 'float', 'float', 'float', 'float']
domainTypes = ['float', 'float']

# Set up the codec for the video file
Writer = animation.writers['ffmpeg']
writer = Writer(fps=25, metadata=dict(artist='Lucas H'), bitrate=1800)
fig = plt.figure(figsize=(7, 7))

# this counts the number of frames
path = "./cmake-build-debug/datas/grain/"
pathBarrel = "./cmake-build-debug/datas/barrel/"
num_files = len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])
print(num_files)
totalFrames = num_files
energy = []
energyGrain = []
energyBarrel = []
potentialEnergy = []
kineticEnergy = []

for i in range(0, totalFrames):
    fileName = path + "grain" + str(i) + ".txt"
    fileNameBarrel = pathBarrel + "barrel" + str(i) + ".txt"
    data = np.genfromtxt(fileName,
                         delimiter=',',
                         dtype=types,
                         names=["ID", 'x', 'y', 'vx', 'vy', 'omega', 'radius'])
    data2 = np.genfromtxt(fileNameBarrel,
                          delimiter=',',
                          dtype=typesBarrel,
                          names=['x', 'y', 'vx', 'vy', 'omega', 'radius'])
    eTotal = 0
    potentialEnergyTot = 0
    kineticEnergyTot = 0
    kineticEnergyRotTot = 0
    eGrain = 0
    eBarrel = 0
    if data.size > 1:
        for j in range(0, data.size):
            mass = math.pi * pow(data[j]['radius'], 2)
            potentialEnergyValue = data[j]['y'] * mass * 9.81
            kineticEnergyValue = .5 * mass * (data[j]['vx'] * data[j]['vx'] + data[j]['vy'] * data[j]['vy'])
            kineticEnergyRotValue = .5 * .5 * mass * pow(data[j]['radius'], 2) * pow(data[j]['omega'], 2)

            eGrain += potentialEnergyValue + kineticEnergyValue + kineticEnergyRotValue
            eTotal += eGrain
            potentialEnergyTot += potentialEnergyValue
            kineticEnergyTot += kineticEnergyValue
            kineticEnergyRotTot += kineticEnergyRotValue
    else:
        mass = math.pi * data['radius'] * data['radius']
        potentialEnergyValue = data['y'] * mass * 9.81
        kineticEnergyValue = .5 * mass * (data['vx'] * data['vx'] + data['vy'] * data['vy'])
        kineticEnergyRotValue = .5 * .5 * mass * pow(data['radius'], 2) * pow(data['omega'], 2)

        eGrain += potentialEnergyValue + kineticEnergyValue + kineticEnergyRotValue
        eTotal += eGrain
        potentialEnergyTot += potentialEnergyValue
        kineticEnergyTot += kineticEnergyValue
        kineticEnergyRotTot += kineticEnergyRotValue

    mass = 2. * math.pi * data2['radius'] * 0.01 * 2
    potentialEnergyValue = data2['y'] * mass * 9.81
    kineticEnergyValue = .5 * mass * (data2['vx'] * data2['vx'] + data2['vy'] * data2['vy'])
    kineticEnergyRotValue = 2 / 3 * mass * pow(data2['radius'], 2) * pow(data2['omega'], 2)
    eBarrel += potentialEnergyValue + kineticEnergyValue + kineticEnergyRotValue
    eTotal += eBarrel
    potentialEnergyTot += potentialEnergyValue
    kineticEnergyTot += kineticEnergyValue
    kineticEnergyRotTot += kineticEnergyRotValue

    energy.insert(i, eTotal)
    energyGrain.insert(i, eGrain)
    energyBarrel.insert(i, eBarrel)
    potentialEnergy.insert(i, potentialEnergyTot)
    kineticEnergy.insert(i, kineticEnergyRotTot + kineticEnergyTot)

plt.xlabel('x')
plt.ylabel('y')
# plt.plot(energy)
plt.plot(energyGrain)
plt.plot(energyBarrel)
plt.plot(potentialEnergy)
plt.plot(kineticEnergy)
# plt.legend(['energy', 'energyGrain', 'energyBarrel'])
plt.legend(['energyGrain', 'energyBarrel', 'potentialEnergy', 'totalKineticEnergy'])
# plt.plot(kineticEnergyRot)


plt.savefig('exports/energy.png')
