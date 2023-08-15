import math
import os.path
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

types = ['float', 'float', 'float', 'float', 'float', 'float', 'float']
typesBall = ['float', 'float', 'float', 'float', 'float', 'float']
domainTypes = ['float', 'float']

# Set up the codec for the video file
Writer = animation.writers['ffmpeg']
writer = Writer(fps=25, metadata=dict(artist='Lucas H'), bitrate=1800)
fig = plt.figure(figsize=(7, 7))

# this counts the number of frames
path = "./cmake-build-debug/datas/grain/"
pathBall = "./cmake-build-debug/datas/ball/"
num_files = len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])
print(num_files)
totalFrames = num_files
energy = []
energyGrain = []
energyBall = []
potentialEnergy = []
kineticEnergy = []
kineticRotEnergy = []

for i in range(0, totalFrames):
    fileName = path + "grain" + str(i) + ".txt"
    fileNameBall = pathBall + "ball" + str(i) + ".txt"
    data = np.genfromtxt(fileName,
                         delimiter=',',
                         dtype=types,
                         names=["ID", 'x', 'y', 'vx', 'vy', 'omega', 'radius'])
    data2 = np.genfromtxt(fileNameBall,
                          delimiter=',',
                          dtype=typesBall,
                          names=['x', 'y', 'vx', 'vy', 'omega', 'radius'])
    eTotal = 0
    potentialEnergyTot = 0
    kineticEnergyTot = 0
    kineticEnergyRotTot = 0
    eGrain = 0
    eBall = 0
    if data.size > 1:
        for j in range(0, data.size):
            mass = math.pi * pow(data[j]['radius'], 2) * 1600
            potentialEnergyValue = data[j]['y'] * mass * 9.81
            kineticEnergyValue = .5 * mass * (data[j]['vx'] * data[j]['vx'] + data[j]['vy'] * data[j]['vy'])
            kineticEnergyRotValue = .5 * .5 * mass * pow(data[j]['radius'], 2) * pow(data[j]['omega'], 2)

            eGrain += potentialEnergyValue + kineticEnergyValue + kineticEnergyRotValue
            potentialEnergyTot += potentialEnergyValue
            kineticEnergyTot += kineticEnergyValue
            kineticEnergyRotTot += kineticEnergyRotValue
    else:
        mass = math.pi * data['radius'] * data['radius']
        potentialEnergyValue = data['y'] * mass * 9.81
        kineticEnergyValue = .5 * mass * (data['vx'] * data['vx'] + data['vy'] * data['vy'])
        kineticEnergyRotValue = .5 * .5 * mass * pow(data['radius'], 2) * pow(data['omega'], 2)

        eGrain += potentialEnergyValue + kineticEnergyValue + kineticEnergyRotValue
        potentialEnergyTot += potentialEnergyValue
        kineticEnergyTot += kineticEnergyValue
        kineticEnergyRotTot += kineticEnergyRotValue

    eTotal += eGrain

    mass = math.pi * pow(data2['radius'], 2) * 1400 * 0.05
    potentialEnergyValue = data2['y'] * mass * 9.81
    kineticEnergyValue = .5 * mass * (data2['vx'] * data2['vx'] + data2['vy'] * data2['vy'])
    kineticEnergyRotValue = 2 / 3 * mass * pow(data2['radius'], 2) * pow(data2['omega'], 2)
    eBall = potentialEnergyValue + kineticEnergyValue + kineticEnergyRotValue
    eTotal += eBall
    potentialEnergyTot += potentialEnergyValue
    kineticEnergyTot += kineticEnergyValue
    kineticEnergyRotTot += kineticEnergyRotValue

    energy.insert(i, eTotal)
    energyGrain.insert(i, eGrain)
    energyBall.insert(i, eBall)
    potentialEnergy.insert(i, potentialEnergyTot)
    kineticEnergy.insert(i, kineticEnergyTot)
    kineticRotEnergy.insert(i, kineticEnergyRotTot)
plt.figure(1)
plt.xlabel('x')
plt.ylabel('y')
# plt.plot(energy)
# plt.plot(energyGrain)
# plt.plot(energyBall)
plt.plot(potentialEnergy)
plt.plot(kineticRotEnergy)
plt.plot(kineticEnergy)
# plt.legend(['energyGrain', 'energyBall'])
plt.legend(['potential energy', 'kinetic rot energy', 'kinetic energy'])
# plt.plot(kineticEnergyRot)


plt.savefig('exports/energy.png')

plt.figure(2)
plt.xlabel('x')
plt.ylabel('y')
# plt.plot(energy)
plt.plot(energyGrain)
plt.plot(energyBall)
# plt.plot(potentialEnergy)
# plt.plot(kineticEnergy)
# plt.legend(['energy'])
plt.legend(['energyGrain', 'energyBall'])
# plt.plot(kineticEnergyRot)
#
#
plt.savefig('exports/GrainBall.png')
