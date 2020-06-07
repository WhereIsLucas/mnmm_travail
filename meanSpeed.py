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
speeds = []
standardDev = []

# this counts the number of frames
for k in range(1, 9) :
    base_path = "./datas/"+ str(k) +"00 pi/"
    path = base_path+"datas/"
    num_files = len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])
    print(num_files)
    totalFrames = num_files - 1
    energy = []

    for i in range(0, totalFrames):
        fileName = path + "grain" + str(i) + ".txt"
        data = np.genfromtxt(fileName,
                             delimiter=',',
                             dtype=types,
                             names=["ID", 'x', 'y', 'vx', 'vy', 'omega', 'radius'])
        velocity = 0
        if data.size > 1:
            for j in range(0, data.size):
                velocity += data[j]['vx'] + data[j]['vy']
        energy.insert(i, velocity)

    name = "velocity | m = " + str(np.mean(energy)) + " sd = " +  str(np.std(energy))
    plt.legend([name])
    # plt.plot(kineticEnergyRot)
    plt.savefig(base_path +'/velocity.png')

    speeds.insert(k, np.mean(energy))
    standardDev.insert(k, np.std(energy))
    plt.clf()
    plt.cla()

plt.xlabel('x')
plt.ylabel('y')
# plt.plot(speeds)
plt.plot(standardDev)
plt.legend(["standardDev"])
# plt.plot(kineticEnergyRot)
plt.savefig("./exports/" +'velocities.png')