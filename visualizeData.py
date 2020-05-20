import os.path
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

types = ['float', 'float', 'float', 'float', 'float', 'float', 'float']
barrelTypes = ['float', 'float', 'float', 'float', 'float', 'float']
domainTypes = ['float', 'float']
lineTypes = ['float', 'float']
data = []
dataGrain = []

# Set up the codec for the video file
Writer = animation.writers['ffmpeg']
writer = Writer(fps=25, metadata=dict(artist='Hugo T'), bitrate=1800)
fig = plt.figure(figsize=(7, 7))

# this counts the number of frames
pathBase = "./cmake-build-debug/datas/"
pathGrain = "./cmake-build-debug/datas/grain/"
pathBarrel = "./cmake-build-debug/datas/barrel/"

# Barrel
num_files = len([f for f in os.listdir(pathBarrel) if os.path.isfile(os.path.join(pathBarrel, f))])
print(num_files)
totalFrames = num_files
for i in range(0, totalFrames):
    fileName = pathBarrel + "barrel" + str(i) + ".txt"
    data.insert(i, np.genfromtxt(fileName,
                                 delimiter=',',
                                 dtype=barrelTypes,
                                 names=['x', 'y', 'vx', 'vy', 'theta', 'radius']))
showingFrame = 0
scat = plt.scatter(data[showingFrame]["x"], data[showingFrame]['y'], alpha=0.5, s=data[0]['radius'] * 2.4 * 1000,
                   facecolors="none", edgecolors="red")

# Grain
num_files = len([f for f in os.listdir(pathGrain) if os.path.isfile(os.path.join(pathGrain, f))])
print(num_files)
totalFrames = num_files
for i in range(0, totalFrames):
    fileName = pathGrain + "grain" + str(i) + ".txt"
    dataGrain.insert(i, np.genfromtxt(fileName,
                                      delimiter=',',
                                      dtype=types,
                                      names=['ID', 'x', 'y', 'vx', 'vy', 'theta', 'radius']))
showingFrame = 0
scatGrain = plt.scatter(dataGrain[showingFrame]["x"], dataGrain[showingFrame]['y'], alpha=0.5,
                        s=dataGrain[0]['radius'] * 2.8 * 100,
                        facecolors="none", edgecolors="blue")

plt.title('Scatter plot test')
# plt.gca().set_aspect('equal', adjustable='box')
plt.axis("equal")
domain = np.genfromtxt(pathBase + "domain.txt",
                       delimiter=',',
                       dtype=domainTypes,
                       names=['x', 'y'])

plt.xlim(-0.2 * domain["x"], domain['x'] * 1.2)
plt.ylim(-0.2 * domain["y"], domain['y'] * 1.2)
plt.xlabel('x')
plt.ylabel('y')
# a_circle = plt.Circle((.3, .3), 1.)
# plt.gca().add_artist(a_circle)

plan = np.genfromtxt(pathBase + "plan.txt",
                     delimiter=',',
                     dtype=lineTypes,
                     names=['m', 'p'])

x = np.linspace(-1 * domain["x"], domain['x'] * 2, 100)
plt.plot(x, (plan["m"] * x) + plan["p"])

# plt.show()
plt.savefig("exports/im.png")


#
def update(frame_number):
    scat.set_offsets(np.c_[data[frame_number]["x"], data[frame_number]["y"]])
    scatGrain.set_offsets(np.c_[dataGrain[frame_number]["x"], dataGrain[frame_number]["y"]])


animation = animation.FuncAnimation(fig, update, interval=40, frames=totalFrames)
# plt.show()
animation.save('exports/im.mp4', writer=writer)
