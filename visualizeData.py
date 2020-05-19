import os.path
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

types = ['float', 'float', 'float', 'float', 'float', 'float', 'float']
barrelTypes = ['float', 'float', 'float', 'float', 'float', 'float']
domainTypes = ['float', 'float']
lineTypes = ['float', 'float']
data = []

# Set up the codec for the video file
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Lucas H'), bitrate=1800)
fig = plt.figure(figsize=(7, 7))

# this counts the number of frames
path = "./cmake-build-debug/datas/"
num_files = len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])
print(num_files)
totalFrames = num_files - 3
for i in range(0, totalFrames):
    fileName = path + "barrel" + str(i) + ".txt"
    data.insert(i, np.genfromtxt(fileName,
                                 delimiter=',',
                                 dtype=barrelTypes,
                                 names=['x', 'y', 'vx', 'vy', 'theta', 'radius']))
showingFrame = 0
scat = plt.scatter(data[showingFrame]["x"], data[showingFrame]['y'], alpha=0.5, s=data[0]['radius'] * 3.5 * 10000)
plt.title('Scatter plot test')
# plt.gca().set_aspect('equal', adjustable='box')
plt.axis("equal")
domain = np.genfromtxt(path + "domain.txt",
                       delimiter=',',
                       dtype=domainTypes,
                       names=['x', 'y'])

plt.xlim(-0.2 * domain["x"], domain['x'] * 1.2)
plt.ylim(-0.2 * domain["y"], domain['y'] * 1.2)
plt.xlabel('x')
plt.ylabel('y')
# a_circle = plt.Circle((.3, .3), 1.)
# plt.gca().add_artist(a_circle)

plan = np.genfromtxt(path + "plan.txt",
                     delimiter=',',
                     dtype=lineTypes,
                     names=['m', 'p'])

x = np.linspace(-0.2 * domain["x"], domain['x'] * 1.2, 100)
plt.plot(x, (plan["m"] * x) + plan["p"])


# plt.show()


def update(frame_number):
    scat.set_offsets(np.c_[data[frame_number]["x"], data[frame_number]["y"]])


animation = animation.FuncAnimation(fig, update, interval=10, frames=totalFrames)
# plt.show()
animation.save('exports/im.mp4', writer=writer)
