import os.path
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

types = ['float', 'float', 'float', 'float', 'float', 'float', 'float']
domainTypes = ['float', 'float']
data = []

# Set up the codec for the video file
Writer = animation.writers['ffmpeg']
writer = Writer(fps=25, metadata=dict(artist='Lucas H'), bitrate=1800)
fig = plt.figure(figsize=(7, 7))

# this counts the number of frames
path = "./cmake-build-debug/datas/"
num_files = len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])
print(num_files)
totalFrames = num_files - 1
for i in range(0, totalFrames):
    fileName = path + "grain" + str(i) + ".txt"
    data.insert(i, np.genfromtxt(fileName,
                                 delimiter=',',
                                 dtype=types,
                                 names=["ID", 'x', 'y', 'vx', 'vy', 'theta', 'radius']))
showingFrame = 0
scat = plt.scatter(data[showingFrame]["x"], data[showingFrame]['y'], alpha=0.5, s=pow(data[0]['radius'],2)*270000, c='red')
plt.title('Scatter plot test')
# plt.gca().set_aspect('equal', adjustable='box')
plt.axis("equal")
domain = np.genfromtxt(path + "/domain.txt",
                       delimiter=',',
                       dtype=domainTypes,
                       names=['x', 'y'])

plt.xlim(-0.2 * domain["x"], domain['x'] * 1.2)
plt.ylim(-0.2 * domain["y"], domain['y'] * 1.2)
plt.xlabel('x')
plt.ylabel('y')
# a_circle = plt.Circle((.3, .4), .005)
# plt.gca().add_artist(a_circle)
b_circle = plt.Circle((.3, .3), .3, fill=False)
plt.gca().add_artist(b_circle)


plt.savefig('exports/test.png')


def update(frame_number):
    scat.set_offsets(np.c_[data[frame_number]["x"], data[frame_number]["y"]])


animation = animation.FuncAnimation(fig, update, interval=40, frames=totalFrames)
# plt.show()
animation.save('exports/im.mp4', writer=writer)
