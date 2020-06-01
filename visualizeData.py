import os.path
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

types = ['float', 'float', 'float', 'float', 'float', 'float', 'float']
domainTypes = ['float', 'float']
data = []

# Set up the codec for the video file
writer = animation.FFMpegWriter(fps=25, metadata=dict(artist='Lucas H'), bitrate=1800)
fig = plt.figure(figsize=(7, 7))

# this counts the number of frames
path = os.path.join("out","x64-Release","datas")
num_files = len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])
print(num_files)
totalFrames = num_files - 1
for i in range(0, totalFrames):
    fileName = path + "grain" + str(i) + ".txt"
    data.insert(i, np.genfromtxt(fileName,
                                 delimiter=',',
                                 dtype=types,
                                 names=["ID", 'x', 'y', 'vx', 'vy', 'theta', 'radius']))
showingFrame = 1
norm = plt.Normalize(0, 20)
# cmap = plt.get_cmap('jet')
# speed =
velocity = np.square(data[showingFrame]["vy"]) + np.square(data[showingFrame]["vx"]) + np.divide(
    np.square(data[showingFrame]["theta"]), 100)
scat = plt.scatter(data[showingFrame]["x"], data[showingFrame]['y'], alpha=0.5,
                   s=pow(data[0]['radius'], 2) * 60 * 10000)
plt.title('Scatter plot test')
# plt.gca().set_aspect('equal', adjustable='box')
plt.axis("equal")
plt.colorbar()
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
b_circle = plt.Circle((.33, .33), .3, fill=False)
plt.gca().add_artist(b_circle)

plt.savefig('exports/test.png')


def update(frame_number):
    scat.set_offsets(np.c_[data[frame_number]["x"], data[frame_number]["y"]])
    # velocityUpdated = np.square(data[frame_number]["vy"]) + np.square(data[frame_number]["vx"]) + np.divide(np.square(
    #     data[frame_number]["theta"]), 100)
    # scat.set_array(velocityUpdated)


animation = animation.FuncAnimation(fig, update, interval=40, frames=totalFrames)
# plt.show()
animation.save('exports/im.mp4', writer=writer)
