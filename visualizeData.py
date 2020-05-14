import os.path
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

types = ['float', 'float', 'float', 'float', 'float', 'float', 'float']
data = []

# Set up the codec for the video file
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Lucas H'), bitrate=1800)
fig = plt.figure(figsize=(7, 7))

# this counts the number of frames
path = "./cmake-build-debug/datas/"
num_files = len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])
print(num_files)
totalFrames = num_files
for i in range(0, totalFrames):
    fileName = path + "grain" + str(i) + ".txt"
    data.insert(i, np.genfromtxt(fileName,
                                 delimiter=',',
                                 dtype=types,
                                 names=["ID", 'x', 'y', 'vx', 'vy', 'theta', 'radius']))
showingFrame = 0
scat = plt.scatter(data[showingFrame]["x"], data[showingFrame]['y'], alpha=0.5, s=data[0]['radius'] * 1 * 100000)
plt.title('Scatter plot test')
# plt.gca().set_aspect('equal', adjustable='box')
plt.axis("equal")
plt.xlim(-.4, .4)
plt.ylim(-.4, .4)
plt.xlabel('x')
plt.ylabel('y')
a_circle = plt.Circle((0, 0), .3, fill=False)
plt.gca().add_artist(a_circle)
# plt.show()


def update(frame_number):
    scat.set_offsets(np.c_[data[frame_number]["x"], data[frame_number]["y"]])


animation = animation.FuncAnimation(fig, update, interval=10, frames=totalFrames)
# plt.show()
animation.save('exports/im.mp4', writer=writer)
