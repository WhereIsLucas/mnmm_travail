import os.path
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

types = ['float', 'float', 'float', 'float', 'float', 'float', 'float']
ballTypes = ['float', 'float', 'float', 'float', 'float', 'float']
domainTypes = ['float', 'float']
infosTypes = ['string', 'string', 'string']
lineTypes = ['float', 'float']
ball_data = []
dataInfos = []
dataGrain = []
dataPlanes = []

# Set up the codec for the video file
Writer = animation.writers['ffmpeg']
writer = Writer(fps=25, metadata=dict(artist='Boris'), bitrate=1800)
fig = plt.figure(figsize=(7, 14))

# this counts the number of frames
base_path = "../cmake-build-debug/data/"
path_grains = "../cmake-build-debug/data/grains/"
path_ball = "../cmake-build-debug/data/ball/"
path_plane = "../cmake-build-debug/data/plane/"

# Ball
num_files = len([f for f in os.listdir(path_ball) if os.path.isfile(os.path.join(path_ball, f))])
print(num_files)
totalFrames = num_files
for i in range(0, totalFrames):
    fileName = path_ball + "ball" + str(i) + ".txt"
    ball_data.insert(i, np.genfromtxt(fileName,
                                 delimiter=',',
                                 dtype=ballTypes,
                                 names=['x', 'y', 'vx', 'vy', 'theta', 'radius']))
showingFrame = 0
scat = plt.scatter(ball_data[showingFrame]["x"], ball_data[showingFrame]['y'], alpha=0.5, s=ball_data[0]['radius'] * 1.3 * 1000000,
                   facecolors="none", edgecolors="red")

# Grain
num_files = len([f for f in os.listdir(path_grains) if os.path.isfile(os.path.join(path_grains, f))])
print(num_files)
totalFramesGrains = num_files
if totalFramesGrains > 0:
    for i in range(0, totalFrames):
        fileName = path_grains + "grain" + str(i) + ".txt"
        dataGrain.insert(i, np.genfromtxt(fileName,
                                          delimiter=',',
                                          dtype=types,
                                          names=['ID', 'x', 'y', 'vx', 'vy', 'theta', 'radius']))
    showingFrame = 0
    scatGrain = plt.scatter(dataGrain[showingFrame]["x"], dataGrain[showingFrame]['y'], alpha=0.5,
                            s=dataGrain[0]['radius'] * 1 * 100, color="blue")

plt.title('Scatter plot test')
plt.axis("equal")
domain = np.genfromtxt(base_path + "domain.txt",
                       delimiter=',',
                       dtype=domainTypes,
                       names=['x', 'y'])

plt.xlim(-1 * domain["x"], domain['x'] * 1.2)
plt.ylim(-0.2 * domain["y"], domain['y'] * 1.2)
plt.xlabel('x')
plt.ylabel('y')
# a_circle = plt.Circle((.3, .3), 1.)
# plt.gca().add_artist(a_circle)

# Plane
num_files = len([f for f in os.listdir(path_plane) if os.path.isfile(os.path.join(path_plane, f))])
print(num_files)
totalFrames = num_files
for i in range(0, totalFrames):
    fileName = path_plane + "plane" + str(i) + ".txt"
    dataPlanes.insert(i, np.genfromtxt(fileName,
                                 delimiter=',',
                                 dtype=lineTypes,
                                 names=['m', 'p']))
showingFrame = 0
x = np.linspace(-1 * domain["x"], domain['x'] * 2, 100)
plane_line, = plt.plot(x, (dataPlanes[showingFrame]["m"] * x) + dataPlanes[showingFrame]["p"], color="black")


# plt.show()
plt.savefig("exports/im.png")


#
def update(frame_number):
    x = np.linspace(-1 * domain["x"]*2, domain['x'] * 2, 100)
    scat.set_offsets(np.c_[ball_data[frame_number]["x"], ball_data[frame_number]["y"]])
    if totalFramesGrains > 0:
        scatGrain.set_offsets(np.c_[dataGrain[frame_number]["x"], dataGrain[frame_number]["y"]])
    if totalFrames > 0:
        plane_line.set_ydata((dataPlanes[frame_number]["m"] * x) + dataPlanes[frame_number]["p"])

animation = animation.FuncAnimation(fig, update, interval=40, frames=totalFrames)
animation.save('exports/im.mp4', writer=writer)
