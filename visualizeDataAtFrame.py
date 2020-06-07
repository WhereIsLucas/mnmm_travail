import os.path
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

types = ['float', 'float', 'float', 'float', 'float', 'float', 'float']
domainTypes = ['float', 'float']

# Set up the codec for the video file
# writer = animation.FFMpegWriter(fps=25, metadata=dict(artist='Lucas H'), bitrate=1800)
fig = plt.figure(figsize=(7, 7))

# this counts the number of frames
for k in range(1, 9) :
    base_path = "./datas/"+ str(k) +"00 pi/"
    path = base_path+"datas/"
    num_files = len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])
    print(num_files)
    data = []
    totalFrames = num_files - 1
    for i in range(0, totalFrames):
        fileName = os.path.join(path,"grain" + str(i) + ".txt")
        data.insert(i, np.genfromtxt(fileName,
                                     delimiter=',',
                                     dtype=types,
                                     names=["ID", 'x', 'y', 'vx', 'vy', 'theta', 'radius']))

    for a in range(0,8):
        showingFrame = a* 4
        # cmap = plt.get_cmap('jet')
        # speed =
        velocity = np.square(data[showingFrame]["vy"]) + np.square(data[showingFrame]["vx"]) + np.divide(
            np.square(data[showingFrame]["theta"]), 100)
        scat = plt.scatter(data[showingFrame]["x"], data[showingFrame]['y'], alpha=0.5,
                           s=pow(data[0]['radius'], 2) * 60 * 10000)
        plt.title('Frame '+str(showingFrame)+' - cohesion parameter = '+str(k)+'00')
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
        b_circle = plt.Circle((.33, .33), .3, fill=False)
        plt.gca().add_artist(b_circle)

        plt.savefig('exports/'+str(k)+'00_screenshot'+str(a)+'.png')
        plt.clf()
        plt.cla()
