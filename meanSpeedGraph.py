import math
import os.path
import matplotlib.animation as animation
import matplotlib.pyplot as plt


import numpy as np
speeds = []

standardDev =[]

speeds.insert(100,214)
speeds.insert(200,2)
speeds.insert(300,2)
speeds.insert(400,2)
speeds.insert(500,2)
speeds.insert(600,2)
speeds.insert(700,2)
speeds.insert(800,2)

standardDev.insert(100,56)
standardDev.insert(200,2)
standardDev.insert(300,2)
standardDev.insert(400,2)
standardDev.insert(500,2)
standardDev.insert(600,2)
standardDev.insert(700,2)
standardDev.insert(800,2)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(speeds)
plt.legend(['speeds'])
# plt.plot(kineticEnergyRot)


plt.savefig('exports/speeds.png')
