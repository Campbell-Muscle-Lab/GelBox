import matplotlib.pyplot as plt
import numpy as np
import pywt

from skimage import (
    data, restoration, util
)
x = pywt.data.ecg()
background = restoration.rolling_ball(x, radius=80)
background2 = restoration.rolling_ball(x, radius=10)
plt.figure()
plt.plot(x, label='original')
utku2 = x - background2
utku1 = x - background

# plt.plot(x - background, label='radius=80')
# plt.plot(x - background2, label='radius=10')
plt.legend()
plt.show()