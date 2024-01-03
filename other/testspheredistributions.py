import matplotlib.pyplot as plt
import numpy as np

fig=plt.figure()
ax = fig.add_subplot(projection="3d")

n = 10000

# biased around poles
# theta1 = np.random.uniform(0, 2*np.pi, n)
# theta2 = np.random.uniform(0, 2*np.pi, n)
# r = np.random.uniform(size=n)**(1/3)

# x = r * np.sin(theta1) * np.cos(theta2)
# y = r * np.sin(theta1) * np.sin(theta2)
# z = r * np.cos(theta1)

# uniform
r = np.random.uniform(size=n)**(1/3)
x = np.random.normal(size=n)
y = np.random.normal(size=n)
z = np.random.normal(size=n)
abs = np.sqrt(x**2 + y**2 + z**2)
x = x * r / abs
y = y * r / abs
z = z * r / abs

ax.scatter(x, y, z, s=1, alpha=0.1)
plt.show()