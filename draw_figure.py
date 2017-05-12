import matplotlib.pyplot as plt
import numpy as np

fp = open('/home/zt/Maillage/sibson/errors.txt', 'r')
errors = fp.readline().split(' ')[:-1]
fp.close()
errors = np.array(errors).reshape((-1, 3)).T

plt.figure(figsize=(15, 8))
plt.plot(errors[0, :-10], color = 'red')
plt.plot(errors[1, :-10], color = 'green')
plt.plot(errors[2, :-10], color = 'blue')
plt.xlabel('iteration')
plt.ylabel('mae')
plt.grid()
plt.title('Training Error')
plt.show()
