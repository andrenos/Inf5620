from cooling import *

T = [37.]
k = 3.7e-5
dt = 500.
T = cooling(T, k, 20., 24*3600, dt, theta = 0.5)
	
#Find time of death
for i in range(len(T))[0:-2]:
	print T[i][0]
	if(T[i][0] > 25.8 and T[i+1][0] <= 25.8):
		print("Death occured between %4.2f and %4.2f hours before 3pm" % (i*dt/3600, (i+1)*dt/3600s))
		break