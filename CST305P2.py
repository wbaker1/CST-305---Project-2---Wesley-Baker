#CST-305 Project 2 - Matthew Powers and Wesley Baker 
#This project is a calculation of the first-order differential of the equation given of ym/(np.exp(xm)-1) using odeint and runge kutta and then comparing them


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

# runge_kutta formula learned in class
def runge_kutta(xm, ym, h):
    # xy equation
    def equation(xm,ym):
        return ym/(np.exp(xm)-1)
    k1 = equation(xm,ym)
    k2 = equation((xm+(h/2)),(ym+(h/2)*k1))
    k3 = equation((xm+(h/2)),(ym+(h/2)*k2))
    k4 = equation((xm+h),(ym+(h*k3)))
    # what xm+h equals after calculations
    xnext = xm + h
    ynext = ym + (h/6)*(k1 +2*k2+2*k3+k4)
    # returns all values above
    return (xnext, ynext, k1, k2, k3, k4)

# number of loops
iterations = 1000
# initial x value
xm = 1
# initial y value
ym = 5
# h (step value)
h = 0.02

runge_kutta_start_time = time.time()
# start with initial x and y values
xmArray = [1]
ymArray = [5]

# loop through runge kutta 1000 times and adds to x and 7 arrays soit can be plotted
print("Total Runge Kutta Computations: " + str(iterations))
for iteration in range(1,iterations+1):
    val = runge_kutta(xm, ym, h)
    xm = val[0]
    xmArray.append(val[0])
    ym = val[1]
    ymArray.append(val[1])
runge_kutta_end_time = time.time()
runge_kutta_total_time = runge_kutta_end_time - runge_kutta_start_time
print("Compute time of 1000 values using Runge Kutta: " + str(runge_kutta_total_time))
# initial condition for odeint
y0 = 5

# function that returns dy/dt
def model(y0, t):
    dydt = y0/(np.exp(t)-1)
    return dydt


# time points
t = np.linspace(1,21,1000)

odeint_start_time = time.time()
# solve ODE for y1
y1 = odeint(model,y0,t)
odeint_end_time = time.time()
odeint_total_time = odeint_end_time - odeint_start_time
print("Compute time of 1000 values using Odeint: " + str(odeint_total_time))

# plot result of y1
odeintPlot = plt
odeintPlot.plot(t,y1,'r-', label='Odeint line')

# label neccessary info on graph
odeintPlot.title('Odeint Graph')
odeintPlot.xlabel('x val')
odeintPlot.ylabel('y val')
odeintPlot.legend()
# Display the graph
odeintPlot.show()

# plot results for runge kutta
rungeKuttaPlot = plt
rungeKuttaPlot.plot(xmArray,ymArray,'b-', label='Runge Kutta line')

# label neccessary info on graph
rungeKuttaPlot.title('Runge Kutta Graph')
rungeKuttaPlot.xlabel('x val')
rungeKuttaPlot.ylabel('y val')
rungeKuttaPlot.legend()
# Display the graph
rungeKuttaPlot.show()

# plot results for runge kutta and odeint over eachother
mixedPlot = plt
mixedPlot.plot(t,y1,'r-', label='Odeint line')
mixedPlot.plot(xmArray,ymArray,'b--', label='Runge Kutta line')


# label neccessary info on graph
mixedPlot.title('Runge Kutta and Odeint Graph')
mixedPlot.xlabel('x val')
mixedPlot.ylabel('y val')
mixedPlot.legend()
# Display the graph
mixedPlot.show()