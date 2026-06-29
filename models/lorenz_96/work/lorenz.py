#!/bin/python

import numpy as np
import matplotlib.pyplot as plt

#Compile and run lorenz 96 on DART too !
#This python file intends to help generate and analyse the lorenz_96 model, the mathematical set of equations created by Edward Lorenz in 1996

def dxdt(x,F):
    """
    This will work if F is a number or a numpy array of equal length to x
    """
    return (np.roll(x,-1) - np.roll(x,2))*np.roll(x,1) - x + F

def vanilla(x="Not Given",F=8,size=100,cycles=100,step=0.05/6):
    #Check if Euler VS RK4 is better -> but in general RK4 is more standard and you get less arguments
    #Also, a step of 0.05 is usually more stable with RK4, otherwise do 0.05/6
    if (type(x) == str) and (x == "Not Given"):
        x = np.ones(size)
    data = [x]
    time = [0]
    for i in range(int(cycles/step)):
        data.append(x + dxdt(x,F)*step)
        time.append(time[-1] + step)
        x = data[-1]
    return [time,data]

def RungeKutta4(forward_function,current_value,current_time,dt):
    order_1 = forward_function(current_value,current_time)
    order_2 = forward_function(current_value + (order_1 * dt / 2), current_time + (dt / 2))
    order_3 = forward_function(current_value + (order_2 * dt / 2), current_time + (dt / 2))
    order_4 = forward_function(current_value + order_3 * dt, current_time + dt)
    return current_value + dt / 6 * (order_1 + (2 * order_2) + (2 * order_3) + order_4)

def vanilla_RK4(x="Not Given",F=8,size=100,cycles=100,step=0.05):
    if (type(x) == str) and (x == "Not Given"):
        x = np.ones(size)
    data = [x]
    time = [0]
    def forward_function(current_value,current_time):
        return dxdt(current_value,F)
    for i in range(int(cycles/step)):
        data.append(RungeKutta4(forward_function,x,time[-1],step))
        time.append(time[-1] + step)
        x = data[-1]
    return [time,data]

def dependent_forcing(F,x=None,size=100,cycles=100,step=1/100):
    """
    F should be a function: F(values,time) --> float/int or numpy array
    """
    if (type(x) == str) and (x == "Not Given"):
        x = np.ones(size)
    data = [x]
    time = [0]
    for i in range(cycles):
        data.append(x + dxdt(x,F(data[-1],time[-1]))*step)
        time.append(time[-1] + step)
        x = data[-1]
    return [time,data]

def quickplot(datalist,title,filename=None):
    """
    Assumes that the title of the graph will be sanitary input for the filename
    """
    filename = filename or title + ".png"
    n_graphs = len(datalist[1][0])
    graphs = np.array(datalist[1])
    n_rows = int(np.ceil(n_graphs/3))
    fig,ax = plt.subplots(nrows=n_rows,ncols=3,figsize=(15,4.8*n_rows))
    fig.suptitle(title)
    for i in range(n_rows):
        for j in range(3):
            if i*3 + j < n_graphs:
                ax[i,j].plot(datalist[0],graphs[:,i*3+j])
                ax[i,j].set_title("Variable index {}".format(i*3+j))
    plt.savefig(filename)
    plt.close()
    return

def run_simple_ensemble(N_experiments, sided = False, x = "Not Given", F = 8, size = 360, cycles = 1000, step = 1/100, style = "Euler"):
    """
    To make life easier, let us assume that only normally distributed differences occur
    To make life even easier, let us assume that only the first index gets perturbed
    """
    RNGesus = np.random.default_rng()
    perturbations = []
    if (type(sided) == bool) and (sided == False):
        for i in range(N_experiments):
            perturbations.append(RNGesus.random() - 0.5)
    elif sided == "right":
        for i in range(N_experiments):
            perturbations.append(RNGesus.random())
    elif sided == "left":
        for i in range(N_experiments):
            perturbations.append(RNGesus.random() - 1)
    data = []
    time = []
    if (type(x) == str) and (x == "Not Given"):
        x = np.ones(size)
    elif isinstance(x, list):
        x = np.array(x)
    elif isinstance(x, np.ndarray) != True:
        print("Unrecognised Data Type for x-vector")
        return []
    for i in range(N_experiments):
        x_i = x.copy()
        x_i[0] += perturbations[i]
        if style == "Euler":
            experiment_data = vanilla(x=x_i,F=F,size=size,cycles=cycles,step=step)
        elif style == "RK4":
            experiment_data = vanilla_RK4(x=x_i,F=F,size=size,cycles=cycles,step=step)
        data.append(experiment_data[1])
        if i == 0:
            time = experiment_data[0]
    return [time,np.array(data)]

