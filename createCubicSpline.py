#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 11:38:25 2020

@author: ahonts

This program for project 1 for the scientific
computing course. The project was to create a program
that takes in four data points and then constructs a
cubic spline that maintains C2 continuity.
The program gives the user the option to view an example
with four data points, constructs a cubic spline from
the data points, and then takes samples in between the
first and last input data points to create a plot. 
In the plot there is the original function that was used
to start with 4 data points and there is also the
interpolated function that was created via a cubic spline.

The user of this program also has the capability to 
input their own data and view the results from a cubic
spline that is created from their data.

The code for the construction of the cubic splines
can be found in the following methods:

create_banded_mat() this creates the matrix A

create_b_mat() this creates the vector b

spline_func() the function used to get outputs 
from the spline

interpolate() decides which slopes to give to
spline_func() based on where the sample was taken

"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve

def create_banded_mat(data):
    '''
    Args:
        data: consists of the (x,y) coordinates that
        are used to build the spline.
    Returns:
        returns the matrix A used to solve
        for the unknown slopes of the splines that
        will be used for data interpolation.
        A is a banded matrix with a bandwith of 2
        where the diagonal of the matrix is constructed
        by 2*(delta_x_i + delta_x_i+1) where i ranges
        from 1 to the number of data points.
        The bottom band of the diagonal consists of
        delta_x_i where i ranges from 3 to the number of 
        data points. The upperband of the diagonal
        consists of delta_x_i where i ranges from
        1 to the number of data points minus 2
    '''
    A = np.zeros((len(data)-2,len(data)-2))
    for i in range(len(data)-2):
        A[i][i] = 2*((data[i+1][0] - data[i][0]) +\
                     (data[i+2][0] - data[i+1][0]))
    for i in range(len(data)-3):
        A[i][i+1] = (data[i+1][0] - data[i][0])
    for i in range(len(data)-3):
        A[i+1][i] = data[i+3][0] - data[i+2][0]
    return A

def create_b_mat(data,s):
    '''
    Args:
        data: data consists of the (x,y) coordinates that
        are used to build the spline.
        s: s is an array that contains the inital and
        end slopes.
    Returns:
        returns the vector b that is used to solve
        for the unknown slopes of the splines used
        in data interpolation. b is constructed by
        3*(delta_x_i+1 * y_i_prime) where i ranges from
        1 to the number of data points - 1.
        If the initial and end slopes are not 0,
        then the first element in the vector has
        delta_x_2 * initial slope subtracted from it
        and the last element in the vectior has
        delta_x_end-1 * end slope subtracted from it
    '''
    b = np.zeros((len(data)-2,1))
    for i in range(b.shape[0]):
        delta_x_2 = data[i+2][0] - data[i+1][0]
        delta_x_1 = data[i+1][0] - data[i][0]
        y_prime_1 = (data[i+1][1] - data[i][1])/delta_x_1
        y_prime_2 = (data[i+2][1] - data[i+1][1])/delta_x_2
        b[i][0] = 3*(delta_x_2*y_prime_1 +\
                       delta_x_1*y_prime_2)
        if i == 0 or i == b.shape[0] - 1:
            b[i][0] -= delta_x_2*s[i]
    return b

def spline_func(x_in,x1,x2,y1,y2,s1,s2):
    '''
    Args:
        x_in: is the input we wish to interpolate
        an output for. This input falls in between
        the given inputs x1 and x2
        x1: is the left given data point that x_in is
        between
        x2: is the right given data point that x_in is
        between
        y1: is the output of the left given data 
        point that x_in is between
        y2: is the output of the right given data 
        point that x_in is between
        s1: is slope of the spline function at x1
        s2: is slope of the spline function at x2
    Returns:
        returns: the interpolation value for the 
        spline function at the input of x_in.
    '''
    delta_x = (x2 - x1)
    y_prime = (y2 - y1)/delta_x
    y_double_prime = (y_prime - s1)/(delta_x)
    y_triple_prime = (s1 - 2*(y_prime) + s2)/(delta_x**2)
    return (y1 + s1*(x_in - x1) +\
           y_double_prime*(x_in - x1)**2 +\
           y_triple_prime*((x_in - x1)**2)*(x_in-x2))

def interpolate(data,samples,s):
    '''
    Args:
        data: consists of the (x,y) coordinates that
        are used to build the spline.
        samples: are the x coordinates we wish to
        have and interpolated output for
        s: is a vector of slopes where s_i is the
        slope of the spline at x_i where i ranges
        from 1 to number of data points
    Returns:
        returns: a list of outputs corresponding to
        each of the the samples that were put 
        into this function
    '''
    outputs = []
    for sample in samples:
        if sample < data[0][0] or sample > data[len(data)-1][0]:
            return None
        if sample == data[len(data)-1][0]:
            i = len(data)-1
            outputs.append(
                    spline_func(sample,data[i][0],
                                data[i+1][0],data[i][1],
                                data[i+1][1],s[i],s[i+1]))
        else:
            for i in range(len(data)-1):
                if sample >= data[i][0] and sample < data[i+1][0]:
                    outputs.append(
                        spline_func(sample,data[i][0],
                                    data[i+1][0],data[i][1],
                                    data[i+1][1],s[i],s[i+1]))
                    break
        
    return outputs

def is_number(n):
    '''
    Args:
        n: a string
    Returns:
        returns True if n is a number or False
        if n is not a number
    '''
    try:
        float(n)
    except ValueError:
        return False
    return True

def get_input(xs):
    '''
    Args:
        xs: is a list of inputs that have been
        given to the program
    Returns:
        False if the user did not enter a number
        or if the input is already used. Otherwise
        it will return the number as a float
    '''
    x = input("please enter an input number: \n")
    if is_number(x) == False:
        print("you did not enter a number.\n")
        return False
    elif float(x) in xs:
        print("You have already used the input %s."%x)
        return False
    return float(x)

def get_output(x):
    '''
    Args:
        x: is the input for the current output that
        we want from the user
    Returns:
        False if the user did not enter a number.
        Otherwise it will return the output as a float
    '''
    y = input("please enter an output number for the input %f: \n"%x)
    if is_number(y) == False:
        print("you did not enter a number.\n")
        return False
    return float(y)

def get_num_data_pts():
    '''
    Args:
        None
    Returns:
        False if the user did not enter a number or if
        the user did not enter at least 3. We need
        at least 3 data points to be able to construct
        the spline. Otherwise it will return the 
        num_data_pts as an int
    '''
    num_data_pts = input("please enter how many data"+ 
                         "points you will enter: \n")
    if is_number(num_data_pts) == False:
        print('You did not enter a valid number \n')
        return False
    if int(num_data_pts) < 3:
        print('Please enter a number larger than 3\n')
        return False
    return(int(num_data_pts))

def get_init_slope():
    '''
    Args:
        None
    Returns:
        False if the user did not enter a number.
        Otherwise it will return the slope given
        by the user as a float
    '''
    s_1 = input("Please enter the slope of the inital point: \n")
    if is_number(s_1) == False:
        print('You did not enter a valid number \n')
        return False
    return(float(s_1))

def get_end_slope():
    '''
    Args:
        None
    Returns:
        False if the user did not enter a number.
        Otherwise it will return the slope given
        by the user as a float
    '''
    s_end = input("Please enter the slope of the end point: \n")
    if is_number(s_end) == False:
        print('You did not enter a valid number \n')
        return False
    return(float(s_end))

def func(t):
    return(5*np.sin(2*t+np.pi)*t + 4)

def four_point_example():
    '''
    This function will display the results of the 
    cubic spline interpolation from the 4 data points
    (0.5,1.896322537980259) , (0.79,0.05016729632997663),
    (1.62,4.795813609335379), (2.48,16.021808673227227)
    
    These data points were collected from the function:
        4sin(3*x)x +4
    '''
    print('\nThis example comes from the four data points\n'+\
          '(0.5,1.896322537980259)\n(0.79,0.05016729632997663)\n'+\
          '(1.62,4.795813609335379)\n(2.48,16.021808673227227)\n\n'+\
          'These data points come from the function\n\n'+\
          '4sin(3*x)x +4\n\nThe slopes of the cubic spline '+\
          'at the end points are both 0\n'+\
          'The graph shows the interpolated '+\
          'function plotted against the true function\n'+\
          'The red stars are the four data points')
    data = [(0.5,1.896322537980259),
            (0.79,0.05016729632997663),
            (1.62,4.795813609335379),
            (2.48,16.021808673227227)]
    init_slope = 0.0
    end_slope = 0.0
    s = np.zeros((len(data),1))
    #Create the matix A to solve for the unknown slopes
    A = create_banded_mat(data)
    #Create the vector b to solve for the unknown slopes
    b = create_b_mat(data,s)
    #Solve for the unknown slopes using numpy's 
    #Linear algebra solver
    s_2 = solve(A,b)
    #Fill the slopes into the s vector
    s[1:len(data)-1][:] = s_2
    #Set the initial slope
    s[0][:] = init_slope
    #Set the end slope
    s[len(s)-1][:] = end_slope
    #Get Samples ranging from the smallest input value
    #to the highest input value
    samples = np.arange(data[0][0],data[len(data)-1][0],.01)
    #Interpolate the sample inputs and store the 
    #interpolated outputs in y_hat
    y_hat = interpolate(data,samples,s)
    
    y_actual = [func(i) for i in samples]
    
    #plot the resulting interpolation
    fig, ax = plt.subplots()
    stars, = ax.plot([x[0] for x in data],[x[1] for x in data],'r*',label='stars')
    interpolated, = ax.plot(samples,y_hat,'g--', label='Line 2')
    actual, = ax.plot(samples,y_actual,'b-', label='Line 1')
    ax.set_xlabel('input axis')
    ax.set_ylabel('output axis')
    plt.legend([stars,interpolated, actual], ['Given Data Points','Interpolated Function', 'Actual Function'])
    plt.show()

def run_example_or_user():
    example = input("input y to see example with 4"+\
                    " data points\ninput n to create"+\
                    " your own data:\n")
    return example


def __main__():
    example = run_example_or_user()
    while example != 'y' and example != 'n':
        example = run_example_or_user()
    if example == 'y':
        four_point_example()
    else:
        #List of inputs
        xs = []
        # List of tuples that represent (x,y) coordinates
        data = []
        #Ask the user for the number of data points
        #that they will enter
        num_data_pts = get_num_data_pts()
        #Ensure that the user entered a valid number
        while num_data_pts == False:
            num_data_pts = get_num_data_pts()
        #Ask the user for an initaial slope
        init_slope = get_init_slope()
        #Ensure that the user entered a valid number
        #0.0 is equivalent to False in python
        #init_slope != 0.0 allows the user to enter that number
        while init_slope == False and init_slope != 0.0:
            init_slope = get_init_slope()
        #Ask the user for an initaial slope
        end_slope = get_end_slope()
        #Ensure that the user entered a valid number
        while end_slope == False and end_slope != 0.0:
            end_slope = get_end_slope()
        #Get all data points from the user
        while len(data) < num_data_pts:
            #Get an input value
            x = get_input(xs)
            while x == False and x != 0.0:
                x = get_input(xs)
            #Get the corresponding output value
            y = get_output(x)
            while y == False and y != 0.0:
                y = get_output(x)
            #add x to the list of used inputs
            xs.append(x)
            #add the coordiante (x,y) to the data list
            data.append((x,y))
            #Ensure the data is sorted by input values
            #in ascending order
            data = sorted(data, key=lambda x: x[0])
        #Create the vector s for the slopes of the splines
        s = np.zeros((len(data),1))
        #Create the Matrix A to solve for the unknown slopes
        A = create_banded_mat(data)
        #Create the vector b to solve for the unknown slopes
        b = create_b_mat(data,s)
        #Solve for the unknown slopes using numpy's 
        #Linear algebra solver
        s_2 = solve(A,b)
        #Fill the slopes into the s vector
        s[1:len(data)-1][:] = s_2
        #Set the initial slope
        s[0][:] = init_slope
        #Set the end slope
        s[len(s)-1][:] = end_slope
        #Get Samples ranging from the smallest input value
        #to the highest input value
        samples = np.arange(data[0][0],data[len(data)-1][0],.01)
        #Interpolate the sample inputs and store the 
        #interpolated outputs in y_hat
        y_hat = interpolate(data,samples,s)
        #plot the resulting interpolation
        fig, ax = plt.subplots()
        stars, = ax.plot([x[0] for x in data],[x[1] for x in data],'r*')
        interpolated, = ax.plot(samples,y_hat,'g--',label='interpolated')
        ax.set_xlabel('input axis')
        ax.set_ylabel('output axis')
        plt.legend([stars,interpolated], ['Given Data Points','Interpolated Function'])
        plt.show()

__main__()