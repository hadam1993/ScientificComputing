# ScientificComputing Project 1

The program for project 1 for the scientific
computing course. The project was to create a program
that takes in four data points and then constructs a
cubic spline that maintains C2 continuity.
The program gives the user the option to view an example
with four data points, constructs a cubic spline from
the data points, and then takes samples in between the
first and last input data points to create a plot. 

The data in the example was taken from the function 4sin(3*x)x +4

The coordinates that were used to generate the spline were the following:
(0.5,1.896322537980259) , (0.79,0.05016729632997663),
(1.62,4.795813609335379), (2.48,16.021808673227227)

A natural spline was used to interpolate data points in between the
control points that were used, meaning that the slopes of the
end points are 0.

In the plot there is the original function that was used
to start with 4 data points and there is also the
interpolated function that was created via a cubic spline.

The user of this program also has the capability to 
input their own data and view the results from a cubic
spline that is created from their data, or can be generated
randomly. The user is asked to put in a seed for reproducibility,
if they choose to have random data points given. The user
also has the ability to define the slopes of their initial and
end points.

The code for the construction of the cubic splines
can be found in the following methods:
create_banded_mat() this creates the matrix A
create_b_mat() this creates the vector b
spline_func() the function used to get outputs 
from the spline
interpolate() decides which slopes to give to
spline_func() based on where the sample was taken
