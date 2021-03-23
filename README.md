# Q-ball-System

## Acknowledgements

This research project was developed jointly with Natalie Ho and our supervisor Michael Kinach.

## Purpose

[`qball.py`](https://github.com/Vismai-Khanderao/Q-ball-System/blob/main/qball.py) is made as part of a research project which would allow further study of Q-ball dynamics in more complex scenarios.

## Methods

I highly recommend reading through the Q-ball scattering section of [Relativistic Scattering of Solitons in
Nonlinear Field Theory](https://laplace.phas.ubc.ca/Members/matt/Doc/Theses/Phd/gutierrez.pdf) as this is the basis of the script in this project.

Q-balls are described by Scalar Fields and hence have a set of ODEs which can model them. 

The system of Q-Balls' parameters are input into [`params.csv`](https://github.com/Vismai-Khanderao/Q-ball-System/blob/main/params.csv), these parameters are used to solve the Boundary-Eigenvalue problem for the ODE using the Shooting Method.

Using the eigenvalue and plot generated from this, an exponential decay tail is added to create the profile of a static Q-ball. The Q-ball profile is then interpolated on the (x,y) grid using Lagrange interpolation  

To take relativistic effects into account, the Scalar Fields are then Lorentz Transformed in the direction the Q-ball is travelling depending on the velocity and direction given in the input file.

The effects of individual Q-balls are imposed creating a system at t=0 consisting of 4 Scalar Fields. 

![Graphs](https://user-images.githubusercontent.com/59114226/112082872-50d2e400-8b43-11eb-9643-807dee077234.png)

A key aspect of the script was to provide readability which is the reason for verbose code in certain areas.

## Outcome

The initial system data generated can then be passed onto a separate existing code which computes the time evolution and generates the simulation you see below.

![Triple Collision](https://user-images.githubusercontent.com/59114226/112082760-2123dc00-8b43-11eb-8fcf-c84b54b8d01f.mp4)

The Python script can now allow further study on new and different types of collisions such as the 3 Q-ball collision.

This was presented at the **UBC Multidisciplinary Undergraduate Research Conference 2021** under the title *The Complex Dynamical Behaviour of Q-balls in Collisions*




