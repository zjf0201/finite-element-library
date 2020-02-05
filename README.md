# finite-element-library
This is a three-dimensional MT finite element modeling program based on model reduction
**Authors:** 
Jifeng Zhang [zjf0201@126.com]

Jiren Liu [2439260700@qq.com]


                                                             Code description:

The code computes the electromagnetic response of 3D model for  magnetotelluric method.
It uses the basic hexahedral mesh and finite element method. Model reduction was adopted to solve the large sparse system to save a lot of runtime. Compared with the finite element 
approach based on direct solver Pardiso,  Caculation speed was greatly increases and runtime is less than one tenth of the original algorithm . So it is an efficient and fast modeling approach in geophysical electromagnetic forward.  

Instructions for three-dimensional MT foward modeling program.

1. please modify the file cmdfile.txt,It includes the path of the grid files and the frequency file.
2. frequency file:Number of frequency points,TM-mode(0)/TE-mode(1)，The number of iterations
                                          ......
                            Magnitude of frequences

            grid file： Grid information is stored
       output files:   4 files
    node_num.txt:  "x:"    Number of grids in the x direction
                              "y:"    Number of grids in the y direction
                              "z:"    Number of grids in the z direction
3. Please arrange the frequencies in order.
4. The mesh boundary must be selected far enough and the mesh generation must be dense enough at the interface of electrical change.
5.  If the results are bad, try changing the grid, frequency range, and number of iterations.
6.  If you find some problems,please don't hesitate to comment.
