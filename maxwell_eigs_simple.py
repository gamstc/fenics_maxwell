"""
First trial in FEniCS to implement the Maxwell eigenvalue problem in 2D/3D using Nedelec elements
(here compared with the Laplace eigenvalue problem using Lagrange elements)

-Nabla x Nabla x u = k^2 u   on the domain
			 n x u = 0 	     on the boundary
     
"""

from dolfin import *


#----------- preliminary settings and checks -----------------------------------------------------------------------------------------------------------------


# Test for PETSc and SLEPc
if not has_linear_algebra_backend("PETSc"):
    print "DOLFIN has not been configured with PETSc. Exiting."
    exit()

if not has_slepc():
    print "DOLFIN has not been configured with SLEPc. Exiting."
    exit()


# Get user input about which problem shall be solved
prob = raw_input("Please specify if you want to solve the Laplace or the Maxwell eigenproblem. ")
prob = prob.capitalize()
allowed_problems = ["Maxwell", "Laplace"]
if prob in allowed_problems:
    print "The ", prob, " eigenvalue problem will be solved."
else:
    print "That was not a valid choice."
    exit()
if prob == "Maxwell":
    ind = 1
else: 
    ind = 0
    
    
# Get user input about how many eigenvalues to compute
n_eigs = raw_input("Please specify the number of eigenvalues you want computed. ")
n = int(n_eigs)
if n > 100:
    print "Sorry, that are too many."
    exit()
elif n < 1:
    print "Sorry, that are too few."
    exit()
else:
    print "The first ", n, " eigenvalues and eigenfunctions will be computed."


# Get user input about how many grid points to use
grid_p = raw_input("Please specify the number of grid points to be used. ")
p = int(grid_p)
if p > 100:
    print "Sorry, that are too many."
    exit()
elif p < 1:
    print "Sorry, that are too few."
    exit()
    
    
#----------- actual calculations -----------------------------------------------------------------------------------------------------------------


# Define mesh, function space
#mesh = UnitCubeMesh(p,p,p)
mesh = UnitSquareMesh(p,p)
#mesh = CircleMesh(Point(0,0,0), 1, 0.05)

# Maxwell or Laplace function space
if ind > 0:
    V = FunctionSpace(mesh, "N1curl", 1)
else:
    V = FunctionSpace(mesh, "Lagrange", 1)


# Define basis and bilinear form
u = TrialFunction(V)
v = TestFunction(V)

# Maxwell or Laplace bilinear form
if ind > 0:
    a = dot(curl(u), curl(v))*dx
else:
    a = dot(grad(u), grad(v))*dx


# Assemble stiffness form
A = PETScMatrix()
assemble(a, tensor=A)

# Create eigensolver
eigensolver = SLEPcEigenSolver(A)

# Compute n eigenvalues of A x = \lambda x
print "Computing the first", n, "eigenvalues. This can take a minute."
eigensolver.solve(n)


#----------- vizualization -----------------------------------------------------------------------------------------------------------------


for i in range(0,n): 
    # Extract ith eigenpair
    r, c, rx, cx = eigensolver.get_eigenpair(i)

    print i, ". eigenvalue: ", r

    # Initialize function and assign eigenvector
    u = Function(V)
    u.vector()[:] = rx

    # Plot eigenfunction
    plot(u)

interactive()
 
