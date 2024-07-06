import numpy
import math
import scipy.linalg
from sklearn import linear_model
from sklearn import preprocessing


# Read data

A = numpy.genfromtxt("A.txt", dtype='float')
b = numpy.genfromtxt("b.txt", dtype='float')


# Do the fit

reg   = linear_model.Lasso(1.0E-9,fit_intercept=True,max_iter=100000)
#reg   = linear_model.Lasso()
reg.fit(A,b)
x     = reg.coef_


######## SVD

eps = 1.0e-5

# Make the scipy call

U,D,VT = scipy.linalg.svd(A,overwrite_a=False)
Dmat   = numpy.array((numpy.transpose(A)))

# Process output

dmax = 0.0

for i in range(0,len(Dmat)):
    if ( abs(D[i]) > dmax ) :
        dmax = abs(D[i])

    for j in range(0,len(Dmat[i])):
        Dmat[i][j]=0.0

# Cut off singular values based on fraction of maximum value as per numerical recipes.

eps_in = eps
eps = eps* dmax
nvars = 0

for i in range(0,len(D)):
    if abs(D[i]) > eps:
        Dmat[i][i]=1.0/D[i]
        nvars += 1

print ("! eps (= args.eps*dmax)          =  %11.4e" % eps_in)
print ("! SVD regularization factor      = %11.4e" % eps)

x = numpy.dot(numpy.transpose(VT),Dmat)
x = numpy.dot(x,numpy.dot(numpy.transpose(U),b))





# Calc/output the fitted values, calc RMSE

y     = numpy.dot(A,x)
yfile = open("force.txt", "w")

z = 0

for a in range(0,len(b)):
    z = z + (y[a] - b[a]) ** 2.0
    yfile.write("%13.6e\n"% y[a])
yfile.close()

# Print the parameters

numpy.savetxt('params.txt', x)

print ("RMSE: ", math.sqrt(z/float(b.shape[0])))


exit(0)
