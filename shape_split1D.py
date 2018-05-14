
# coding: utf-8

# In[2]:

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint
import math

#get_ipython().magic('matplotlib inline')


# In[ ]:




# In[3]:

def Getx(nx, xmin, xmax):
    dx = np.abs(xmax - xmin)/float(nx)
    x  = xmin + np.arange(nx+1)*dx
    return x


# In[4]:

# reduce is not recognized by python 3.5.2
#def fact(n):return reduce(lambda x,y:x*y,[1]+range(1,n+1))

def fact(n):return math.factorial(n)

# In[5]:

''' Calculate binomial coefficient nCk = n! / (k! (n-k)!)
'''


def binomial(n, k):
    
    if n<0 or k<0 or k>n:
        binom=0
    else:
        binom = fact(n) // fact(k) // fact(n - k)
    
    return binom


#Print Pascal's triangle to test binomial()
def pascal(m):
    for x in range(m + 1):
        print([binomial(x, y) for y in range(x + 1)])



# In[6]:

pascal(5)


# In[7]:

def trinomial(n,k):
    
    coef = 0
    
    j_l = np.arange(n+1)
    #print("j_l",j_l)
       
    for j in j_l:
        val1 = binomial(n, j)
        val2 = binomial(2*n-2*j, n-k-j)
        coef = coef + (-1)**j * val1 * val2
        
        #print("n = %d, k = %d, j = %d, coef = %d"%(n, k, j, coef))
        #print("binomial(n, j)           = %d"%val1)
        #print("binomial(2*n-2*j, n-k-j) = %d"%val2)
        
    return coef

#Print Trinomial's triangle to test trinomial()
def tri_tri(m):
    for x in range(m + 1):
        #print("Row %d"%x)
        #print(range(-x, x+1))
        print([trinomial(x, y) for y in range(-x, x+1)])


# In[8]:

tri_tri(10)


# In[ ]:




# In[ ]:




# In[ ]:




# In[9]:

def GetOrd1(x,x0,dx):
    
    xi=(x-x0)/dx

    izeros = np.argwhere(abs(xi) > 1.)
    ileft  = np.argwhere(xi < 0.)
    iright = np.argwhere(xi >= 0.)
    
    shape = np.zeros(len(x))
    
    for ik in ileft:
        shape[ik] = xi[ik]+1.
    
    for ik in iright:
        shape[ik] = 1.-xi[ik]
        
    for ik in izeros:
        shape[ik] = 0.
    
    return (shape)


# In[10]:

def GetOrd2(x,x0,dx):
    
    xi=(x-x0)/dx

    iequ1  = np.argwhere( abs(xi) < 0.5  )
    iequ2  = np.argwhere( abs(xi) >= 0.5 ) 
    izeros = np.argwhere( abs(xi) > 1.5  )
    
    shape = np.zeros(len(x))
    
    for ik in iequ1:
        shape[ik] = 0.75 - abs(xi[ik])**2.
    
    for ik in iequ2:
        shape[ik] = 0.5*(1.5 - abs(xi[ik]))**2.
        
    for ik in izeros:
        shape[ik] = 0.
    
    return (shape)


# In[11]:

def Spline2(x,x0,dx):
    
    xi=(x-x0)/dx
   
    val = 0.
    if not isinstance(x, float) : 
        print("Cette fonction renvoie un scalaire !")
        return(val)
    
    if abs(xi)<0.5:
        val = 0.75 - abs(xi)**2.
    elif abs(xi)>=0.5 and abs(xi)<1.5:
        val = 0.5*(1.5 - abs(xi))**2.
    else:
        val = 0.
    
    return (val)


# In[12]:

x0=3.
Spline2(x0+0.5, x0+0., 0.5)


# In[13]:

def GetOrd3(x,x0,dx):
    
    xi=(x-x0)/dx

    iequ1  = np.argwhere( abs(xi) < 1.  )
    iequ2  = np.argwhere( abs(xi) >= 1. ) 
    izeros = np.argwhere( abs(xi) > 2.  )
    
    shape = np.zeros(len(x))
    
    for ik in iequ1:
        shape[ik] = 0.5*abs(xi[ik])**3. - xi[ik]**2. + 2./3.
    
    for ik in iequ2:
        shape[ik] = 4./3.*(1. - 0.5*abs(xi[ik]))**3.
        
    for ik in izeros:
        shape[ik] = 0.
    
    return (shape)


# In[14]:

def GetOrd4(x,x0,dx):
    
    xi=(x-x0)/dx

    iequ1  = np.argwhere( abs(xi) <  0.5 )
    iequ2  = np.argwhere( abs(xi) >= 0.5 ) 
    iequ3  = np.argwhere( abs(xi) >= 1.5 ) 
    izeros = np.argwhere( abs(xi) >  2.5 )
    
    shape = np.zeros(len(x))
    
    for ik in iequ1:
        shape[ik] = 1./192.*( 115. - 120.*xi[ik]**2. + 48.*xi[ik]**4. )
    
    for ik in iequ2:
        shape[ik] = 1./96.*( 55. + 20.*abs(xi[ik]) - 120.*xi[ik]**2. + 80.*abs(xi[ik])**3. - 16.*xi[ik]**4.)
        
    for ik in iequ3:
        shape[ik] = 1./24.*(2.5 - abs(xi[ik]))**4.
        
    for ik in izeros:
        shape[ik] = 0.
    
    return (shape)


# In[15]:

def getShape(order, x, x0, dx_L0):
    return {
        1: GetOrd1(x, x0, dx_L0),
        2: GetOrd2(x, x0, dx_L0),
        3: GetOrd3(x, x0, dx_L0),
        4: GetOrd4(x, x0, dx_L0),
    }.get(order, GetOrd1(x, x0, dx_L0))


# In[ ]:




# In[ ]:




# In[ ]:




# In[9]:

x = Getx(10000,10)
shape = GetOrd1(x, 4., 1.)
bmin = 2.
bmax = 6.
plt.plot(x,shape)
plt.xlim(bmin, bmax)


# In[ ]:




# In[10]:

x = Getx(10000,10)
shape = GetOrd2(x, 4., 1.)
bmin = 2.
bmax = 6.
plt.plot(x,shape)
plt.xlim(bmin, bmax)


# In[11]:

x = Getx(10000,10)
shape = GetOrd3(x, 4., 1.)
bmin = 2.
bmax = 6.
plt.plot(x,shape)
plt.xlim(bmin, bmax)


# In[12]:

x = Getx(10000,10)
shape = GetOrd4(x, 4., 1.)
bmin = 1.
bmax = 7.
plt.plot(x,shape)
plt.xlim(bmin, bmax)


# In[ ]:




# In[ ]:




# In[ ]:




# In[13]:

NbPoints = 10000
Lmax = 10.
dx = Lmax/NbPoints
print("dx = %f" %dx)

x = Getx(NbPoints,Lmax)

x0 = 4.
dx_L0 = 1.
dx_L1 = dx_L0/2.
RF = dx_L0/dx_L1

shape = GetOrd1(x, 4., 1.)
w1=0.25 * RF
shape_split = w1*GetOrd1(x, x0-dx_L1, dx_L1) + w1*GetOrd1(x, x0+dx_L1, dx_L1)             + 2.*w1*GetOrd1(x, x0, dx_L1)

# Compute the area using the composite Simpson's rule.
area_ref = scint.simps(shape, dx=dx)
area_split = scint.simps(shape_split, dx=dx)
print("area ref=", area_ref)
print("area split=", area_split)
    
bmin = 2.
bmax = 6.
plt.plot(x,shape)
plt.plot(x,shape_split,'red')
plt.xlim(bmin, bmax)


# In[ ]:




# In[17]:

NbPoints = 10000
Lmax = 10.
dx = Lmax/NbPoints
print("dx = %f" %dx)

x = Getx(NbPoints,Lmax)

x0 = 4.
dx_L0 = 1.
dx_L1 = dx_L0/2.
RF = dx_L0/dx_L1

dta = 1.5*dx_L1

w1 = 0.25
w2 = w1
w3 = 2.*0.5

print("w1 = %.2f, w2 = %.2f, w3 = %.2f"%(w1, w2, w3))

wtot = w1 + w2 + w3
print("wtot = %.2f"%wtot)

w1 = w1/wtot * RF
w2 = w2/wtot * RF
w3 = w3/wtot * RF
print("w1 = %.2f, w2 = %.2f, w3 = %.2f"%(w1, w2, w3))

shape = GetOrd2(x, x0, dx_L0)

s1 = w1*GetOrd2(x, x0-dta, dx_L1)
s2 = w2*GetOrd2(x, x0+dta, dx_L1)
s3 = w3*GetOrd2(x, x0, dx_L1)

shape_split = s1 +s2 + s3
        
# Compute the area using the composite Simpson's rule.
area_ref = scint.simps(shape, dx=dx)
area_split = scint.simps(shape_split, dx=dx)
print("area ref=", area_ref)
print("area split=", area_split)
        
bmin = 2.
bmax = 6.
plt.plot(x,shape)
plt.plot(x,shape_split,'red')
plt.plot(x,s1,'r--')
plt.plot(x,s2,'r--')
plt.plot(x,s3,'r--')
plt.xlim(bmin, bmax)


# In[ ]:




# In[ ]:




# In[35]:

#
#  APPROXIMATE TEST
#
    
NbPoints = 10000
xmin = -5.
xmax =  5.
Lmax = np.abs(xmax - xmin)
dx = Lmax/NbPoints
print("dx = %f" %dx)

x = Getx(NbPoints, xmin, xmax)

x0 = 2.
dx_L0 = 1.
dx_L1 = dx_L0/2.
RF = dx_L0/dx_L1

dta = 1.*dx_L1

shape = GetOrd2(x, x0, dx_L0)
order = 2
nbpts = (2*int(RF) - 1) + (order-1)

# just to have a taste of things
coef1 = 1.125 # 0.9
coef2 = 1.125 # 0.9
coef3 = 1.

w1 = coef1*Spline2(x0-dta, x0, dx_L0)
w2 = coef2*Spline2(x0+dta, x0, dx_L0)
w3 = coef3*Spline2(x0, x0, dx_L0)
print("w1 = %.2f, w2 = %.2f, w3 = %.2f"%(w1, w2, w3))

wtot = w1 + w2 + w3
print("wtot = %.2f"%wtot)

w1 = w1/wtot * RF
w2 = w2/wtot * RF
w3 = w3/wtot * RF
print("w1 = %.2f, w2 = %.2f, w3 = %.2f"%(w1, w2, w3))

s1 = w1*GetOrd2(x, x0-dta, dx_L1)
s2 = w2*GetOrd2(x, x0+dta, dx_L1)
s3 = w3*GetOrd2(x, x0, dx_L1)

shape_split = s1 + s2 + s3

# Compute the area using the composite Simpson's rule.
area_ref = scint.simps(shape, dx=dx)

As1 = scint.simps(s1, dx=dx)
As2 = scint.simps(s2, dx=dx)
As3 = scint.simps(s3, dx=dx)

print("As1 = %f, As2 = %f, As3 = %f"%(As1, As2, As3))

area_split = As1 + As2 +As3

print("area ref=", area_ref)
print("area split=", area_split)

#s1 = s1/area_split
#s2 = s2/area_split
#s3 = s3/area_split
#shape_split = shape_split/area_split

print("After normalization :")
area_split = scint.simps(shape_split, dx=0.001)
print("area split=", area_split)

bmin = x0 - 0.5*nbpts*dx_L1 - (order+1)*0.5*dx_L1
bmax = x0 + 0.5*nbpts*dx_L1 + (order+1)*0.5*dx_L1
plt.plot(x,shape)
plt.plot(x,shape_split,'red')
plt.plot(x,s1,'r--')
plt.plot(x,s2,'r--')
plt.plot(x,s3,'r--')
plt.xlim(bmin, bmax)


# In[ ]:

#
#  EXACT TEST FOR 1st ORDER
#

NbPoints = 10000
xmin = -5.
xmax =  5.
Lmax = np.abs(xmax - xmin)
dx = Lmax/NbPoints
print("dx = %f" %dx)

x = Getx(NbPoints, xmin, xmax)

x0 = 2.
dx_L0 = 1.
dx_L1 = dx_L0/2.
RF = dx_L0/dx_L1

shape = GetOrd1(x, x0, dx_L0)
order = 2
nbpts = (2*int(RF) - 1) + (order-1)

w1 = 1./2. 
w2 = 1. 
w3 = 1./2. 

print("w1 = %.2f"%(w1))
print("w2 = %.2f"%(w2))
print("w3 = %.2f"%(w3))

wtot = w1 + w2 + w3 
print("wtot = %.2f"%wtot)

print("RF = %.2f"%RF)

s1 = w1*GetOrd1(x, x0-dx_L1, dx_L1)
s2 = w2*GetOrd1(x, x0      , dx_L1)
s3 = w3*GetOrd1(x, x0+dx_L1, dx_L1)

shape_split = s1 + s2 + s3 

# Compute the area using the composite Simpson's rule.
area_ref = scint.simps(shape, dx=dx)

As1 = scint.simps(s1, dx=dx)
As2 = scint.simps(s2, dx=dx)
As3 = scint.simps(s3, dx=dx)

print("As1 = %f"%(As1))
print("As2 = %f"%(As2))
print("As3 = %f"%(As3))

area_split = As1 + As2 + As3 

print("area ref=", area_ref)
print("area split=", area_split)

print("After normalization :")
area_split = scint.simps(shape_split, dx=0.001)
print("area split=", area_split)

bmin = x0 - 0.5*nbpts*dx_L1 - (order+1)*0.5*dx_L1
bmax = x0 + 0.5*nbpts*dx_L1 + (order+1)*0.5*dx_L1
plt.plot(x,shape)
plt.plot(x,shape_split,'red')
plt.plot(x,s1,'r--')
plt.plot(x,s2,'r--')
plt.plot(x,s3,'r--')

plt.xlim(bmin, bmax)
    
 
    

# In[83]:

#
#  EXACT TEST FOR 2nd ORDER
#

NbPoints = 10000
xmin = -5.
xmax =  5.
Lmax = np.abs(xmax - xmin)
dx = Lmax/NbPoints
print("dx = %f" %dx)

x = Getx(NbPoints, xmin, xmax)

x0 = 2.
dx_L0 = 1.
dx_L1 = dx_L0/2.
RF = dx_L0/dx_L1

dta  = 0.5 *dx_L1
dta2 = 1.5*dx_L1

shape = GetOrd2(x, x0, dx_L0)
order = 2
nbpts = (2*int(RF) - 1) + (order-1)

w1 = 3./4. # x0-dta
w2 = 3./4. # x0-dta
w3 = 1./4. # x0-dta2
w4 = 1./4. # x0-dta2

print("w1 = %.2f, w2 = %.2f"%(w1, w2))
print("w3 = %.2f, w4 = %.2f"%(w3, w4))

wtot = w1 + w2 + w3 + w4 
print("wtot = %.2f"%wtot)

print("RF = %.2f"%RF)
print("w1 = %.2f, w2 = %.2f"%(w1, w2))
print("w3 = %.2f, w4 = %.2f"%(w3, w4))

s1 = w1*GetOrd2(x, x0-dta, dx_L1)
s2 = w2*GetOrd2(x, x0+dta, dx_L1)
s3 = w3*GetOrd2(x, x0-dta2, dx_L1)
s4 = w4*GetOrd2(x, x0+dta2, dx_L1)

shape_split = s1 + s2 + s3 + s4

# Compute the area using the composite Simpson's rule.
area_ref = scint.simps(shape, dx=dx)

As1 = scint.simps(s1, dx=dx)
As2 = scint.simps(s2, dx=dx)
As3 = scint.simps(s3, dx=dx)
As4 = scint.simps(s4, dx=dx)

print("As1 = %f, As2 = %f"%(As1, As2))
print("As3 = %f, As4 = %f"%(As3, As4))

area_split = As1 + As2 + As3 + As4  

print("area ref=", area_ref)
print("area split=", area_split)

print("After normalization :")
area_split = scint.simps(shape_split, dx=0.001)
print("area split=", area_split)

bmin = x0 - 0.5*nbpts*dx_L1 - (order+1)*0.5*dx_L1
bmax = x0 + 0.5*nbpts*dx_L1 + (order+1)*0.5*dx_L1
plt.plot(x,shape)
plt.plot(x,shape_split,'red')
plt.plot(x,s1,'r--')
plt.plot(x,s2,'r--')
plt.plot(x,s3,'r--')
plt.plot(x,s4,'r--')
#plt.plot(x,s5,'r--')
plt.xlim(bmin, bmax)


# In[ ]:




# In[ ]:




# In[85]:

NbPoints = 10000
xmin = -5.
xmax =  5.
Lmax = np.abs(xmax - xmin)
dx = Lmax/NbPoints
print("dx = %f" %dx)

x = Getx(NbPoints, xmin, xmax)

x0 = 2.
dx_L0 = 1.
dx_L1 = dx_L0/2.
RF = dx_L0/dx_L1

shape = GetOrd3(x, x0, dx_L0)
order = 2
nbpts = (2*int(RF) - 1) + (order-1)

w1 = 1./8. # x0 - 2*dx1
w2 = 4./8. # x0 - dx1
w3 = 6./8. # x0 
w4 = 4./8. # x0 + dx1
w5 = 1./8. # x0 + 2*dx1

print("w1 = %.2f, w2 = %.2f"%(w1, w2))
print("w3 = %.2f"%(w3))
print("w4 = %.2f, w5 = %.2f"%(w4, w5))

wtot = w1 + w2 + w3 + w4 + w5
print("wtot = %.2f"%wtot)

#w1 = w1/wtot * RF
#w2 = w2/wtot * RF 
#w3 = w3/wtot * RF 
#w4 = w4/wtot * RF 
#w5 = w5/wtot * RF 

print("RF = %.2f"%RF)
print("w1 = %.2f, w2 = %.2f"%(w1, w2))
print("w3 = %.2f"%(w3))
print("w4 = %.2f, w5 = %.2f"%(w4, w5))

s1 = w1*GetOrd3(x, x0 - 2.*dx_L1, dx_L1)
s2 = w2*GetOrd3(x, x0 -    dx_L1, dx_L1)
s3 = w3*GetOrd3(x, x0           , dx_L1)
s4 = w4*GetOrd3(x, x0 +    dx_L1, dx_L1)
s5 = w5*GetOrd3(x, x0 + 2.*dx_L1, dx_L1)

shape_split = s1 + s2 + s3 + s4 + s5

# Compute the area using the composite Simpson's rule.
area_ref = scint.simps(shape, dx=dx)

As1 = scint.simps(s1, dx=dx)
As2 = scint.simps(s2, dx=dx)
As3 = scint.simps(s3, dx=dx)
As4 = scint.simps(s4, dx=dx)
As5 = scint.simps(s5, dx=dx)

print("As1 = %f, As2 = %f"%(As1, As2))
print("As3 = %f"%(As3))
print("As4 = %f, As5 = %f"%(As4, As5))

area_split = As1 + As2 + As3 + As4 + As5 

print("area ref=", area_ref)
print("area split=", area_split)

print("After normalization :")
area_split = scint.simps(shape_split, dx=0.001)
print("area split=", area_split)

bmin = x0 - 0.5*nbpts*dx_L1 - (order+1)*0.5*dx_L1
bmax = x0 + 0.5*nbpts*dx_L1 + (order+1)*0.5*dx_L1
plt.plot(x,shape)
plt.plot(x,shape_split,'red')
plt.plot(x,s1,'r--')
plt.plot(x,s2,'r--')
plt.plot(x,s3,'r--')
plt.plot(x,s4,'r--')
plt.plot(x,s5,'r--')
plt.xlim(bmin, bmax)





# In[ ]:

# GENERAL SCRIPT FOR 1st order 
# AND ARBITRARY REFINEMENT FACTOR (RF) 

NbPoints = 10000
xmin = -5.
xmax =  5.
Lmax = np.abs(xmax - xmin)
dx = Lmax/NbPoints
print("dx = %f" %dx)

x = Getx(NbPoints, xmin, xmax)

x0 = 1.
RF = 3.
dx_L0 = 1.
dx_L1 = dx_L0/RF

# 1st order
order=1
shape = GetOrd1(x, x0, dx_L0)
nbpts = (2*int(RF) - 1) + (order-1)

print("RF = %.2f"%RF)
print("order = %d"%order)

print("nbpts = (2*RF - 1) + (order-1)")
print("nbr of points = %d"%nbpts)


w=np.zeros(nbpts)

# Calcul des coefficients
ik_l = np.arange(nbpts)
for ik in np.nditer(ik_l): print("%d "%ik, end='')
print("")

inum=1
for ik in ik_l:
    w[ik] = float(inum)/float(RF)
    
    if ik-RF+1 < 0:
        inum=inum+1
    else:
        inum=inum-1

for wi in np.nditer(w): print("%f "%wi, end='')
print("")
        
wtot = np.sum(w)
print("wtot = %.2f"%wtot)

split_tab = np.zeros((nbpts, len(x)))

inum=-(RF-1)
for ik in ik_l:
    split_tab[ik,:] = w[ik]*GetOrd1(x, x0 + float(inum)*dx_L1, dx_L1)
    print(inum)
    
    inum=inum+1

shape_split = split_tab.sum(axis=0)

# Compute the area using the composite Simpson's rule.
area_ref = scint.simps(shape, dx=dx)

Asplit_tab = np.zeros(nbpts)

for ik in ik_l:
    Asplit_tab[ik] = scint.simps(split_tab[ik,:], dx=dx)
    print("Asplit_tab[%d] = %f"%(ik, Asplit_tab[ik]))

print("area ref=", area_ref)
    
area_split = np.sum(Asplit_tab)
print("Sum split areas =", area_split)

area_split = scint.simps(shape_split, dx=0.001)
print("Sum split areas =", area_split)

bmin = x0 - 0.5*nbpts*dx_L1 - (order+1)*0.5*dx_L1
bmax = x0 + 0.5*nbpts*dx_L1 + (order+1)*0.5*dx_L1
plt.plot(x,shape)
plt.plot(x,shape_split,'red')

for ik in ik_l:
    plt.plot(x,split_tab[ik,:],'r--')

plt.xlim(bmin, bmax)
plt.ylim(0, 1)

plt.xlabel(r"$\xi$",size=20) ; 
plt.ylabel(r"$S(\xi)$",size=20) 


# In[ ]:



# In[376]:

# GENERAL SCRIPT FOR 
# ARBITRARY ORDER AND RF = 2

NbPoints = 10000
xmin = -5.
xmax =  5.
Lmax = np.abs(xmax - xmin)
dx = Lmax/NbPoints
print("dx = %f" %dx)

x = Getx(NbPoints, xmin, xmax)

# We set the refinement ratio
# and the shape function order
RF = 2.
x0 = 1.5
dx_L0 = 1.
dx_L1 = dx_L0/RF

# Arbitrary order
order = 2
shape = getShape(order, x, x0, dx_L0)

print("RF = %.2f"%RF)
print("order = %d"%order)

nbpts = order+2
print("nbpts = 2*(order+1)+1")
print("nbr of points = %d"%nbpts)

nt = order+1
# Calcul des coefficients
w = np.array([float(binomial(nt, y)) for y in range(0, nt+1)])

w[:] = w[:]/RF**order
for wi in np.nditer(w): print("%f "%wi, end='')
print("")

ik_l = np.arange(nbpts)
for ik in np.nditer(ik_l): print("%d "%ik, end='')
print("")
        
wtot = np.sum(w)
print("wtot = %.2f"%wtot)

split_tab = np.zeros((nbpts, len(x)))

inum=-(nbpts-1)/2.
for ik in ik_l:
    split_tab[ik,:] = w[ik]*getShape(order, x, x0 + float(inum)*dx_L1, dx_L1)
    print(inum)
    
    inum=inum+1

shape_split = split_tab.sum(axis=0)

# Compute the area using the composite Simpson's rule.
area_ref = scint.simps(shape, dx=dx)

Asplit_tab = np.zeros(nbpts)

for ik in ik_l:
    Asplit_tab[ik] = scint.simps(split_tab[ik,:], dx=dx)
    print("Asplit_tab[%d] = %f"%(ik, Asplit_tab[ik]))

print("area ref=", area_ref)
    
area_split = np.sum(Asplit_tab)
print("Sum split areas =", area_split)

area_split = scint.simps(shape_split, dx=0.001)
print("Sum split areas =", area_split)

bmin = x0 - 0.5*nbpts*dx_L1 - (order+1)*0.5*dx_L1
bmax = x0 + 0.5*nbpts*dx_L1 + (order+1)*0.5*dx_L1
plt.plot(x,shape)
plt.plot(x,shape_split,'red')

for ik in ik_l:
    plt.plot(x,split_tab[ik,:],'r--')

plt.xlim(bmin, bmax)
#plt.legend([r"$L_0$, %.2f dx : $\rho$"%delta, r"$L_1$, %.2f dx : $\rho$"%delta])  

plt.xlabel(r"$\xi$",size=20) ; 
plt.ylabel(r"$S(\xi)$",size=20) 


# In[ ]:




# In[ ]:




# In[374]:

# GENERAL SCRIPT FOR 
# ARBITRARY ORDER AND RF = 3

NbPoints = 10000
xmin = -5.
xmax =  5.
Lmax = np.abs(xmax - xmin)
dx = Lmax/NbPoints
print("dx = %f" %dx)

x = Getx(NbPoints, xmin, xmax)

# We set RF=3
# Specialized script for RF=3
# and arbitrary order
RF = 3.
x0 = 2.
dx_L0 = 1.
dx_L1 = dx_L0/RF

# Arbitrary order
order=3
shape = getShape(order, x, x0, dx_L0)

print("RF = %.2f"%RF)
print("order = %d"%order)

nbpts = 2*(order+1)+1
print("nbpts = 2*(order+1)+1")
print("nbr of points = %d"%nbpts)

nt = order+1
# Calcul des coefficients
w = np.array([float(trinomial(nt, y)) for y in range(-nt, nt+1)])
print(w)

w[:] = w[:]/RF**order
print(w)

ik_l = np.arange(nbpts)
print(ik_l)

        
wtot = np.sum(w)
print("wtot = %.2f"%wtot)

split_tab = np.zeros((nbpts, len(x)))

inum=-(nbpts-1)/2.
for ik in ik_l:
    split_tab[ik,:] = w[ik]*getShape(order, x, x0 + float(inum)*dx_L1, dx_L1)
    print(inum)
    
    inum=inum+1

shape_split = split_tab.sum(axis=0)

# Compute the area using the composite Simpson's rule.
area_ref = scint.simps(shape, dx=dx)

Asplit_tab = np.zeros(nbpts)

for ik in ik_l:
    Asplit_tab[ik] = scint.simps(split_tab[ik,:], dx=dx)
    print("Asplit_tab[%d] = %f"%(ik, Asplit_tab[ik]))

print("area ref=", area_ref)
    
area_split = np.sum(Asplit_tab)
print("Sum split areas =", area_split)

area_split = scint.simps(shape_split, dx=0.001)
print("Sum split areas =", area_split)

bmin = x0 - 0.5*nbpts*dx_L1 - (order+1)*0.5*dx_L1
bmax = x0 + 0.5*nbpts*dx_L1 + (order+1)*0.5*dx_L1
plt.plot(x,shape)
plt.plot(x,shape_split,'red')

for ik in ik_l:
    plt.plot(x,split_tab[ik,:],'r--')

plt.xlim(bmin, bmax)

plt.xlabel(r"$\xi$",size=20) ; 
plt.ylabel(r"$S(\xi)$",size=20) 

# In[139]:




# In[ ]:




# In[346]:

def f(x):
    return {
        'a': 1,
        'b': 2,
    }.get(x, 9)


# In[348]:

f('g')


# In[349]:




# In[ ]:




# In[203]:




# In[ ]:




# In[ ]:




# In[229]:

ik_l = np.arange(5)+1
print(ik_l)


# In[254]:

j_l = np.arange(n+1)
print(j_l)


# In[ ]:




# In[ ]:




# In[238]:

a=range(-3, 3+1)
print(a)




# In[ ]:




# In[204]:

n=2

for k in range(2*(n+1)):
    print( trinomial(n,k) )


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[92]:

aa=np.zeros(5)


# In[93]:

print(aa)


# In[98]:

np.size(aa)


# In[101]:

ind=np.arange(5)+1


# In[106]:

ind


# In[108]:

np.sum(ind)


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



