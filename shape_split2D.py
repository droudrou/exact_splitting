
# coding: utf-8

# In[1]:

import mpl_toolkits.mplot3d as mplt3d
import matplotlib.pyplot as plt
import scipy.integrate as scint
import math
import numpy as np
#%matplotlib inline


# In[ ]:




# In[19]:

def Getx(nx, xmin, xmax):
    dx = abs(xmax-xmin)/float(nx)
    x  = np.arange(nx+1)*dx + xmin
    return x


# In[59]:

# reduce is not recognized by python 3.5.2
#def fact(n):return reduce(lambda x,y:x*y,[1]+range(1,n+1))

def fact(n):return math.factorial(n)

# In[60]:

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



# In[61]:

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


# In[ ]:




# In[ ]:




# In[20]:

def GetOrd1(x,x0,dx):
    
    xi=(x-x0)/dx

    izeros = np.argwhere(abs(xi) > 1.)
    ileft  = np.argwhere(xi < 0.)
    iright = np.argwhere(xi >= 0.)
    
    shape = np.zeros(len(x))
    
    shape[ileft] = xi[ileft]+1.

    shape[iright] = 1.-xi[iright]
    
    shape[izeros] = 0.
    
    return (shape)


# In[21]:

def GetOrd2(x,x0,dx):
    
    xi=(x-x0)/dx

    iequ1  = np.argwhere( abs(xi) < 0.5  )
    iequ2  = np.argwhere( abs(xi) >= 0.5 ) 
    izeros = np.argwhere( abs(xi) > 1.5  )
    
    shape = np.zeros(len(x))
       
    shape[iequ1] = 0.75 - abs(xi[iequ1])**2.
    shape[iequ2] = 0.5*(1.5 - abs(xi[iequ2]))**2.
    shape[izeros] = 0.
    
    return (shape)


# In[22]:

def GetOrd3(x,x0,dx):
    
    xi=(x-x0)/dx

    iequ1  = np.argwhere( abs(xi) < 1.  )
    iequ2  = np.argwhere( abs(xi) >= 1. ) 
    izeros = np.argwhere( abs(xi) > 2.  )
    
    shape = np.zeros(len(x))
    
    shape[iequ1] = 0.5*abs(xi[iequ1])**3. - xi[iequ1]**2. + 2./3.
    
    shape[iequ2] = 4./3.*(1. - 0.5*abs(xi[iequ2]))**3.
    
    shape[izeros] = 0.
    
    return (shape)


# In[23]:

def GetOrd4(x,x0,dx):
    
    xi=(x-x0)/dx

    iequ1  = np.argwhere( abs(xi) <  0.5 )
    iequ2  = np.argwhere( abs(xi) >= 0.5 ) 
    iequ3  = np.argwhere( abs(xi) >= 1.5 ) 
    izeros = np.argwhere( abs(xi) >  2.5 )
    
    shape = np.zeros(len(x))
    
    shape[iequ1] = 1./192.*( 115. - 120.*xi[iequ1]**2. + 48.*xi[iequ1]**4. )
    
    shape[iequ2] = 1./96.*( 55. + 20.*abs(xi[iequ2]) - 120.*xi[iequ2]**2. + 80.*abs(xi[iequ2])**3. - 16.*xi[iequ2]**4.)
    
    shape[iequ3] = 1./24.*(2.5 - abs(xi[iequ3]))**4.
    
    shape[izeros] = 0.
    
    return (shape)


# In[ ]:




# In[25]:

bmin = 2.
bmax = 6.

x = Getx(10000, bmin, bmax)
shape = GetOrd2(x, 4., 1.)

plt.plot(x,shape)
#plt.xlim(bmin, bmax)

plt.show()


# In[ ]:




# In[26]:

fig = plt.figure()
ax = mplt3d.Axes3D(fig)
X = np.arange(-4, 4, 0.25)
Y = np.arange(-4, 4, 0.25)
X, Y = np.meshgrid(X, Y)

R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)

ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='hot')
plt.show()


# In[ ]:




# In[9]:

X = np.arange(-4, 4, 0.25)
print(X)


# In[ ]:


# In[49]:

def GetOrd2_2d(x_vect, y_vect, x0, y0, dx, dy):
    
    Xshape = GetOrd2(x_vect, x0, dx).reshape(np.size(x_vect), 1  )
    Yshape = GetOrd2(y_vect, y0, dy).reshape(1  , np.size(y_vect))

    shape = Xshape * Yshape
    
    return (shape)

# In[ ]:



# In[50]:

fig = plt.figure( figsize=(18,10) )
ax = mplt3d.Axes3D(fig)

bmin = 2.
bmax = 6.

order = 2
dx_L0 = 1.
dx_L1 = dx_L0/2.
RF = dx_L0/dx_L1

x0 = 4.
y0 = 4.

nbcellx = 100
nbcelly = 100

nbx = nbcellx+1
nby = nbcelly+1

x_vect = Getx(nbcellx, bmin, bmax)
y_vect = Getx(nbcelly, bmin, bmax)

x2, y2 = np.meshgrid(x_vect, y_vect)
print(x2.shape)
print(y2.shape)

#Xshape = GetOrd2(x_vect, x0, dx_L0).reshape(nbx, 1  )
#Yshape = GetOrd2(y_vect, y0, dx_L0).reshape(1  , nby)
#shape = Xshape * Yshape

shape = GetOrd2_2d(x_vect, y_vect, x0, y0, dx_L0, dx_L0)

print(shape.shape)

ax.plot_surface(x2, y2, shape, rstride=1, cstride=1, cmap='hot')

plt.show(fig)


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:







# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[83]:

bmin = 2.
bmax = 6.

nbcellx = 100
nbcelly = 100

order = 2
dx_L0 = 1.
dx_L1 = dx_L0/2.
RF = dx_L0/dx_L1

x0 = 4.
y0 = 4.

nbx = nbcellx+1
nby = nbcelly+1

x = Getx(nbcellx, bmin, bmax)
y = Getx(nbcelly, bmin, bmax)

xm, ym = np.meshgrid(x, y)
print(xm.shape)
print(ym.shape)

shape = GetOrd2_2d(x, y, x0, y0, dx_L0, dx_L0)

nt = order+1
# Calcul des coefficients
w = np.array([float(binomial(nt, bk)) for bk in range(0, nt+1)])
print(w)

w_2d = w.reshape(order+2, 1) * w.reshape(1, order+2)
print(w_2d)

w_2d[:,:] = w_2d[:,:]/RF**(2*order)
print(w_2d)

wtot = np.sum(w)
print("wtot = %.2f"%wtot)

nbpts = int( (2*RF - 1) + (RF-1)*(order-1) )
print("nbpts = %d"%nbpts)
ik_l = np.arange(nbpts)
print(ik_l)
jk_l = np.arange(nbpts)
print(jk_l)
iNew_l = np.arange(nbpts**2)

nb_2d = nbpts**2
print("nb_2d = %d"%nb_2d)

xvec = np.array( np.zeros(nbpts) )
yvec = np.array( np.zeros(nbpts) )

ix = -(nbpts-1)/2.
for ik in ik_l:
    xvec[ik] = x0 + ix*dx_L1
    yvec[ik] = y0 + ix*dx_L1
    ix = ix + 1.

print(xvec)
print(yvec)


split_tab = np.zeros((nb_2d, len(x), len(y)))

shape_split = np.zeros((len(x), len(y)))

iNew = 0
for ik in ik_l:
    for jk in jk_l:
        #shape_split[:, :] = shape_split[:,:] + w_2d[ik,jk]*GetOrd2_2d(x, y, xvec[ik], yvec[jk], dx_L1, dx_L1)
        split_tab[iNew, :, :] = w_2d[ik,jk]*GetOrd2_2d(x, y, xvec[ik], yvec[jk], dx_L1, dx_L1)
        shape_split[:, :] = shape_split[:,:] + split_tab[iNew, :, :]
        iNew = iNew+1

#split_tab[ik, :, :] = w_2d[ik,jk]*GetOrd2_2d(x, y, xvec[ik], yvec[ik], dx_L1, dx_L1)


delta_shape = np.zeros((len(x), len(y)))
delta_shape[:, :] = shape[:, :] - shape_split[:, :] 


fig = plt.figure( figsize=(16,10) )
ax = mplt3d.Axes3D(fig)
ax.plot_surface(xm, ym, shape, rstride=1, cstride=1, cmap='hot')

fig = plt.figure( figsize=(16,10) )
ax = mplt3d.Axes3D(fig)
ax.plot_surface(xm, ym, delta_shape, rstride=1, cstride=1, cmap='hot')

#i2d_l = np.arange(nb_2d)
#print(i2d_l)
#fig = plt.figure( figsize=(18,10) )
#for i2d in i2d_l:
##    fig = plt.figure( figsize=(18,10) )
#    ax = mplt3d.Axes3D(fig)
#    ax.plot_surface(xm, ym, split_tab[i2d, :, :], rstride=1, cstride=1, cmap='hot')
#
#plt.show()


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



