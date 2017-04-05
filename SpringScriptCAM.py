
# coding: utf-8

# # CAM-Compression Spring-Follower
# ### INME 4012
# #### Christian J. Lagares Nieves
# 
# Illustrated below is a circular disk cam assembly with a helical compression spring. The spring has plain ground ends. The cam is a circle with an eccentricity e = 0.75 in. The total follower displacement is 2e or 1.5 in. The length of the spring when the follower is fully extended is 6.75 in. A rod with a 1.5 in diameter passes through the center of the helical compression spring. The spring force in the BDC (Bottom Dead Center) load position is 150 lb and the spring force in the TDC (Top Dead Center) load position is 300 lb. The disk cam rotates at 1,200 RPM. The spring material is made of chrome-vanadium wire that has been shot peened.
# 
# List all assumptions and determine:
# 
# * Wire diameter d
# * Mean Coil Diameter D
# * Number of active turns N
# * Fundamental Frequency f
# * Spring rate k
# * Free length of spring L<sub>f</sub>
# * Solid length of spring L<sub>s</sub>
# * Fatigue safety factor

# ![CAM-Spring Figure](http://localhost:8888/notebooks/Desktop/Screen%20Shot%202017-04-05%20at%209.43.23%20AM.png)

# ## Modules Import

# In[1]:

import numpy as np
#get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
import pandas as pd


# ## Main Parameters

# In[2]:

Fmax = 300     # lbf
Fmin = 150     # lbf
e = 0.75       # in
omega = 1200   # RPM
minD = 1.5     # in
maxComp = 6.75 # in


# ## Procedure

# In order to determine k, we can simply find the relation between the maximum displacement/force and minimum, or initial, displacement/force.
# 
# $$F = k \delta$$
# 
# $$k = \frac{F}{\delta}$$

# In[3]:

dF = Fmax-Fmin # lbf
delta = 2*e    # in
k = dF/delta   # lbf/in
print("k =",k,"lbf/in")


# In[4]:

drange = np.zeros((1000,),dtype=np.float32)
sultrange = np.zeros((1000,),dtype=np.float32)


# In[5]:

# Chrome Vanadium Wire
mind = 0.032
maxd = 0.437
A = 169000 
m = 0.168
G = 11.2 # Mpsi

drange = np.array(range(len(drange)))*((maxd-mind)/len(drange))+0.032
sultrange = A/(drange**m)
syrange = 0.45*sultrange


# In[6]:

ultplot = plt.plot(drange,sultrange,label="$S_{ult}$")
yieplot = plt.plot(drange,syrange,label="$S_y$")
plt.legend()
plt.xlabel("d [in]")
plt.ylabel("S [ksi]")
plt.grid(True, which='major',linestyle='--',color='k')
plt.show()


# ## Determining C and D
# 
# Given the provided parameters, D will be computed parting from an allowance and the follower's diameter.
# 
# $$D = d_{rod}+d+allowance$$
# 
# $$C = \frac{D}{d}$$

# In[7]:

allowance = 0.125
Drange = minD + drange + allowance

C = Drange/drange


# In[8]:

Cplot = plt.plot(drange,C,label="$C$")
plt.legend()
plt.xlabel("d")
plt.grid(True, which='major',linestyle='--',color='k')
plt.show()


# In[9]:

KB = ((4*C)+2)/((4*C)-3)
notComply = []
for ii, Cs in enumerate(C):
    if Cs > 12:
        notComply.append(ii)
    elif Cs < 4:
        notComply.append(ii)

d = np.delete(drange,notComply)
D = np.delete(Drange,notComply)
KB = np.delete(KB, notComply)
sult = np.delete(sultrange,notComply)
sy = np.delete(syrange,notComply)
C = np.delete(C,notComply)


# ## Computing Stress 
# 
# $$\tau_{s} = \frac{8K_B(1+\epsilon)F_{max}D}{\pi d^3} $$
# 

# In[10]:

epsilon = 0.15

taus = (8*KB*(1+epsilon)*Fmax*D)/(np.pi*(d**3))

ns = sy/taus


# Deleting all values along the arrays corresponding to a safety factor below 1.2.

# In[11]:

notComply = []
for ii, Ns in enumerate(ns):
    if Ns < 1.2:
        notComply.append(ii)

d = np.delete(drange,notComply)
D = np.delete(Drange,notComply)
KB = np.delete(KB, notComply)
sult = np.delete(sultrange,notComply)
sy = np.delete(syrange,notComply)
C = np.delete(C,notComply)
taus = np.delete(taus, notComply)
ns = np.delete(ns,notComply)


# ## Determine number of Active Turns
# 
# From the previously determined value of $k$, one can determine the required number of active turns.
# 
# $$k = \frac{d^4 G}{8D^3N}$$
# 
# Solving for N,
# 
# $$N_a = \frac{d^4 G}{8D^3k}$$

# In[12]:

Na = ((d**4)*G)/(8*(D**3)*k)


# In[13]:

# Plain and Ground
Ne = 1
Nt = Na + Ne

Ls = d*Nt

a = 0.5
clearance = 0.25

Lo = 2.63*(D/a) - clearance
p = Lo/(Na+1)


# In[14]:

delta = Lo-maxComp
Ferror = 100*(np.absolute((delta*k)-Fmax)/Fmax)


# In[15]:

Dpossible = []
dpossible = []
KBpossible = []
sultpossible = []
sypossible = []
Cpossible = []
tauspossible = []
nspossible = []
Napossible = []
Ntpossible = []
Lopossible = []
Lspossible = []
ppossible = []
Ferrorpossible = []

for ii, error in enumerate(Ferror):
    if error < 1:
        Dpossible.append(D[ii])
        dpossible.append(d[ii])
        KBpossible.append(KB[ii])
        sultpossible.append(sult[ii])
        sypossible.append(sy[ii])
        Cpossible.append(C[ii])
        tauspossible.append(taus[ii])
        nspossible.append(ns[ii])
        Napossible.append(Na[ii])
        Ntpossible.append(Nt[ii])
        Lopossible.append(Lo[ii])
        Lspossible.append(Ls[ii])
        ppossible.append(p[ii])
        Ferrorpossible.append(Ferror[ii])


# In[16]:

PossibleData = pd.DataFrame({
    'D': Dpossible,
    'd': dpossible,
    'Kb': KBpossible,
    'Sult': sultpossible,
    'Sy': sypossible,
    'C': Cpossible,
    'Tau': tauspossible,
    'ns': nspossible,
    'Na': Napossible,
    'Nt': Ntpossible,
    'Lo': Lopossible,
    'Ls': Lspossible,
    'p': ppossible,
    'F error': Ferrorpossible
})
PossibleData


# In[17]:

Dmin = []
dmin = []
KBmin = []
sultmin = []
symin = []
Cmin = []
tausmin = []
nsmin = []
Namin = []
Ntmin = []
Lomin = []
Lsmin = []
pmin = []
Ferrormin = []

for ii, error in enumerate(Ferror):
    if error == min(Ferror):
        Dmin.append(D[ii])
        dmin.append(d[ii])
        KBmin.append(KB[ii])
        sultmin.append(sult[ii])
        symin.append(sy[ii])
        Cmin.append(C[ii])
        tausmin.append(taus[ii])
        nsmin.append(ns[ii])
        Namin.append(Na[ii])
        Ntmin.append(Nt[ii])
        Lomin.append(Lo[ii])
        Lsmin.append(Ls[ii])
        pmin.append(p[ii])
        Ferrormin.append(Ferror[ii])


# In[18]:

OptimalData = pd.DataFrame({
    'D': Dmin,
    'd': dmin,
    'Kb': KBmin,
    'Sult': sultmin,
    'Sy': symin,
    'C': Cmin,
    'Tau': tausmin,
    'ns': nsmin,
    'Na': Namin,
    'Nt': Ntmin,
    'Lo': Lomin,
    'Ls': Lsmin,
    'p': pmin,
    'F Error': Ferrormin
})
print("Taking as optimum the Free Length minimizing the resulting error force, \n")
OptimalData


# ## Dynamic Analysis
# 
# The fundamental frequency for the current scenario can be found as:
# 
# $$f = \frac{1}{4} \sqrt{\frac{kg}{W}}$$
# 
# where W can be found as,
# 
# $$W = \rho (\pi N_t D) \frac{\pi}{4}d^2$$
# 

# In[19]:

density = 0.284 #lb/in^3
Wopt = density * (np.pi*np.array(Ntmin)*np.array(Dmin)) *((np.pi/4)*np.array(dmin)**2)
Wpos = density * (np.pi*np.array(Ntpossible)*np.array(Dpossible)) *((np.pi/4)*np.array(dpossible)**2)
print("Weight of optimal solution:",Wopt[0], "lb")


# In[20]:

g = 386.0886 #in/s^2
fopt = 0.25 * np.sqrt((k*g)/Wopt)
print("Fundamental Frequency:",fopt[0],"Hz")
print("Operating Frequency:", omega//60, "Hz")


# In[21]:

taua = np.array(KBmin)*((8*((Fmax-Fmin)/2)*np.array(Dmin))/(np.pi*(np.array(dmin)**3)))
taum = np.array(KBmin)*((8*((Fmax+Fmin)/2)*np.array(Dmin))/(np.pi*(np.array(dmin)**3)))


Ssa = 57500
Ssm = 77500
Ssu = 0.67*np.array(sultmin)
Sse = Ssa/(1-(Ssm/Ssu)**2)
r = taua/taum

Ssa = (((r**2) * (Ssu**2))/(2*Sse))*(-1+np.sqrt(1+(((2*Sse)/(r*Ssu))**2)))
Ssa
print("Ssa =",Ssa[0]/1000,"kpsi")


# In[22]:

nf = Ssa/taua
print("nf =",nf[0])


# In[23]:

print("k =",k,"lbf/in")
print("d =",dmin[0],"in")
print("D =",Dmin[0],"in")
print("Lf =",Lomin[0],"in")
print("Ls =",Lsmin[0],"in")
print("Na =",Namin[0],"turns")
print("Nt =",Ntmin[0],"turns")
print("f =",fopt[0],"Hz")
print("nf =", nf[0])

