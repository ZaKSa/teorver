import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.stats import norm
from scipy.integrate import quad
from numpy import *
from random import *
from math import *


gamma2=6.76
gamma = sqrt(gamma2)
n = 350
a = 2.2
alfa = 0.001
y=0.99

xi = [0 for i in range(n)]

input = open(r"C:\Users\Ksenia\terver.txt", "r")
for i in range(n):
    xi[i] = double(input.readline())
input.close()

text = open(r"C:\Users\Ksenia\output.txt","w")

#количество интервалов 4
N = floor(1+ 3.32*log(n, 10)) + 1
text.write("N> "+str(N)+"\n")
#Подсчитываем частоты Vk попадания выборочных значений в k-ый интервал kJ
u=(max(xi)+1-min(xi))/N
text.write("u> "+str(u)+"\n")
for i in range (N):
    text.write("> "+str(min(xi)+u*i)+"\n")
text.write(str(max(xi))+"\n")


vk = [0]*N
for k in xi:
    l = (k-min(xi))//u
    vk[int(l)] += 1
text.write("chastota> "+str(vk)+"\n")


pk=[0]*N
for i in range(N):
    pk[i]=vk[i]/n
text.write("otnosit chastota> "+str(pk) +"\n")


hk=[0]*N
for i in range(N):
    hk[i]=vk[i]/(u*n)
text.write("visota> "+str(hk) +"\n")


for i in range (N):
    mid = u/2*(i+1)
    f = math.exp((mid-a)**2/(2*gamma2))/(gamma*sqrt(2*math.pi))
    text.write("teor znachenie f> "+str(f) + "\n")


x = np.linspace(a-8, a+8, 1000)
pdf = norm.pdf(x, loc = a, scale=gamma)

plt.hist(xi, bins = N, density=True)
plt.plot(x, pdf)
plt.show()
#plt.hist(xi, bins=np.arange(min(xi), max(xi)+1))
#plt.show()

text.write("\n")
x = 0
for i in range(n):#выборочное среднее
    x += xi[i] 
x/=n
text.write("mx1= " + str(x)+"\n")


dx = 0
for i in range(n):#выборочная дисперсия смещенная
    dx+=(xi[i]-x)**2
dx/=n
text.write("dx1= " +str(dx)+"\n")

dx_no = 0
for i in range(n):#выборочная дисперсия несмещенная
    dx_no+=(xi[i]-x)**2
dx_no/=(n-1)
text.write("dx1 nesmeshennaya= " +str(dx_no)+"\n")


x2 = 0
for i in range(N):#
    mid = ((min(xi)+u*i)+min(xi)+u*(i+1))/2
    x2 += vk[i]*mid
x2/=n
text.write("mx2= " +str(x2)+"\n")


dx2 = 0
for i in range(N):#
    mid = ((min(xi)+u*i)+min(xi)+u*(i+1))/2
    dx2+=(mid-x2)**2*vk[i]
dx2/=n
text.write("dx2= " + str(dx2)+"\n")


text.write("mx= " + str(a)+"\n")
text.write("dx= "+str(gamma2)+"\n")

#1.3
mx3 = 0
for i in range(n):
    mx3+=xi[i]
mx3/=n
text.write("mx3= "+str(mx3)+"\n")

dx3 = 0
for i in range(n):
    dx3+=(xi[i]-mx3)**2
dx3/=n
text.write("dx3= "+str(dx3)+"\n")


#1.4 (стр 15)
#дов интервал для М(неизв М, изв Д)
c=2.576
t1_1 = x-c*gamma/sqrt(n)
t2_1 = x+c*gamma/sqrt(n)
text.write("дов интервал для М(неизв М, изв Д)\n"+"t1_1= "+str(t1_1)+"\n")
text.write("t2_1= "+str(t2_1)+"\n")

#дов интервал для М(неизв М, неизв Д)
t=c
t1_2=x - t*dx/sqrt(n-1)
t2_2=x + t*dx/sqrt(n-1)
text.write("дов интервал для М(неизв М, неизв Д)\n"+"t1_2= "+str(t1_2)+"\n")
text.write("t2_2= "+str(t2_2)+"\n")

#дов интервал для Д(изв М) для наблюдаемой св
#p1 =0.995, c= 2.576, 
#p2=0.005, p2=-0.995
hi1=((c+sqrt(2*n-1))**2)/2
hi2=((-c+sqrt(2*n-1))**2)/2
t1_3=0
t2_3=0
for i in range(n):
    t1_3+=(xi[i]-a)**2
t1_3/=hi1

for i in range(n):
    t2_3+=(xi[i]-a)**2
t2_3/=hi2
text.write("дов интервал для Д(изв М) для наблюдаемой св\n"+"t1_3= "+str(t1_3)+"\n")
text.write("t2_3= "+str(t2_3)+"\n")

#дов интервал для Д(неизв М) для наблюдаемой св
hi3=((c+sqrt(2*(n-1)-1))**2)/2
hi4=((-c+sqrt(2*(n-1)-1))**2)/2
t1_4=n*dx/hi3
t2_4=n*dx/hi4
text.write("дов интервал для Д(неизв М) для наблюдаемой св\n"+"t1_4= "+str(t1_4)+"\n")
text.write("t2_4= "+str(t2_4)+"\n")

###################################
#приближенный дов интервал для Д
t1_5=n*dx/(n-1+c*sqrt(2*(n-1)))
t2_5=n*dx/(n-1-c*sqrt(2*(n-1)))
text.write("приближенный дов интервал для Д\n"+"t1_5= "+str(t1_5)+"\n")
text.write("t2_5= "+str(t2_5)+"\n")
#приближенный дов для М
t1_6=x - c*dx/sqrt(n)
t2_6=x + c*dx/sqrt(n)
text.write("приближенный дов для М\n"+"t1_6= "+str(t1_6)+"\n")
text.write("t2_6= "+str(t2_6)+"\n")
#приближенный дов для Д
M4=0
for i in range(n):
    M4+=(xi[i]-x)**4
M4/=n
t1_7=dx -c*sqrt(M4-dx**2)/sqrt(n)
t2_7=dx +c*sqrt(M4-dx**2)/sqrt(n)
text.write("приближенный дов для Д\n"+"t1_7= "+str(t1_7)+"\n")
text.write("t2_7= "+str(t2_7)+"\n")
################################################

#1.5
#привести значения статистики критерия и порога
#порог
#1-alfa=0.999, N-1=8
hi2=27.877
#1-alfa=0.999, N-1-r=10-1-2=6
hi2_un=24.322

#объединеним последний и предпоследний интервалы
for m in vk,pk:
    m[-2]+=m[-1]
    m=m[0:-1]
    text.write(str(m) +"\n")


#статистика критерия
def F(u):
    return norm.cdf(u, loc = a, scale = gamma)
def F_un(u):
    return norm.cdf(u, loc = x, scale = sqrt(dx))

p = [0]*(N-1)
p_un = [0]*(N-1)
hi_stat=0
hi_stat_unknown=0

for i in range(N-1):
     p[i]=F(min(xi)+u*(i+1))-F(min(xi)+u*i)
     p_un[i]=F_un(min(xi)+u*(i+1))-F_un(min(xi)+u*i)

text.write("veroyatnost' "+(str(p))+ "\n")
text.write("veroyatnost' un "+(str(p_un))+ "\n")


for k in range(N-1):
    hi_stat += (vk[k]-n*p[k])**2/(n*p[k])
    hi_stat_unknown += (vk[k]-n*p_un[k])**2/(n*p_un[k])

text.write(str(hi_stat)+"\n")
text.write(str(hi_stat_unknown)+"\n")

if hi_stat >= hi2:
    text.write("гипотезу H0 отвергают\n")
else:
    text.write("гипотезу H0 принимают\n")

if hi_stat_unknown >= hi2_un:
    text.write("гипотезу H0 отвергают\n")
else:
    text.write("гипотезу H0 принимают\n")