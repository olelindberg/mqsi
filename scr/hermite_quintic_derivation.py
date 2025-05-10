from sympy import Symbol
from sympy import *

t = Symbol('t')

H   = Matrix([
1 -                 10*t**3 +  15*t**4 -   6*t**5,     
    t -              6*t**3 +   8*t**4 -   3*t**5,       
        0.5*t**2 - 1.5*t**3 + 1.5*t**4 - 0.5*t**5, 
                    10*t**3 -  15*t**4 +   6*t**5,       
                 -   4*t**3 +   7*t**4 -   3*t**5,        
                   0.5*t**3 -     t**4 + 0.5*t**5])


Ht   = diff(H,t)
Htt  = diff(Ht,t)
Httt = diff(Htt,t)


print("H")
print(H)

print("Ht")
print(Ht)

print("Htt")
print(Htt)

print("Httt")
print(Httt)
