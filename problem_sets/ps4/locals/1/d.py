from sympy import *
x1,x2,x3,x4,x5,x6,x7 = symbols('x1 x2 x3 x4 x5 x6 x7') 
print(
    integrate((x1+x2+x3+x4+x5+x6+x7)*(x1+x2+x3+x4+x5+x6+x7),
          (x1,0,1),(x2,0,1),(x3,0,1),(x4,0,1),(x5,0,1),(x6,0,1),(x7,0,1)
        )
)
