import cmath
import math
#argument_radians = math.radians(argument_degrees)

p0 = complex(-6.3,0.4)
p1= complex(-8.8,0.4)
p2 = complex(-6.4,2.1)

p0c = p0.conjugate()
p1c= p1.conjugate()
p2c = p2.conjugate()

t1 = 35
t2 = 150
t3 = 268

a1 = p1*cmath.rect(1,t2-t1)-p0
a1c = a1.conjugate()
a2 = p2*cmath.rect(1,t3-t1)-p0
a2c = a2.conjugate()


l2_0 = (p0*p0c*(a1-a2)+a2*p1*p1c-a1*p2*p2c)/(a1c*a2-a2c*a1)
print(l2_0)
print(abs(l2_0))

