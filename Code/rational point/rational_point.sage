# This example always has a rational point
F=QQ

T.<x,y,z>=PolynomialRing(F)

#INPUT:
limit = 100;

a = randint(-1*limit,limit);
b = randint(-1*limit,limit);
c = randint(-1*limit,limit);
e = randint(-1*limit,limit);
f = randint(-1*limit,limit);
g = randint(-1*limit,limit);

Q1 = 2*x^2 - 9*y^2 + e*x*z + a*z^2;
Q2 = x^2 - 6*y^2 + f*y*z + b*z^2;
Q3 = x^2 - 8*y^2 + g*x*z + c*z^2;

print(Q1);
print(Q2);
print(Q3);

f= Q1*Q3 - Q2^2;

print(f);
