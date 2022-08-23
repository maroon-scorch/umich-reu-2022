R.<x, y, z> = QQ[]

q1= x^2 + y^2 - 3*z^2;
q2= -3*x^2 - 2*x*y - 4*y^2 + 2*z^2;
q3= -4*x^2 - 5*y^2;

Delta = q1*q3 - q2^2;

tangent = (-0.36271048855667953)*x + 0.21909666531808614*y + z;

a = 2.0254111 + 0.84737361*I
b = -1.2111648 + 1.4028114*I
c = 1;
d = -0.049945115 - 0.34526257*I;
e = 0.48560133 + 3.3567845*I

print(Delta(a, b, c));
print(tangent(a, b, c));
print(q2(a, b, c));
print(d*e);