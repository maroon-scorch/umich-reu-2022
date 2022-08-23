F=QQ

T.<x,y,z>=PolynomialRing(F)
R.<u, v> = PolynomialRing(F)

#INPUT:
a = 13;
b = -3;
c = 6;

Q1 = 2*x^2 - 9*y^2 + a*z^2;
Q2 = x^2 - 6*y^2 + b*z^2;
Q3 = x^2 - 8*y^2 + c*z^2;

Delta = Q1*Q3 - Q2^2;

print(Delta(3, 1, 0));
print(Q1(3, 1, 0));
print(Q3(3, 1, 0));
print("------------------");
print(Q2(3, 1, 0));
print(sqrt(Q1(3, 1, 0)*Q3(3, 1, 0)))