# Performs a PGL_2 action on Q1, Q2, Q3
# The resulting Y is invariant under this action
# This is Theorem 2.6 of FJSVV's paper "Curve classes on conic bundle threefolds and applications to rationality"

F=QQ;
T.<x,y,z>=PolynomialRing(F);
R.<u>=PolynomialRing(QQbar);

Q1 =  x^2 - y^2 - z^2;
Q2 =  y^2;
Q3 =  -16*x^2 + 6*y^2 + z^2;

a = 1;
b = -2;
c = 3;
d = 4;

M = Matrix([
    [b^2, 2*b*d, d^2],
    [a*b, a*d + b*c, c*d],
    [a^2, 2*a*c, c^2]     
           ])
v = Matrix([[Q1], [Q2], [Q3]]);

result = M*v;

Q1 = result[0][0];
Q2 = result[1][0];
Q3 = result[2][0];

print(Q1);
print(Q2);
print(Q3);

Delta = Q1*Q3 - Q2^2;

print(Delta);

h = 1
quartic = Delta(u, h, 0);

print(quartic);
print(Delta(u, -1*h, 0));

root = list(map(lambda x: x[0], quartic.roots()));

M = Matrix([[root[0], h, sqrt(Q1(root[0], h, 0)), sqrt(Q3(root[0], h, 0))],
            [root[1], h, sqrt(Q1(root[1], h, 0)), sqrt(Q3(root[1], h, 0))],
            [root[2], h, sqrt(Q1(root[2], h, 0)), sqrt(Q3(root[2], h, 0))],
            [root[3], h, sqrt(Q1(root[3], h, 0)), sqrt(Q3(root[3], h, 0))]]);

print(M)
M.determinant()
