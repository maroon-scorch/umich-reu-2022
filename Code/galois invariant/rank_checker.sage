F=QQ;
T.<x,y,z>=PolynomialRing(F);
R.<u>=PolynomialRing(QQbar);

Q1 = 7/17*x^2 + 4/13*x*y + 7/13*y^2 - 4/19*x*z + 9*y*z + 5/13*z^2;
Q2 = 10/17*x^2 + 13/19*x*y + 13/34*y^2 + 19/13*x*z + 17/3*y*z + 7/19*z^2;
Q3 = -8*x^2 - 8*x*y - 10/3*y^2 + 8/17*x*z - 12*y*z - 4*z^2;

Delta = Q1*Q3 - Q2^2;

h = 1;
quartic = Delta(u, h, 0);
root = list(map(lambda x: x[0], quartic.roots()));

M = Matrix([[root[0], h, sqrt(Q1(root[0], h, 0)), sqrt(Q3(root[0], h, 0))],
            [root[1], h, sqrt(Q1(root[1], h, 0)), sqrt(Q3(root[1], h, 0))],
            [root[2], h, sqrt(Q1(root[2], h, 0)), sqrt(Q3(root[2], h, 0))],
            [root[3], h, sqrt(Q1(root[3], h, 0)), sqrt(Q3(root[3], h, 0))]]);

print(M)
M.determinant()