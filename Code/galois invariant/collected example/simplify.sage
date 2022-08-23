T.<x,y,z>=PolynomialRing(QQ);

r = sqrt(5*(289 - 24*sqrt(145)))*I;
Q1 = -x^2 - y^2 - z^2;
Q3 = -110/361*x^2 - 14/19*y^2;

print(Q3(r, I, 0));
print(sqrt(Q3(r, I, 0)));


# Current Example for nicer 3 ovals:

T.<x,y,z>=PolynomialRing(QQ);

h = I;
r = sqrt(10 - 4*sqrt(6))*I;
Q1 =  -x^2 - y^2 - z^2;
Q3 =  -24*x^2 + 4*y^2 - 24*z^2;

print(Q3(r, I, 0));
print(sqrt(Q3(r, I, 0)));

root = [sqrt(10 + 4*sqrt(6))*I, -sqrt(10 + 4*sqrt(6))*I, sqrt(10 - 4*sqrt(6))*I, -sqrt(10 - 4*sqrt(6))*I];

M = Matrix([[root[0], h, sqrt(Q1(root[0], h, 0)), sqrt(Q3(root[0], h, 0))],
            [root[1], -1*h, sqrt(Q1(root[1], -1*h, 0)), sqrt(Q3(root[1], -1*h, 0))],
            [root[2], h, sqrt(Q1(root[2], h, 0)), sqrt(Q3(root[2], h, 0))],
            [root[3], -1*h, sqrt(Q1(root[3], -1*h, 0)), sqrt(Q3(root[3], -1*h, 0))]]);

print(M.determinant());