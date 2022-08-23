import numpy as np
F = QQ
T.<x,y,z>=PolynomialRing(F)

bound = 0.000001;

Q1 = -x^2 - y^2 - 3*z^2
Q2 = 3*x^2 + 5*y^2
Q3 = -7*x^2 - 23*y^2 - 12*z^2

Delta = Q1*Q3 - Q2^2;

M = [
[0.36065792 + 0.95391314*I, 0.31967104 - 0.47695657*I, 1.0000000, -0.13181154 + 1.4533410*I, 0.55179280 + 1.9908512*I],
[0.36065792 - 0.95391314*I, 0.31967104 + 0.47695657*I, 1.0000000, -0.13181154 - 1.4533410*I, 0.55179280 - 1.9908512*I],
[-0.71437360 + 0.73499412*I, -0.28562640 - 0.73499412*I, 1.0000000, -0.19732388 - 1.5970022*I, 0.82981304 - 1.3895333*I],
[-0.71437360 - 0.73499412*I, -0.28562640 + 0.73499412*I, 1.0000000, -0.19732388 + 1.5970022*I, 0.82981304 + 1.3895333*I]
];


for i in range(4):
    a = M[i][0]
    b = M[i][1]
    c = M[i][2]
    r = M[i][3]
    s = M[i][4]
    
    assert abs(Delta(a, b, c)) < bound;
    assert abs(Q1(a, b, c) - r^2) < bound;
    assert abs(Q3(a, b, c) - s^2) < bound;
    assert abs(Q2(a, b, c) - r*s) < bound;

print("-----------------");
mat = np.matrix(M, dtype=complex);
print(Matrix(M));
mat_rk = np.linalg.matrix_rank(mat);
print("Rank of Matrix: " + str(mat_rk));

print("------------------");

for i in range(5):
    new_mat = np.delete(mat, i, 1);
    new_det = np.linalg.det(new_mat);
    print("Determinant with Column " + str(i) + " Removed: " + str(new_det));
    print("Close to Zero: " + str(abs(new_det) < bound));
    print("-------------------------------");