x, y = var('x, y')



T1 = (-2.140100881975361)*x + (-0.5740020094011272)*y + 1;

Q1 = x^2 - 2*x*y - y^2 - 2*x*1 + 2*y*1 - 1^2;
Q2 = -17/2*x^2 + 16*x*y + 7*y^2 + 15*x*1 - 15*y*1 - 17/2*1^2;
Q3 = -60*x^2 + 34*x*y - 33*y^2 - 108*y*1 - 52*1^2;

Delta = Q1*Q3 - Q2^2;

roots = solve([Delta == 0, T1 == 0], x, y)
# print(roots);
roots.pop(3);
roots.pop(1);
print("Bitangents intersect on z = 1");
print(roots);
print("-------------------------");

gal_matrix = [];

for r in roots:
    c1 = r[0] + I - I;
    c2 = r[1] + I - I;
    c3 = sqrt(Q1.subs(c1, c2))
    c4 = sqrt(Q3.subs(c1, c2))
    row1 = [c1.rhs(), c2.rhs(), 1, c3, c4]
    row2 = [c1.rhs(), c2.rhs(), 1, -c3, -c4]
    gal_matrix.append(row1);
    gal_matrix.append(row2);

gal_matrix = Matrix(gal_matrix);
print("Matrix: \n" + str(gal_matrix));
print("Rank of Matrix: " + str(gal_matrix.rank()));

