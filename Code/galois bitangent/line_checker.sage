#INPUT: Q1, Q2, Q3, appropriate plane conics with Rational Coefficients, and a list of homogeneous lines with Rational Coefficients
#OUTPUT: The 5 x 4 matrices formed by points of DeltaTilde over the intersection
# of Delta := Q1*Q3 - Q2^2 and each line in the list given. Checks whether or not they are full-rank.
import numpy as np

F=QQ

# Bound on Similar Roots
bound = 0.001;

T.<x,y,z>=PolynomialRing(F)

#-----------------------------------------
# Helper Functions:
#-----------------------------------------
# Given two expressions (x == x_1, y == x_2), (x == y_1, y == y_2), computes the
# distance between the points they represent
def difference(r, s):
    x1 = r[0].rhs();
    x2 = r[1].rhs();
    
    y1 = s[0].rhs();
    y2 = s[1].rhs();
    
    diff = (x1 - y1).abs()^2 + (x2 - y2).abs()^2
    
    return sqrt(diff);

# Given an expression and a list of expressions, checks if the expressions
# is close to any item of the list
def is_similar_in_list(elt, lst):
    for item in lst:
        if difference(elt, item) < bound:
            return True
    return False

# Due to Sage's implementations, using solve may produce multiple solutions that are close to
# to each other but are actually the same solution, this function is here to prune out these duplicate roots.
def prune_similar(lst):
    unique = [];
    
    for elt in lst:
        if not is_similar_in_list(elt, unique):
            unique.append(elt);
    return unique;

#check whether the discriminant curve is non-singular
def is_delta_smooth(f):
    Grad=ideal(f,diff(f,x),diff(f,y),diff(f,z))
    if Grad.dimension()==0:
        print("The discriminant curve is smooth.")
        return True
    else:
        print("The discriminant curve is not smooth.")
        return False

#check whether the double cover (delta tilde) is non-singular
def is_double_cover_smooth(q1, q2, q3):
    R.<x,y,z,r,s> = ProjectiveSpace(F, 4);

    q1_t = q1 - r^2;
    q2_t = q2 - r*s;
    q3_t = q3 - s^2;
    
    delta_tilde = R.curve([q1_t, q2_t, q3_t]);

    if delta_tilde.is_singular():
        print("The double cover is not smooth.");
        return False;
    else:
        print("The double cover is also smooth.");
        return True;

# Takes in the standard 5x4 Gal(C/R)-matrix (Double Array) and checks whether it's full rank
def process_matrix(M):
    gal_matrix = np.matrix(M, dtype=complex);
    gal_rk = np.linalg.matrix_rank(gal_matrix);
    print("Rank of Matrix: " + str(gal_rk));
    
    # Note that M is full-rank if and only if at least one of its 4 x 4 submatrix has
    # a non-zero determinant
    for i in range(5):
        new_mat = np.delete(gal_matrix, i, 1);
        new_det = np.linalg.det(new_mat);
        print("Determinant with Column " + str(i) + " Removed: " + str(new_det));
        # print("Close to Zero: " + str(abs(new_det) < bound));
        # print("--------");

# Given Q1, Q2, Q3, and a line, and a z-value
# Produces the 5x4 matrix with the 4 points (typically) on DeltaTilde given
# by the 4 intersections (typically) of Delta and the line.
# The matrix is checked to be full-rank or not
# Returns True if processed successfully, False otherwise
def galois_matrix(q1, q2, q3, line, z_val):
    x, y = var('x, y')
    Q1 = q1 + x - x;
    Q2 = q2 + x - x;
    Q3 = q3 + x - x;
    T = line + x - x;
    
    Delta = Q1*Q3 - Q2^2;

    print("Intersection on Delta and Line on z = 1");

    # Remove duplicate roots resulted from imprecision of Sage's solve function
    roots = solve([Delta == 0, T == 0], x, y)
    roots = prune_similar(roots);
    print(roots);
    
    if len(roots) < 4:
        print("The line does not have at least 4 intersections with Delta");
        return False;
    print("-------------------------");

    gal_matrix = [];
    
    # For each root in the list, construct the appropriate row in the matrix and add it in
    for r in roots:
        c1 = (r[0] + I - I).rhs();
        c2 = (r[1] + I - I).rhs();
        
        if abs(c1.imag()) < 0.0001:
            c1 = c1.real()
        
        if abs(c2.imag()) < 0.0001:
            c2 = c2.real()
        
        c1 = c1.n(50) # This is x
        c2 = c2.n(50) # This is y
        c3 = sqrt(Q1.subs(x == c1, y == c2).n(50)) # This is r
        c4 = sqrt(Q3.subs(x == c1, y == c2).n(50)) # This is s
        
        q2_val = Q2.subs(x == c1, y == c2);
        
        print("Q2 : " + str(q2_val));
        
        # Adjusting r and s based on Q2(x, y, z)
        if abs(q2_val - c3*c4) < 0.00001:
            row = [c1, c2, z_val, c3, c4]
        else:
            row = [c1, c2, z_val, -c3, c4]

        gal_matrix.append(row);

    gal_mat = Matrix(gal_matrix);
    print("Matrix: \n" + str(gal_mat));
    process_matrix(gal_matrix);
    
    return True;

#-----------------------------------------
# Main Body of the Code:
#-----------------------------------------

def main():
    # inputs:
    q1 = 2*x^2 - 9*y^2 + x*z - 82*z^2;
    q2 = x^2 - 6*y^2 + 74*y*z - 49*z^2;
    q3 = x^2 - 8*y^2 + 37*x*z + 79*z^2;
    lines = [x + y - 2*z, -2/5*x + y + 1/10*z];
    
    delta = q1*q3 - q2^2;
    z_val = 1;
    
    print("Q1 := " + str(q1) + ";")
    print("Q2 := " + str(q2) + ";")
    print("Q3 := " + str(q3) + ";")
    print("Delta: " + str(delta))
    print("------------------------------------------")
    
    if is_delta_smooth(delta) and is_double_cover_smooth(q1, q2, q3):
        i = 1;
        for line in lines:
            print("Iteration " + str(i));
            print("Line: " + str(line));
            
            gq1 = q1(x, y, z_val)
            gq2 = q2(x, y, z_val)
            gq3 = q3(x, y, z_val)
            gline = line(x, y, z_val)
            
            galois_matrix(gq1, gq2, gq3, gline, z_val);
            i = i + 1;
            print("------------------------------------------")
    

        
if __name__ == '__main__':
    main()