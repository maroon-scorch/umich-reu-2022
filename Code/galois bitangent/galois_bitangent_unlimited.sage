# This code used PSV
#INPUT: Q1, Q2, Q3, appropriate plane conics with Rational Coefficients
#OUTPUT: The 5 x 4 matrices formed by points of DeltaTilde over  the intersection
# of the real bitangents of Delta := Q1*Q3 - Q2^2 with Delta. Checks whether or not they are full-rank.
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

# Given Q1, Q2, Q3, the bitangent line, and a z-value
# Produces the 5x4 matrix with the 4 points on DeltaTilde given by the bitangent line and Delta
# The matrix is checked to be full-rank or not
# Returns True if processed successfully, False otherwise
def galois_matrix(q1, q2, q3, bitangent, z_val):
    x, y = var('x, y')
    Q1 = q1 + x - x;
    Q2 = q2 + x - x;
    Q3 = q3 + x - x;
    T = bitangent + x - x;
    
    Delta = Q1*Q3 - Q2^2;

    print("Intersection on Delta and Bitangent on z = 1");

    # Remove duplicate roots resulted from imprecision of Sage's solve function
    roots = solve([Delta == 0, T == 0], x, y)
    roots = prune_similar(roots);
    print(roots);
    
    if len(roots) != 2:
        print("Something went wrong in the computation, the list should have 2 roots. Try switching another chart!");
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
            row1 = [c1, c2, z_val, c3, c4]
            row2 = [c1, c2, z_val, -c3, -c4]
        else:
            row1 = [c1, c2, z_val, -c3, c4]
            row2 = [c1, c2, z_val, c3, -c4]

        gal_matrix.append(row1);
        gal_matrix.append(row2);

    gal_mat = Matrix(gal_matrix);
    print("Matrix: \n" + str(gal_mat));
    process_matrix(gal_matrix);
    
    return True;

##############################################################

# The following Sage code that computes the real bitangents of the quartic curve Delta was written by Daniel Plaumann, Bernd Sturmfels, Cynthia Vinzant, and is from the code accompanying their paper "Quartic curves and their bitangents."
def real_bitangents(f):
    RealBitangents = []
    return RealBitangents

##############################################################

# Next, using the real bitangents computed from the code of Plaumann, Sturmfels, and Vinzant, we find the four points of Deltatilde lying above the intersection of each bitangent with Delta.

#-----------------------------------------
# Main Body of the Code:
#-----------------------------------------

def main():
    q1 = -x^2 - y^2 - 3*z^2
    q2 = 3*x^2 + 5*y^2
    q3 = -7*x^2 - 23*y^2 - 12*z^2
    
    delta = q1*q3 - q2^2;
    
    z_val = 1;
    
    print("Q1 := " + str(q1) + ";")
    print("Q2 := " + str(q2) + ";")
    print("Q3 := " + str(q3) + ";")
    print("Delta: " + str(delta))
    
    
    if is_delta_smooth(delta) and is_double_cover_smooth(q1, q2, q3):
        try:
            real_bts = real_bitangents(delta)
        except Exception as e:
            print(e)
    print("------------------------------------------")
    
    i = 1;
    for bt in real_bts:
        print("Iteration " + str(i));
        print("Bitangent: " + str(bt));
        
        gq1 = q1(x, y, z_val)
        gq2 = q2(x, y, z_val)
        gq3 = q3(x, y, z_val)
        gbt = bt(x, y, z_val, 0, 0, 0, 0, 0, 0, 0)
        
        galois_matrix(gq1, gq2, gq3, gbt, z_val);
        i = i + 1;
        print("------------------------------------------")
        
if __name__ == '__main__':
    main()