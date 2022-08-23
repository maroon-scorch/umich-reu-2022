# This code is a re-implementation of FJSVV code that checks the signature of the fibers of
# the \pi_1
import sys
import csv

F=QQ

limit = 20;
iter = 1;
DEBUG = False;

T.<x,y,z>=PolynomialRing(F);
S.<t>=PolynomialRing(F);
G.<u>=PolynomialRing(RR);
#-----------------------------------------
# Helper Functions:
#-----------------------------------------
def signature(M):
    eigen = M.eigenvalues();
    
    pos_count = 0;
    neg_count = 0;
    
    for e in eigen:
        if e > 0:
            pos_count += 1;
        elif e < 0:
            neg_count += 1;
    
    return (pos_count, neg_count);

def make_matrix(M1, M2, M3, p):
    M = M1*p^2 + 2*M2*p + M3;
    A = matrix(QQ, 1, [-1])
    return block_diagonal_matrix(M, A);

def print_signature(M1, M2, M3):
    signature_list = [];
    
    char_poly =det(M1*u^2 + 2*M2*u + M3);
    root_list = list(map(lambda x: x[0], char_poly.roots()));
    length = len(root_list);
    print(root_list);   
    
    if length == 0:
        print("No Real Roots");
        num_sig = signature(make_matrix(M1, M2, M3, 1));
        signature_list.append(num_sig);
        print(num_sig);
    else:
        for i in range(0, length):
            pt = root_list[i] + 1 if i == length - 1 else (root_list[i] + root_list[i+1])/2
            print("Evaluated at " + str(pt));
            num_sig = signature(make_matrix(M1, M2, M3, pt));
            signature_list.append(num_sig);
            print(num_sig);
            
    return signature_list;

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

def conic_to_matrix(qt):
    a = Rational(qt.coefficient(x^2))
    b = Rational(qt.coefficient(y^2))
    c = Rational(qt.coefficient(z^2))
    d = Rational(qt.coefficient(x*y))
    e = Rational(qt.coefficient(y*z))
    f = Rational(qt.coefficient(x*z))
    
    M = Matrix([
    [a, d/2, f/2],
    [d/2, b, e/2],
    [f/2, e/2, c]
    ]);
    
    return M;
    
#-----------------------------------------
# Main Body of the Code:
#-----------------------------------------

def main():
    q1 = -2*y^2 - y*z - 2*z^2;
    q2 = 3*x^2 - 10*y^2 - y*z + 7*z^2;
    q3 = 9*x^2 + 4*y^2 - 8*y*z + 8*z^2;
    
    delta = q1*q3 - q2^2;
    
    M1 = conic_to_matrix(q1);
    M2 = conic_to_matrix(q2);
    M3 = conic_to_matrix(q3);
    
    if is_delta_smooth(delta) and is_double_cover_smooth(q1, q2, q3):
        try:
            signature_list = print_signature(M1, M2, M3);
        except Exception as e:
            print(e)
    print("------------------------------------------")

if __name__ == '__main__':
    main()