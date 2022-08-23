# This code used FJSVV
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

Q1 = 2*x^2 + 2/3*x*y + 2/3*y^2 + 2/3*x*z + 4/3*y*z + 1/2*z^2;
Q2 = 4/3*x^2 + 4*x*y + 3/2*y^2 + 3*x*z + 2*y*z + 1/2*z^2;
Q3 = 4/3*x^2 + x*y + 2/3*y^2 + 5*x*z + 10/3*y*z + 5/3*z^2;
delta = Q1*Q3 - Q2^2;

if is_delta_smooth(delta) and is_double_cover_smooth(Q1, Q2, Q3):
    signature_list = print_signature(M1, M2, M3);