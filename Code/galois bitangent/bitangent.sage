# This code used PSV only to check smoothness
import sys
import csv

F=QQ

limit = 5;
iter = 50;
threshold = 0.0001;

T.<x,y,z>=PolynomialRing(F)

#-----------------------------------------
# Helper Functions:
#-----------------------------------------
def difference(r, s):
    x1 = r[0].rhs();
    x2 = r[1].rhs();
    
    y1 = s[0].rhs();
    y2 = s[1].rhs();
    
    diff = (x1 - y1).abs()^2 + (x2 - y2).abs()^2
    
    return sqrt(diff);

def is_similar_in_list(elt, lst):
    for item in lst:
        if difference(elt, item) < threshold:
            return True
    return False

def prune_similar(lst):
    unique = [];
    
    for elt in lst:
        if not is_similar_in_list(elt, unique):
            unique.append(elt);
    return unique;

def real_bitangents(f):
    RealBitangents = []
    return RealBitangents

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

def galois_matrix(f, bitangent, z_val):
    x, y = var('x, y')
    Delta = f + x - x;
    T = bitangent + x - x;

    print("Intersection on Delta and Bitangent on z = 1");

    roots = solve([Delta == 0, T == 0], x, y)
    # print(roots);
    roots = prune_similar(roots);
    print(roots);
    
    if len(roots) != 2:
        print("Something went wrong in the computation, the list should have 2 roots");
        return;
    print("-------------------------");

#-----------------------------------------
# Main Body of the Code:
#-----------------------------------------

def main():
    f = 10*x^4+15*x^3*y-17*x^2*y^2+15*x*y^3+10*y^4+15*x^3*z-833*x^2*y*z-833*x*y^2*z+15*y^3*z-17*x^2*z^2-833*x*y*z^2-17*y^2*z^2 + 15*x*z^3 + 15*y*z^3 + 10*z^4;
    
    z_val = 1;

    print("Delta: " + str(f))
    
    if is_delta_smooth(f):
        try:
            real_bts = real_bitangents(f)
        except Exception as e:
            print(e)
    print("------------------------------------------")
    
    i = 1;
    for bt in real_bts:
        print("Iteration " + str(i));
        print("Bitangent: " + str(bt));
        
        gf = f(x, y, z_val);
        gbt = bt(x, y, z_val, 0, 0, 0, 0, 0, 0, 0)
        
        galois_matrix(gf, gbt, z_val);
        i = i + 1;
        print("------------------------------------------")
        
if __name__ == '__main__':
    main()