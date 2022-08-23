# This uses PSV only for smoothness check
import sys
import csv

F=QQ

limit = 5;
iter = 1;
DEBUG = False;

T.<x,y,z>=PolynomialRing(F)

#-----------------------------------------
# Helper Functions:
#-----------------------------------------
def symmetric():
    a_00 = randint(1,limit)/random_prime(limit)
    a_01 = randint(1,limit)/random_prime(limit)
    a_02 = randint(1,limit)/random_prime(limit)
    a_11 = randint(1,limit)/random_prime(limit)
    a_12 = randint(1,limit)/random_prime(limit)
    a_22 = randint(1,limit)/random_prime(limit)

    result = matrix([[a_00, a_01, a_02], [a_01, a_11, a_12], [a_02, a_12, a_22]])

    return result

# Generates a homogenous polynomial of degree 2 in 3 variables:
def generate_q1():
    # Coeff=matrix([(randint(1,limit)/random_prime(limit)) for i in [0..2]])
    Coeff = symmetric()
    v = matrix([x, y, z])
    # f = (Coeff*transpose(v))[0,0]
    f = (v*Coeff*transpose(v))[0,0]
    return f

def generate_q2():
    return generate_q1()

def generate_q3():
    return generate_q1()

def curve():
    q1 = generate_q1()
    q2 = generate_q2()
    q3 = generate_q3()

    delta = q2^2 - q1*q3
    
    return delta, q1, q2, q3;

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

# Checks if the Double Cover has any real solutions
def has_no_solution(q1, q2, q3):
    R.<x,y,z,r,s> = ProjectiveSpace(F, 4);

    q1_new = q1 + r - r;
    q2_new = q2 + r - r;
    q3_new = q3 + r - r;
    
    q1_t = r^2;
    q2_t = r*s;
    q3_t = s^2;

    C = R.curve([q1_new, q2_new, q3_new]);
    D = R.curve([q1_t, q2_t, q3_t]);

    V = C.intersection(D);
    solutions = V.irreducible_components();
    print("Solutions: ");
    print(solutions);

    return solutions == [];
#-----------------------------------------
# Main Body of the Code:
#-----------------------------------------

def main():
    if DEBUG:
        sys.stdout = open('output.txt','wt')

    for i in range(iter):
        #INPUT:
        print("Iteration " + str(i + 1) + ":")
        delta, q1, q2, q3 = curve()
        
        print("Q1: " + str(q1))
        print("Q2: " + str(q2))
        print("Q3: " + str(q3))
        print("Delta: " + str(delta))

        if is_delta_smooth(delta) and is_double_cover_smooth(q1, q2, q3):
            try:
                if has_no_solution(q1, q2, q3):
                    print("The Double Cover has no real solutions");
            except Exception as e:
                print(e)
        print("------------------------------------------")

if __name__ == '__main__':
    main()