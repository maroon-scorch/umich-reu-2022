# This code used FJSVV and PSV
import sys
import csv

F=QQ

limit = 5;
iter = 50;
DEBUG = False;

T.<x,y,z>=PolynomialRing(F)
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

def symmetric():
    a_00 = randint(1,limit)/random_prime(limit)
    a_01 = randint(1,limit)/random_prime(limit)
    a_02 = randint(1,limit)/random_prime(limit)
    a_11 = randint(1,limit)/random_prime(limit)
    a_12 = randint(1,limit)/random_prime(limit)
    a_22 = randint(1,limit)/random_prime(limit)

    result = matrix([[a_00, 0, 0], [0, a_11, 0], [0, 0, a_22]])

    return result

# Generates a homogenous polynomial of degree 2 in 3 variables:
def generate_q1():
    # Coeff=matrix([(randint(1,limit)/random_prime(limit)) for i in [0..2]])
    Coeff = symmetric()
    v = matrix([x, y, z])
    # f = (Coeff*transpose(v))[0,0]
    f = (v*Coeff*transpose(v))[0,0]
    return f, Coeff;

def generate_q2():
    return generate_q1()

def generate_q3():
    return generate_q1()

def curve():
    q1, M1 = generate_q1()
    q2, M2 = generate_q2()
    q3, M3 = generate_q3()

    delta = q1*q3 - q2^2;
    
    return delta, q1, q2, q3, M1, M2, M3;

# This function originally used PSV's code, but has since been replaced with a blank function
def classify(f):
    classification = "N/A"
    return classification

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



#-----------------------------------------
# Main Body of the Code:
#-----------------------------------------

def main():
    if DEBUG:
        sys.stdout = open('output.txt','wt')
    
    csv_file = open('output.csv', 'w')
    writer = csv.writer(csv_file)

    header = ['Q1', 'Q2', 'Q3', 'Delta', 'Topological Type']
    writer.writerow(header)

    for i in range(iter):
        #INPUT:
        print("Iteration " + str(i + 1) + ":")
        delta, q1, q2, q3, M1, M2, M3 = curve()
        
        print("Q1 := " + str(q1) + ";")
        print("Q2 := " + str(q2) + ";")
        print("Q3 := " + str(q3) + ";")
        print("Delta: " + str(delta))

        if is_delta_smooth(delta) and is_double_cover_smooth(q1, q2, q3):
            try:
                signature_list = print_signature(M1, M2, M3);
                if signature_list.count((0, 4)) == 0:
                    print("\pi_1 is surjective! ---------------------------------");
                    classification = classify(delta)
                    data = [str(q1), str(q2), str(q3), str(delta), classification]
                    writer.writerow(data)
            except Exception as e:
                print(e)
        print("------------------------------------------")

    csv_file.close()

if __name__ == '__main__':
    main()