# This code used PSV
import sys
import csv

F=QQ

# If I set limit to 20, there's a lot of 1 ovals
limit = 10;
iter = 25;
DEBUG = False;

T.<x,y,z>=PolynomialRing(F)

#-----------------------------------------
# Helper Functions:
#-----------------------------------------
def symmetric():
    a_00 = randint(1,limit)/random_prime(limit)
    a_11 = randint(1,limit)/random_prime(limit)
    a_22 = randint(1,limit)/random_prime(limit)

    # result = matrix([[a_00, a_01, a_02], [a_01, a_11, a_12], [a_02, a_12, a_22]])
    result = matrix([[a_00, 0, 0], [0, a_11, 0], [0, 0, a_22]])

    return result

# Generates a homogenous polynomial of degree 2 in 3 variables:
def generate_q1():
    Coeff = symmetric()
    v = matrix([x, y, z])
    f = (v*Coeff*transpose(v))[0,0]
    return f

def generate_q2():
    return generate_q1()

def generate_q3():
    return generate_q1()

def curve():
    q1 = -1*generate_q1()
    q2 = generate_q2()
    q3 = generate_q3()

    delta = q2^2 - q1*q3
    
    return delta, q1, q2, q3;

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
        delta, q1, q2, q3 = curve()
        
        print("Q1 := " + str(q1) + ";")
        print("Q2 := " + str(q2) + ";")
        print("Q3 := " + str(q3) + ";")
        print("Delta: " + str(delta))

        if is_delta_smooth(delta) and is_double_cover_smooth(q1, q2, q3):
            try:
                classification = classify(delta)
                data = [str(q1), str(q2), str(q3), str(delta), classification]
                writer.writerow(data)
            except Exception as e:
                print(e)
        print("------------------------------------------")

    csv_file.close()

if __name__ == '__main__':
    main()