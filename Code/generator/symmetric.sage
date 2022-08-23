# The smoothness check here is not related to PSV
F=QQ
T.<x,y,z,r,s>=PolynomialRing(F)

limit=10

# Generates a 3x3 symmetric matrix
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
def generate():
    # Coeff=matrix([(randint(1,limit)/random_prime(limit)) for i in [0..2]])
    Coeff = symmetric()
    v = matrix([x, y, z])
    # f = (Coeff*transpose(v))[0,0]
    f = (v*Coeff*transpose(v))[0,0]
    return f

q1 = generate()
q2 = generate()
q3 = generate()

delta = q2^2 - q1*q3

q1_tilde = q1 - r^2
q2_tilde = q2 - r*s
q3_tilde = q3 - s^2

m1 = jacobian(delta, (x,y,z,r,s))
if (m1.rank() == 1):
    print("The discriminant curve is smooth")
else:
    sys.exit("The discriminant curve is not smooth")
    
m2 = jacobian((q1_tilde, q2_tilde, q3_tilde), (x, y, z, r, s))
if (m2.rank() == 3):
    print("The discriminant double cover is smooth")
else:
    sys.exit("The discriminant double cover is not smooth")
    
print("Q1: " + str(q1))
print("Q2: " + str(q2))
print("Q3: " + str(q3))
print("Delta: " + str(delta))