# This uses PSV only for smoothness check
F=QQ
T.<x,y,z>=PolynomialRing(F)

# Generates the relevant example here
def curve():
    v = matrix([x^2, y^2, z^2])
    
    q1 = (matrix([8/7, 3, 2])*transpose(v))[0,0]
    q2 = (matrix([7/5, 3/5, 10/7])*transpose(v))[0,0]
    q3 = (matrix([2, 1, 8/5])*transpose(v))[0,0]

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


delta, q1, q2, q3 = curve();

print("Q1: " + str(q1))
print("Q2: " + str(q2))
print("Q3: " + str(q3))
print("Delta: " + str(delta))

is_delta_smooth(delta);
is_double_cover_smooth(q1, q2, q3);

# Delta is the quartic: -57/175*x^4 - 956/175*x^2*y^2 - 66/25*y^4 - 64/35*x^2*z^2 - 178/35*y^2*z^2 - 284/245*z^4
# When plugging this into PSV, the topological type of Delta is no real points. This implies that delta tilde has no real solutions.

# Output from PSV:
# -57/175*x^4 - 956/175*x^2*y^2 - 66/25*y^4 - 64/35*x^2*z^2 - 178/35*y^2*z^2 - 284/245*z^4
# Yes, it is smooth.
# The quartic has 28 bitangets, stored in 'Bitangents', and 4 real bitangents, stored in 'RealBitangents'.
# We have found:
# 63 complex Gram matrices of rank 3. Stored in 'V'.
# 15 real Gram matrices. Stored in 'Vreal'.
# 0 psd Gram matrices. Stored in 'Vpsd'.
# 8 nsd Gram matrices. Stored in 'Vnsd'.
# The real variety of the quartic is empty.