# This smoothness check is an attempt to check if DeltaTilde is smooth and is not related to PSV
F=QQ
T.<x,y,z,r,s>=PolynomialRing(F)

def is_smooth(q_1, q_2, q_3):
    q1_tilde = q1 - r^2
    q2_tilde = q2 - r*s
    q3_tilde = q3 - s^2

    delta_tilda = (q1_tilde, q2_tilde, q3_tilde);

    Grad=ideal(delta_tilda,diff(delta_tilda,x),diff(delta_tilda,y),diff(delta_tilda,z), diff(delta_tilda, r), diff(delta_tilda, s));
    if Grad.dimension()==0:
        print("Yes, it is smooth.")
    else:
        print("No, it is not smooth.")


P.<x,y,z,r,s> = ProjectiveSpace(QQ, 4);

q1 = x^2 + 2*x*y + 2/3*y^2 + 8/3*x*z + 2*y*z + 2/3*z^2;
q2 = x^2 + 4*x*y + 5/3*y^2 + 5*x*z + y*z + 3/2*z^2;
q3 = 1/2*x^2 + 5*x*y + 5/2*y^2 + 2/3*x*z + 4*y*z + 1/2*z^2;

q1_t = q1 - r^2;
q2_t = q2 - r*s;
q3_t = q3 - s^2;

C = P.curve([q1_t, q2_t, q3_t]);

C.is_singular();

# R.<x,y,z> = ProjectiveSpace(QQ, 2);
# Delta = R.curve([q2*q2 - q1*q3]);
# Delta.is_singular();
