loadPackage "Resultants";
R = QQ[s];
S = R[x, y, z];
T = R[t];

-- Helper Function --------------------------------------------------------------
-- Given a quadratic form of 3 variables, returns its associated Symmetric Matrix
Symmetric = (Q) -> (
(M,C) := coefficients(Q, Monomials=>{x^2, y^2, z^2, x*y, x*z, y*z});
matrix {{C_(0, 0), C_(3, 0)/2, C_(5, 0)/2},
{C_(3, 0)/2, C_(1, 0), C_(4, 0)/2},
{C_(5, 0)/2, C_(4,0)/2, C_(2, 0)}});

-- Input ------------------------------------------------------------------------
Q1 = s*x^2 + y^2 - z^2;
Q2 = -43/57*x^2 - 93/14*y^2 + 85/39*z^2;
Q3 = 8/57*x^2 - 221/14*y^2 + 50/13*z^2;

-- Q1 = s*x^2 - 5*y^2 + 4*x*z - 5*y*z - 2*z^2;
-- Q2 = x^2 + y^2 + 10*x*z + 2*y*z + 3*z^2;
-- Q3 = -x^2 - y^2 - 10*x*z - 2*y*z - 4*z^2;

-- Main Body of Code ------------------------------------------------------------
Delta = Q1*Q3 - Q2^2;
<< "Delta: " << Delta

-- Constructs the Gamma Curve
M1 = Symmetric(Q1);
M2 = Symmetric(Q2);
M3 = Symmetric(Q3);

M = M1*x^2 + 2*M2*x + M3;
Gammabranch = -1*(det M);
G = sub(Gammabranch, {x => t, y => 0, z => 0});
<< "Gamma: " << G

-- Finds the discriminant and its roots
Disc = affineDiscriminant G;
<< "Discriminant Polynomial: " << Disc
<< "Roots of the Discriminant Polynomial: " << roots Disc << "."

-- Computing the partial derivative of Delta_s with respect to x, y, z
partial_disc = diff(vars S, Delta);
<< "Partial Derivative of Delta_s (to x): " << partial_disc_(0, 0)
<< "Partial Derivative of Delta_s (to y): " << partial_disc_(0, 1)
<< "Partial Derivative of Delta_s (to z): " << partial_disc_(0, 2)

-- Computing Res_3(\parital_x \Delta_s, \parital_y \Delta_s, \parital_s \Delta_s)
F = {partial_disc_(0, 0), partial_disc_(0, 1), partial_disc_(0, 2)};
P = resultant F;

<< "Resultant: " << F

-- Checks if the Discriminant divides the Resultant:
rem = P % Disc;
<< "Remainder of Resultant Divided by Discriminant: " << rem

-- Find the Root of the Resultant
<< "Roots of the Resultant: " << roots P << "."