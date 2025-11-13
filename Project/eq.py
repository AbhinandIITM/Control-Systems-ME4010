import sympy as sp

# Define symbols
x, x_dot, theta, theta_dot, u = sp.symbols('x x_dot theta theta_dot u')
Ir, Mb, Mr, Ib, h, r, d, g = sp.symbols('Ir Mb Mr Ib h r d g')

# --- Mass matrix terms ---
M11 = Ir + Mb*(h**2 + 2*h*r + r**2)
M12 = -Mb*(h + r)
M21 = -Mb*(h + r)
M22 = Ib/r**2 + Mb

# Mass matrix
M = sp.Matrix([[M11, M12],
               [M21, M22]])

# --- Gravitational and coupling terms (linearized) ---
# For small theta, cos(theta) ~ 1 and sin(theta) ~ theta
G1 = (Mr*d*g + Mb*g*h + Mb*g*r)*theta  # torque-like restoring term
G2 = -g*Mb*(h + r)*theta               # force-like coupling

# --- Input matrix ---
B = sp.Matrix([[0],
               [Ib/r**2 * u]])

# --- Acceleration equations: M * [theta_ddot, x_ddot]^T = -G + B ---
rhs = -sp.Matrix([[G1],
                  [G2]]) + B

# Solve for [theta_ddot, x_ddot]
acc = M.inv() * rhs
theta_ddot = sp.simplify(acc[0])
x_ddot = sp.simplify(acc[1])

# --- Define state vector ---
# x1 = x, x2 = x_dot, x3 = theta, x4 = theta_dot
x1, x2, x3, x4 = x, x_dot, theta, theta_dot

# State derivatives
dx1 = x2
dx2 = x_ddot.subs(theta, x3)
dx3 = x4
dx4 = theta_ddot.subs(theta, x3)

# State vector and derivatives
X = sp.Matrix([x1, x2, x3, x4])
dX = sp.Matrix([dx1, dx2, dx3, dx4])

# --- Linearize: take Jacobians w.r.t. states and input ---
A = dX.jacobian(X)
B_mat = dX.jacobian(sp.Matrix([u]))

# Simplify the matrices
A_simplified = sp.simplify(A)
B_simplified = sp.simplify(B_mat)

sp.pprint(A_simplified)
print("\nB matrix:")
sp.pprint(B_simplified)
print(sp.latex(A_simplified))
print(sp.latex(B_simplified))



from sympy import symbols, Matrix, Identity, pprint, simplify

N_x = [m_b, m_b*(h+r), (I_w/(R**2)) + M_r + M_w + m_b]
N_theta = [h*m_b, m_b*((h+r)**2) + I_r, m_b*(h+r)]
N_x_w = [(I_b/(r**2)) + m_b, m_b*(h+r), m_b]

N = Matrix([N_x, N_theta, N_x_w])

# 3. Define the I matrix
I_x = [0, 0, 0]
I_theta = [0, -M_r*d*g - m_b*(h+r)*g, 0]
I_x_w = [0, -m_b*g, 0]

I = Matrix([I_x, I_theta, I_x_w])

# 4. Define the Identity matrix
Id = Identity(3)

# 5. Perform the calculation: M = (Id - I) * N.inv()

# First, calculate the numerator: (Id - I)
Numerator = Id - I

# Next, calculate the matrix inverse of N
N_inverse = N.inv()

# Finally, perform matrix multiplication
M = Matrix(Numerator * N_inverse)

# 6. Print the resulting matrix
print("Resulting Matrix M = (Id - I) * N.inv()")

# Note: The symbolic inverse is extremely large.
# We will use simplify() to make the output
# as clean as possible, but it will still be complex.
M_inv = simplify(M.inv())

from sympy import symbols, Matrix, Identity, pprint, simplify

N_x = [m_b, m_b*(h+r)]
N_theta = [h*m_b, m_b*((h+r)**2) + I_r]


N = Matrix([N_x, N_theta])

# 3. Define the I matrix
I_x = [0, 0]
I_theta = [-g*m_b, -M_r*d*g - m_b*(h+r)*g]

M = Matrix([I_x, I_theta])

G = Matrix([[(1/R) + ((R/I_w)*(M_r + M_w +m_b)) ], [(h+r)*m_b]])

# 4. Define the Identity matrix
Id = Identity(2)

# 5. Perform the calculation: M = (Id - I) * N.inv()

# First, calculate the numerator: (Id - I)
I1 = Id - M

# Next, calculate the matrix inverse of N
N_inverse = N.inv()

# Finally, perform matrix multiplication
M1 = Matrix(N_inverse*I1)


M2 = Matrix(N_inverse*G)
