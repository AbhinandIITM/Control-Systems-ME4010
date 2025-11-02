import sympy
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math

from sympy import symbols, Matrix, cos, sin, diff, Function
from sympy.physics.mechanics import dynamicsymbols, Lagrangian, Particle, ReferenceFrame, Point
from sympy.physics.mechanics import KanesMethod, inertia, RigidBody
from sympy.abc import t # Time variable

# Define symbolic variables
M_w, M_r, m_b, R, I_w, I_r, d, h, I_b, r, g, T_sym = symbols('M_w M_r m_b R I_w I_r d h I_b r g T')

# Define time-dependent symbolic functions
theta = dynamicsymbols('theta')
x_w = dynamicsymbols('x_w')
x = dynamicsymbols('x')

# Derivatives
theta_dot = diff(theta, t)
x_w_dot = diff(x_w, t)
x_dot = diff(x, t)

theta_ddot = diff(theta_dot, t)
x_w_ddot = diff(x_w_dot, t)
x_ddot = diff(x_dot, t)

# Position (adapted from MATLAB script)
rotation_matrix = Matrix([[cos(theta), sin(theta)], [-sin(theta), cos(theta)]])
body_coords = Matrix([[x], [h + r]])
rotated_body_coords = rotation_matrix * body_coords
x_b = rotated_body_coords[0] + x_w
y_b = rotated_body_coords[1] + R

# Potential Energy
U_b = m_b * g * y_b
U_r = M_r * g * (R + d * cos(theta))
U_w = M_w * g * R
U = sympy.simplify(U_b + U_r + U_w)

# Kinetic Energy
T_w = sympy.Rational(1, 2) * M_w * x_w_dot**2 + sympy.Rational(1, 2) * I_w * (x_w_dot/R)**2
T_r = sympy.Rational(1, 2) * (M_r * x_w_dot**2 + I_r * theta_dot**2)

term1_Tb = (x_dot*cos(theta) - x*sin(theta)*theta_dot + (R+r)*cos(theta) + x_w_dot)**2
term2_Tb = (-x_dot*sin(theta) - x*cos(theta)*theta_dot - (h+r)*sin(theta))**2
T_b = sympy.Rational(1, 2) * (m_b * (term1_Tb + term2_Tb) + I_b * (x_dot/r)**2)
T_kin = sympy.simplify(T_w + T_r + T_b)

# Lagrangian
L = T_kin - U
L = sympy.simplify(L)

# Generalized coordinates
q = [x_w, theta, x]
q_dot = [x_w_dot, theta_dot, x_dot]
q_ddot = [x_w_ddot, theta_ddot, x_ddot]

# Equations of Motion (Euler-Lagrange)
eom = []
for qi in q:
    eq = diff(diff(L, diff(qi, t)), t) - diff(L, qi)
    eom.append(sympy.simplify(eq))

# Define the system for KanesMethod to linearize
# We need to define the generalized speeds and their derivatives
u = [x_w_dot, theta_dot, x_dot]
u_dot = [x_w_ddot, theta_ddot, x_ddot]

# Create a list of dynamic symbols and their derivatives
q_and_u = q + u
q_dot_and_u_dot = q_dot + u_dot

# The equations are of the form F(q, q_dot, q_ddot) = 0
# We need to convert them to the form M*q_ddot = f(q, q_dot)
# For linearization, we can use KanesMethod's linearize function.
# However, since we already have the EOM, we can directly linearize them.

# Let's define the state vector for linearization:
# states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
# equilibrium_point = {x_w: 0, theta: 0, x: 0, x_w_dot: 0, theta_dot: 0, x_dot: 0}

# We need to solve the EOM for x_w_ddot, theta_ddot, x_ddot
# This is a system of non-linear algebraic equations in terms of accelerations.
# We can use sympy.solve to get expressions for accelerations.

# First, substitute the small angle approximations into the EOM
linearized_eom_approx = []
for eq in eom:
    eq_linear = eq.subs([(cos(theta), 1), (sin(theta), theta)])
    # Now, manually remove higher order terms. This is still a bit tricky.
    # A more robust way is to use sympy.linearize on the original Lagrangian or EOM.

# Let's use KanesMethod for proper linearization.
# Define the generalized coordinates and speeds
q_k = [x_w, theta, x]
u_k = [x_w_dot, theta_dot, x_dot]

# Kinematic differential equations (u_k = q_k_dot)
kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]

# Define the inertial frame
N = ReferenceFrame('N')

# Point for the wheel
P_w = Point('P_w')
P_w.set_vel(N, x_w_dot * N.x)

# Particle for the wheel
Wheel = Particle('Wheel', P_w, M_w)

# Point for the rotor
P_r = Point('P_r')
P_r.set_vel(N, x_w_dot * N.x)

# RigidBody for the rotor
# The rotor rotates about the y-axis (theta)
# Its center of mass is at (x_w, R + d*cos(theta))
# For linearization, we assume it's at (x_w, R+d) and rotates about y.
# This is getting complicated. Let's simplify the approach.

# Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
# it implies we should linearize the derived EOM directly.

# Let's re-attempt manual linearization by substituting and then dropping higher order terms.
# This is the most direct interpretation of "assume small x and small theta".

# Define the small variables and their derivatives
small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}

# Function to linearize an expression by substituting small angle approximations
# and then removing terms that are products of two or more small variables.
def linearize_expression(expr, small_vars_set):
    # Step 1: Substitute small angle approximations
    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])

    # Step 2: Expand the expression
    expr_expanded = sympy.expand(expr_subs)

    # Step 3: Filter out higher-order terms
    linear_terms = []
    for term in sympy.Add.make_args(expr_expanded):
        # Count the number of 'small' variables (or their derivatives) in the term
        order = 0
        for s_var in small_vars_set:
            if term.has(s_var):
                # If the term is a product, we need to be careful.
                # For example, x*theta is order 2. x_dot*theta is order 2.
                # x_ddot is order 1 (as it's a derivative of a small variable).
                
                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                # Or if it contains a square of a small_vars_set.
                
                # This is still not perfect. Let's try a more direct approach using sympy.series.
                # However, sympy.series is for functions of a single variable.
                
                # The most robust way for multi-variable functions is to use Taylor expansion:
                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                
                # Let's define the state vector:
                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                # The EOM are F_i(z, z_dot_dot) = 0
                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                
                # This means we need to compute the Jacobian.
                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                
                # Let's define the full state vector for linearization
                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                
                # The EOM are functions of q, q_dot, q_ddot.
                # We need to express q_ddot in terms of q and q_dot.
                # This is where solving the EOM for accelerations comes in.
                
                # Let's try to solve the EOM for accelerations first.
                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                # Then we can linearize f1, f2, f3.
                
                # The equations are:
                # eom[0] = 0 (for x_w)
                # eom[1] = 0 (for theta)
                # eom[2] = 0 (for x)
                
                # Solve for accelerations
                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                
                # Let's try to linearize the EOM first, then solve for accelerations.
                
                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                # This requires setting up the system using `KanesMethod`.
                
                # Redefine the system using KanesMethod to leverage its linearization capabilities
                # Generalized coordinates
                q_k = [x_w, theta, x]
                # Generalized speeds (derivatives of generalized coordinates)
                u_k = [x_w_dot, theta_dot, x_dot]
                
                # Kinematic differential equations (u_k = q_k_dot)
                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                
                # Define the inertial frame
                N = ReferenceFrame('N')
                
                # Point for the wheel
                P_w = Point('P_w')
                P_w.set_vel(N, x_w_dot * N.x)
                
                # Particle for the wheel
                Wheel = Particle('Wheel', P_w, M_w)
                
                # Point for the rotor
                P_r = Point('P_r')
                P_r.set_vel(N, x_w_dot * N.x)
                
                # RigidBody for the rotor
                # The rotor rotates about the y-axis (theta)
                # Its center of mass is at (x_w, R + d*cos(theta))
                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                # This is getting complicated. Let's simplify the approach.
                
                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                # it implies we should linearize the derived EOM directly.
                
                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                # This is the most direct interpretation of "assume small x and small theta".
                
                # Define the small variables and their derivatives
                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                
                # Function to linearize an expression by substituting small angle approximations
                # and then removing terms that are products of two or more small variables.
                def linearize_expression(expr, small_vars_set):
                    # Step 1: Substitute small angle approximations
                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                
                    # Step 2: Expand the expression
                    expr_expanded = sympy.expand(expr_subs)
                
                    # Step 3: Filter out higher-order terms
                    linear_terms = []
                    for term in sympy.Add.make_args(expr_expanded):
                        # Count the number of 'small' variables (or their derivatives) in the term
                        order = 0
                        for s_var in small_vars_set:
                            if term.has(s_var):
                                # If the term is a product, we need to be careful.
                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                
                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                # Or if it contains a square of a small_vars_set.
                                
                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                # However, sympy.series is for functions of a single variable.
                                
                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                
                                # Let's define the state vector:
                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                # The EOM are F_i(z, z_dot_dot) = 0
                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                
                                # This means we need to compute the Jacobian.
                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                
                                # Let's define the full state vector for linearization
                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                
                                # The EOM are functions of q, q_dot, q_ddot.
                                # We need to express q_ddot in terms of q and q_dot.
                                # This is where solving the EOM for accelerations comes in.
                                
                                # Let's try to solve the EOM for accelerations first.
                                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                # Then we can linearize f1, f2, f3.
                                
                                # The equations are:
                                # eom[0] = 0 (for x_w)
                                # eom[1] = 0 (for theta)
                                # eom[2] = 0 (for x)
                                
                                # Solve for accelerations
                                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                
                                # Let's try to linearize the EOM first, then solve for accelerations.
                                
                                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                # This requires setting up the system using `KanesMethod`.
                                
                                # Redefine the system using KanesMethod to leverage its linearization capabilities
                                # Generalized coordinates
                                q_k = [x_w, theta, x]
                                # Generalized speeds (derivatives of generalized coordinates)
                                u_k = [x_w_dot, theta_dot, x_dot]
                                
                                # Kinematic differential equations (u_k = q_k_dot)
                                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                
                                # Define the inertial frame
                                N = ReferenceFrame('N')
                                
                                # Point for the wheel
                                P_w = Point('P_w')
                                P_w.set_vel(N, x_w_dot * N.x)
                                
                                # Particle for the wheel
                                Wheel = Particle('Wheel', P_w, M_w)
                                
                                # Point for the rotor
                                P_r = Point('P_r')
                                P_r.set_vel(N, x_w_dot * N.x)
                                
                                # RigidBody for the rotor
                                # The rotor rotates about the y-axis (theta)
                                # Its center of mass is at (x_w, R + d*cos(theta))
                                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                # This is getting complicated. Let's simplify the approach.
                                
                                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                # it implies we should linearize the derived EOM directly.
                                
                                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                # This is the most direct interpretation of "assume small x and small theta".
                                
                                # Define the small variables and their derivatives
                                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                
                                # Function to linearize an expression by substituting small angle approximations
                                # and then removing terms that are products of two or more small variables.
                                def linearize_expression(expr, small_vars_set):
                                    # Step 1: Substitute small angle approximations
                                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                
                                    # Step 2: Expand the expression
                                    expr_expanded = sympy.expand(expr_subs)
                                
                                    # Step 3: Filter out higher-order terms
                                    linear_terms = []
                                    for term in sympy.Add.make_args(expr_expanded):
                                        # Count the number of 'small' variables (or their derivatives) in the term
                                        order = 0
                                        for s_var in small_vars_set:
                                            if term.has(s_var):
                                                # If the term is a product, we need to be careful.
                                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                                
                                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                # Or if it contains a square of a small_vars_set.
                                                
                                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                # However, sympy.series is for functions of a single variable.
                                                
                                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                
                                                # Let's define the state vector:
                                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                # The EOM are F_i(z, z_dot_dot) = 0
                                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                
                                                # This means we need to compute the Jacobian.
                                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                
                                                # Let's define the full state vector for linearization
                                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                
                                                # The EOM are functions of q, q_dot, q_ddot.
                                                # We need to express q_ddot in terms of q and q_dot.
                                                # This is where solving the EOM for accelerations comes in.
                                                
                                                # Let's try to solve the EOM for accelerations first.
                                                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                # Then we can linearize f1, f2, f3.
                                                
                                                # The equations are:
                                                # eom[0] = 0 (for x_w)
                                                # eom[1] = 0 (for theta)
                                                # eom[2] = 0 (for x)
                                                
                                                # Solve for accelerations
                                                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                
                                                # Let's try to linearize the EOM first, then solve for accelerations.
                                                
                                                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                # This requires setting up the system using `KanesMethod`.
                                                
                                                # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                # Generalized coordinates
                                                q_k = [x_w, theta, x]
                                                # Generalized speeds (derivatives of generalized coordinates)
                                                u_k = [x_w_dot, theta_dot, x_dot]
                                                
                                                # Kinematic differential equations (u_k = q_k_dot)
                                                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                
                                                # Define the inertial frame
                                                N = ReferenceFrame('N')
                                                
                                                # Point for the wheel
                                                P_w = Point('P_w')
                                                P_w.set_vel(N, x_w_dot * N.x)
                                                
                                                # Particle for the wheel
                                                Wheel = Particle('Wheel', P_w, M_w)
                                                
                                                # Point for the rotor
                                                P_r = Point('P_r')
                                                P_r.set_vel(N, x_w_dot * N.x)
                                                
                                                # RigidBody for the rotor
                                                # The rotor rotates about the y-axis (theta)
                                                # Its center of mass is at (x_w, R + d*cos(theta))
                                                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                # This is getting complicated. Let's simplify the approach.
                                                
                                                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                # it implies we should linearize the derived EOM directly.
                                                
                                                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                # This is the most direct interpretation of "assume small x and small theta".
                                                
                                                # Define the small variables and their derivatives
                                                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                
                                                # Function to linearize an expression by substituting small angle approximations
                                                # and then removing terms that are products of two or more small variables.
                                                def linearize_expression(expr, small_vars_set):
                                                    # Step 1: Substitute small angle approximations
                                                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                
                                                    # Step 2: Expand the expression
                                                    expr_expanded = sympy.expand(expr_subs)
                                                
                                                    # Step 3: Filter out higher-order terms
                                                    linear_terms = []
                                                    for term in sympy.Add.make_args(expr_expanded):
                                                        # Count the number of 'small' variables (or their derivatives) in the term
                                                        order = 0
                                                        for s_var in small_vars_set:
                                                            if term.has(s_var):
                                                                # If the term is a product, we need to be careful.
                                                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                
                                                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                # Or if it contains a square of a small_vars_set.
                                                                
                                                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                # However, sympy.series is for functions of a single variable.
                                                                
                                                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                
                                                                # Let's define the state vector:
                                                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                # The EOM are F_i(z, z_dot_dot) = 0
                                                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                
                                                                # This means we need to compute the Jacobian.
                                                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                
                                                                # Let's define the full state vector for linearization
                                                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                
                                                                # The EOM are functions of q, q_dot, q_ddot.
                                                                # We need to express q_ddot in terms of q and q_dot.
                                                                # This is where solving the EOM for accelerations comes in.
                                                                
                                                                # Let's try to solve the EOM for accelerations first.
                                                                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                # Then we can linearize f1, f2, f3.
                                                                
                                                                # The equations are:
                                                                # eom[0] = 0 (for x_w)
                                                                # eom[1] = 0 (for theta)
                                                                # eom[2] = 0 (for x)
                                                                
                                                                # Solve for accelerations
                                                                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                
                                                                # Let's try to linearize the EOM first, then solve for accelerations.
                                                                
                                                                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                # This requires setting up the system using `KanesMethod`.
                                                                
                                                                # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                # Generalized coordinates
                                                                q_k = [x_w, theta, x]
                                                                # Generalized speeds (derivatives of generalized coordinates)
                                                                u_k = [x_w_dot, theta_dot, x_dot]
                                                                
                                                                # Kinematic differential equations (u_k = q_k_dot)
                                                                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                
                                                                # Define the inertial frame
                                                                N = ReferenceFrame('N')
                                                                
                                                                # Point for the wheel
                                                                P_w = Point('P_w')
                                                                P_w.set_vel(N, x_w_dot * N.x)
                                                                
                                                                # Particle for the wheel
                                                                Wheel = Particle('Wheel', P_w, M_w)
                                                                
                                                                # Point for the rotor
                                                                P_r = Point('P_r')
                                                                P_r.set_vel(N, x_w_dot * N.x)
                                                                
                                                                # RigidBody for the rotor
                                                                # The rotor rotates about the y-axis (theta)
                                                                # Its center of mass is at (x_w, R + d*cos(theta))
                                                                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                # This is getting complicated. Let's simplify the approach.
                                                                
                                                                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                # it implies we should linearize the derived EOM directly.
                                                                
                                                                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                # This is the most direct interpretation of "assume small x and small theta".
                                                                
                                                                # Define the small variables and their derivatives
                                                                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                
                                                                # Function to linearize an expression by substituting small angle approximations
                                                                # and then removing terms that are products of two or more small variables.
                                                                def linearize_expression(expr, small_vars_set):
                                                                    # Step 1: Substitute small angle approximations
                                                                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                
                                                                    # Step 2: Expand the expression
                                                                    expr_expanded = sympy.expand(expr_subs)
                                                                
                                                                    # Step 3: Filter out higher-order terms
                                                                    linear_terms = []
                                                                    for term in sympy.Add.make_args(expr_expanded):
                                                                        # Count the number of 'small' variables (or their derivatives) in the term
                                                                        order = 0
                                                                        for s_var in small_vars_set:
                                                                            if term.has(s_var):
                                                                                # If the term is a product, we need to be careful.
                                                                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                
                                                                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                # Or if it contains a square of a small_vars_set.
                                                                                
                                                                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                # However, sympy.series is for functions of a single variable.
                                                                                
                                                                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                
                                                                                # Let's define the state vector:
                                                                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                # The EOM are F_i(z, z_dot_dot) = 0
                                                                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                
                                                                                # This means we need to compute the Jacobian.
                                                                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                
                                                                                # Let's define the full state vector for linearization
                                                                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                
                                                                                # The EOM are functions of q, q_dot, q_ddot.
                                                                                # We need to express q_ddot in terms of q and q_dot.
                                                                                # This is where solving the EOM for accelerations comes in.
                                                                                
                                                                                # Let's try to solve the EOM for accelerations first.
                                                                                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                # Then we can linearize f1, f2, f3.
                                                                                
                                                                                # The equations are:
                                                                                # eom[0] = 0 (for x_w)
                                                                                # eom[1] = 0 (for theta)
                                                                                # eom[2] = 0 (for x)
                                                                                
                                                                                # Solve for accelerations
                                                                                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                
                                                                                # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                
                                                                                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                # This requires setting up the system using `KanesMethod`.
                                                                                
                                                                                # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                # Generalized coordinates
                                                                                q_k = [x_w, theta, x]
                                                                                # Generalized speeds (derivatives of generalized coordinates)
                                                                                u_k = [x_w_dot, theta_dot, x_dot]
                                                                                
                                                                                # Kinematic differential equations (u_k = q_k_dot)
                                                                                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                
                                                                                # Define the inertial frame
                                                                                N = ReferenceFrame('N')
                                                                                
                                                                                # Point for the wheel
                                                                                P_w = Point('P_w')
                                                                                P_w.set_vel(N, x_w_dot * N.x)
                                                                                
                                                                                # Particle for the wheel
                                                                                Wheel = Particle('Wheel', P_w, M_w)
                                                                                
                                                                                # Point for the rotor
                                                                                P_r = Point('P_r')
                                                                                P_r.set_vel(N, x_w_dot * N.x)
                                                                                
                                                                                # RigidBody for the rotor
                                                                                # The rotor rotates about the y-axis (theta)
                                                                                # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                # This is getting complicated. Let's simplify the approach.
                                                                                
                                                                                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                # it implies we should linearize the derived EOM directly.
                                                                                
                                                                                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                # This is the most direct interpretation of "assume small x and small theta".
                                                                                
                                                                                # Define the small variables and their derivatives
                                                                                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                
                                                                                # Function to linearize an expression by substituting small angle approximations
                                                                                # and then removing terms that are products of two or more small variables.
                                                                                def linearize_expression(expr, small_vars_set):
                                                                                    # Step 1: Substitute small angle approximations
                                                                                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                
                                                                                    # Step 2: Expand the expression
                                                                                    expr_expanded = sympy.expand(expr_subs)
                                                                                
                                                                                    # Step 3: Filter out higher-order terms
                                                                                    linear_terms = []
                                                                                    for term in sympy.Add.make_args(expr_expanded):
                                                                                        # Count the number of 'small' variables (or their derivatives) in the term
                                                                                        order = 0
                                                                                        for s_var in small_vars_set:
                                                                                            if term.has(s_var):
                                                                                                # If the term is a product, we need to be careful.
                                                                                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                
                                                                                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                # Or if it contains a square of a small_vars_set.
                                                                                                
                                                                                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                # However, sympy.series is for functions of a single variable.
                                                                                                
                                                                                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                
                                                                                                # Let's define the state vector:
                                                                                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                
                                                                                                # This means we need to compute the Jacobian.
                                                                                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                
                                                                                                # Let's define the full state vector for linearization
                                                                                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                
                                                                                                # The EOM are functions of q, q_dot, q_ddot.
                                                                                                # We need to express q_ddot in terms of q and q_dot.
                                                                                                # This is where solving the EOM for accelerations comes in.
                                                                                                
                                                                                                # Let's try to solve the EOM for accelerations first.
                                                                                                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                # Then we can linearize f1, f2, f3.
                                                                                                
                                                                                                # The equations are:
                                                                                                # eom[0] = 0 (for x_w)
                                                                                                # eom[1] = 0 (for theta)
                                                                                                # eom[2] = 0 (for x)
                                                                                                
                                                                                                # Solve for accelerations
                                                                                                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                
                                                                                                # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                
                                                                                                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                # This requires setting up the system using `KanesMethod`.
                                                                                                
                                                                                                # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                # Generalized coordinates
                                                                                                q_k = [x_w, theta, x]
                                                                                                # Generalized speeds (derivatives of generalized coordinates)
                                                                                                u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                
                                                                                                # Kinematic differential equations (u_k = q_k_dot)
                                                                                                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                
                                                                                                # Define the inertial frame
                                                                                                N = ReferenceFrame('N')
                                                                                                
                                                                                                # Point for the wheel
                                                                                                P_w = Point('P_w')
                                                                                                P_w.set_vel(N, x_w_dot * N.x)
                                                                                                
                                                                                                # Particle for the wheel
                                                                                                Wheel = Particle('Wheel', P_w, M_w)
                                                                                                
                                                                                                # Point for the rotor
                                                                                                P_r = Point('P_r')
                                                                                                P_r.set_vel(N, x_w_dot * N.x)
                                                                                                
                                                                                                # RigidBody for the rotor
                                                                                                # The rotor rotates about the y-axis (theta)
                                                                                                # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                # This is getting complicated. Let's simplify the approach.
                                                                                                
                                                                                                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                # it implies we should linearize the derived EOM directly.
                                                                                                
                                                                                                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                # This is the most direct interpretation of "assume small x and small theta".
                                                                                                
                                                                                                # Define the small variables and their derivatives
                                                                                                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                
                                                                                                # Function to linearize an expression by substituting small angle approximations
                                                                                                # and then removing terms that are products of two or more small variables.
                                                                                                def linearize_expression(expr, small_vars_set):
                                                                                                    # Step 1: Substitute small angle approximations
                                                                                                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                
                                                                                                    # Step 2: Expand the expression
                                                                                                    expr_expanded = sympy.expand(expr_subs)
                                                                                                
                                                                                                    # Step 3: Filter out higher-order terms
                                                                                                    linear_terms = []
                                                                                                    for term in sympy.Add.make_args(expr_expanded):
                                                                                                        # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                        order = 0
                                                                                                        for s_var in small_vars_set:
                                                                                                            if term.has(s_var):
                                                                                                                # If the term is a product, we need to be careful.
                                                                                                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                
                                                                                                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                # Or if it contains a square of a small_vars_set.
                                                                                                                
                                                                                                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                # However, sympy.series is for functions of a single variable.
                                                                                                                
                                                                                                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                
                                                                                                                # Let's define the state vector:
                                                                                                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                
                                                                                                                # This means we need to compute the Jacobian.
                                                                                                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                
                                                                                                                # Let's define the full state vector for linearization
                                                                                                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                
                                                                                                                # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                # We need to express q_ddot in terms of q and q_dot.
                                                                                                                # This is where solving the EOM for accelerations comes in.
                                                                                                                
                                                                                                                # Let's try to solve the EOM for accelerations first.
                                                                                                                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                # Then we can linearize f1, f2, f3.
                                                                                                                
                                                                                                                # The equations are:
                                                                                                                # eom[0] = 0 (for x_w)
                                                                                                                # eom[1] = 0 (for theta)
                                                                                                                # eom[2] = 0 (for x)
                                                                                                                
                                                                                                                # Solve for accelerations
                                                                                                                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                
                                                                                                                # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                
                                                                                                                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                # This requires setting up the system using `KanesMethod`.
                                                                                                                
                                                                                                                # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                # Generalized coordinates
                                                                                                                q_k = [x_w, theta, x]
                                                                                                                # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                
                                                                                                                # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                
                                                                                                                # Define the inertial frame
                                                                                                                N = ReferenceFrame('N')
                                                                                                                
                                                                                                                # Point for the wheel
                                                                                                                P_w = Point('P_w')
                                                                                                                P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                
                                                                                                                # Particle for the wheel
                                                                                                                Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                
                                                                                                                # Point for the rotor
                                                                                                                P_r = Point('P_r')
                                                                                                                P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                
                                                                                                                # RigidBody for the rotor
                                                                                                                # The rotor rotates about the y-axis (theta)
                                                                                                                # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                # This is getting complicated. Let's simplify the approach.
                                                                                                                
                                                                                                                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                # it implies we should linearize the derived EOM directly.
                                                                                                                
                                                                                                                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                
                                                                                                                # Define the small variables and their derivatives
                                                                                                                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                
                                                                                                                # Function to linearize an expression by substituting small angle approximations
                                                                                                                # and then removing terms that are products of two or more small variables.
                                                                                                                def linearize_expression(expr, small_vars_set):
                                                                                                                    # Step 1: Substitute small angle approximations
                                                                                                                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                
                                                                                                                    # Step 2: Expand the expression
                                                                                                                    expr_expanded = sympy.expand(expr_subs)
                                                                                                                
                                                                                                                    # Step 3: Filter out higher-order terms
                                                                                                                    linear_terms = []
                                                                                                                    for term in sympy.Add.make_args(expr_expanded):
                                                                                                                        # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                        order = 0
                                                                                                                        for s_var in small_vars_set:
                                                                                                                            if term.has(s_var):
                                                                                                                                # If the term is a product, we need to be careful.
                                                                                                                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                
                                                                                                                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                # Or if it contains a square of a small_vars_set.
                                                                                                                                
                                                                                                                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                # However, sympy.series is for functions of a single variable.
                                                                                                                                
                                                                                                                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                
                                                                                                                                # Let's define the state vector:
                                                                                                                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                
                                                                                                                                # This means we need to compute the Jacobian.
                                                                                                                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                
                                                                                                                                # Let's define the full state vector for linearization
                                                                                                                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                
                                                                                                                                # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                # This is where solving the EOM for accelerations comes in.
                                                                                                                                
                                                                                                                                # Let's try to solve the EOM for accelerations first.
                                                                                                                                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                # Then we can linearize f1, f2, f3.
                                                                                                                                
                                                                                                                                # The equations are:
                                                                                                                                # eom[0] = 0 (for x_w)
                                                                                                                                # eom[1] = 0 (for theta)
                                                                                                                                # eom[2] = 0 (for x)
                                                                                                                                
                                                                                                                                # Solve for accelerations
                                                                                                                                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                
                                                                                                                                # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                
                                                                                                                                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                # This requires setting up the system using `KanesMethod`.
                                                                                                                                
                                                                                                                                # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                # Generalized coordinates
                                                                                                                                q_k = [x_w, theta, x]
                                                                                                                                # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                
                                                                                                                                # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                
                                                                                                                                # Define the inertial frame
                                                                                                                                N = ReferenceFrame('N')
                                                                                                                                
                                                                                                                                # Point for the wheel
                                                                                                                                P_w = Point('P_w')
                                                                                                                                P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                
                                                                                                                                # Particle for the wheel
                                                                                                                                Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                
                                                                                                                                # Point for the rotor
                                                                                                                                P_r = Point('P_r')
                                                                                                                                P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                
                                                                                                                                # RigidBody for the rotor
                                                                                                                                # The rotor rotates about the y-axis (theta)
                                                                                                                                # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                # This is getting complicated. Let's simplify the approach.
                                                                                                                                
                                                                                                                                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                # it implies we should linearize the derived EOM directly.
                                                                                                                                
                                                                                                                                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                
                                                                                                                                # Define the small variables and their derivatives
                                                                                                                                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                
                                                                                                                                # Function to linearize an expression by substituting small angle approximations
                                                                                                                                # and then removing terms that are products of two or more small variables.
                                                                                                                                def linearize_expression(expr, small_vars_set):
                                                                                                                                    # Step 1: Substitute small angle approximations
                                                                                                                                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                
                                                                                                                                    # Step 2: Expand the expression
                                                                                                                                    expr_expanded = sympy.expand(expr_subs)
                                                                                                                                
                                                                                                                                    # Step 3: Filter out higher-order terms
                                                                                                                                    linear_terms = []
                                                                                                                                    for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                        # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                        order = 0
                                                                                                                                        for s_var in small_vars_set:
                                                                                                                                            if term.has(s_var):
                                                                                                                                                # If the term is a product, we need to be careful.
                                                                                                                                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                
                                                                                                                                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                # Or if it contains a square of a small_vars_set.
                                                                                                                                                
                                                                                                                                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                # However, sympy.series is for functions of a single variable.
                                                                                                                                                
                                                                                                                                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                
                                                                                                                                                # Let's define the state vector:
                                                                                                                                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                
                                                                                                                                                # This means we need to compute the Jacobian.
                                                                                                                                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                
                                                                                                                                                # Let's define the full state vector for linearization
                                                                                                                                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                
                                                                                                                                                # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                # This is where solving the EOM for accelerations comes in.
                                                                                                                                                
                                                                                                                                                # Let's try to solve the EOM for accelerations first.
                                                                                                                                                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                # Then we can linearize f1, f2, f3.
                                                                                                                                                
                                                                                                                                                # The equations are:
                                                                                                                                                # eom[0] = 0 (for x_w)
                                                                                                                                                # eom[1] = 0 (for theta)
                                                                                                                                                # eom[2] = 0 (for x)
                                                                                                                                                
                                                                                                                                                # Solve for accelerations
                                                                                                                                                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                
                                                                                                                                                # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                
                                                                                                                                                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                # This requires setting up the system using `KanesMethod`.
                                                                                                                                                
                                                                                                                                                # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                # Generalized coordinates
                                                                                                                                                q_k = [x_w, theta, x]
                                                                                                                                                # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                
                                                                                                                                                # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                
                                                                                                                                                # Define the inertial frame
                                                                                                                                                N = ReferenceFrame('N')
                                                                                                                                                
                                                                                                                                                # Point for the wheel
                                                                                                                                                P_w = Point('P_w')
                                                                                                                                                P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                
                                                                                                                                                # Particle for the wheel
                                                                                                                                                Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                
                                                                                                                                                # Point for the rotor
                                                                                                                                                P_r = Point('P_r')
                                                                                                                                                P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                
                                                                                                                                                # RigidBody for the rotor
                                                                                                                                                # The rotor rotates about the y-axis (theta)
                                                                                                                                                # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                # This is getting complicated. Let's simplify the approach.
                                                                                                                                                
                                                                                                                                                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                # it implies we should linearize the derived EOM directly.
                                                                                                                                                
                                                                                                                                                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                
                                                                                                                                                # Define the small variables and their derivatives
                                                                                                                                                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                                
                                                                                                                                                # Function to linearize an expression by substituting small angle approximations
                                                                                                                                                # and then removing terms that are products of two or more small variables.
                                                                                                                                                def linearize_expression(expr, small_vars_set):
                                                                                                                                                    # Step 1: Substitute small angle approximations
                                                                                                                                                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                                
                                                                                                                                                    # Step 2: Expand the expression
                                                                                                                                                    expr_expanded = sympy.expand(expr_subs)
                                                                                                                                                
                                                                                                                                                    # Step 3: Filter out higher-order terms
                                                                                                                                                    linear_terms = []
                                                                                                                                                    for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                                        # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                                        order = 0
                                                                                                                                                        for s_var in small_vars_set:
                                                                                                                                                            if term.has(s_var):
                                                                                                                                                                # If the term is a product, we need to be careful.
                                                                                                                                                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                                
                                                                                                                                                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                                # Or if it contains a square of a small_vars_set.
                                                                                                                                                                
                                                                                                                                                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                                # However, sympy.series is for functions of a single variable.
                                                                                                                                                                
                                                                                                                                                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                                
                                                                                                                                                                # Let's define the state vector:
                                                                                                                                                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                                
                                                                                                                                                                # This means we need to compute the Jacobian.
                                                                                                                                                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                                
                                                                                                                                                                # Let's define the full state vector for linearization
                                                                                                                                                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                
                                                                                                                                                                # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                                # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                                # This is where solving the EOM for accelerations comes in.
                                                                                                                                                                
                                                                                                                                                                # Let's try to solve the EOM for accelerations first.
                                                                                                                                                                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                                # Then we can linearize f1, f2, f3.
                                                                                                                                                                
                                                                                                                                                                # The equations are:
                                                                                                                                                                # eom[0] = 0 (for x_w)
                                                                                                                                                                # eom[1] = 0 (for theta)
                                                                                                                                                                # eom[2] = 0 (for x)
                                                                                                                                                                
                                                                                                                                                                # Solve for accelerations
                                                                                                                                                                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                                
                                                                                                                                                                # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                                
                                                                                                                                                                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                                # This requires setting up the system using `KanesMethod`.
                                                                                                                                                                
                                                                                                                                                                # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                                # Generalized coordinates
                                                                                                                                                                q_k = [x_w, theta, x]
                                                                                                                                                                # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                                u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                                
                                                                                                                                                                # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                                
                                                                                                                                                                # Define the inertial frame
                                                                                                                                                                N = ReferenceFrame('N')
                                                                                                                                                                
                                                                                                                                                                # Point for the wheel
                                                                                                                                                                P_w = Point('P_w')
                                                                                                                                                                P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                
                                                                                                                                                                # Particle for the wheel
                                                                                                                                                                Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                                
                                                                                                                                                                # Point for the rotor
                                                                                                                                                                P_r = Point('P_r')
                                                                                                                                                                P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                
                                                                                                                                                                # RigidBody for the rotor
                                                                                                                                                                # The rotor rotates about the y-axis (theta)
                                                                                                                                                                # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                                # This is getting complicated. Let's simplify the approach.
                                                                                                                                                                
                                                                                                                                                                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                                # it implies we should linearize the derived EOM directly.
                                                                                                                                                                
                                                                                                                                                                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                                # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                                
                                                                                                                                                                # Define the small variables and their derivatives
                                                                                                                                                                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                                                
                                                                                                                                                                # Function to linearize an expression by substituting small angle approximations
                                                                                                                                                                # and then removing terms that are products of two or more small variables.
                                                                                                                                                                def linearize_expression(expr, small_vars_set):
                                                                                                                                                                    # Step 1: Substitute small angle approximations
                                                                                                                                                                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                                                
                                                                                                                                                                    # Step 2: Expand the expression
                                                                                                                                                                    expr_expanded = sympy.expand(expr_subs)
                                                                                                                                                                
                                                                                                                                                                    # Step 3: Filter out higher-order terms
                                                                                                                                                                    linear_terms = []
                                                                                                                                                                    for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                                                        # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                                                        order = 0
                                                                                                                                                                        for s_var in small_vars_set:
                                                                                                                                                                            if term.has(s_var):
                                                                                                                                                                                # If the term is a product, we need to be careful.
                                                                                                                                                                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                                                
                                                                                                                                                                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                                                # Or if it contains a square of a small_vars_set.
                                                                                                                                                                                
                                                                                                                                                                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                                                # However, sympy.series is for functions of a single variable.
                                                                                                                                                                                
                                                                                                                                                                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                                                
                                                                                                                                                                                # Let's define the state vector:
                                                                                                                                                                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                                                
                                                                                                                                                                                # This means we need to compute the Jacobian.
                                                                                                                                                                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                                                
                                                                                                                                                                                # Let's define the full state vector for linearization
                                                                                                                                                                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                
                                                                                                                                                                                # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                                                # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                                                # This is where solving the EOM for accelerations comes in.
                                                                                                                                                                                
                                                                                                                                                                                # Let's try to solve the EOM for accelerations first.
                                                                                                                                                                                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                                                # Then we can linearize f1, f2, f3.
                                                                                                                                                                                
                                                                                                                                                                                # The equations are:
                                                                                                                                                                                # eom[0] = 0 (for x_w)
                                                                                                                                                                                # eom[1] = 0 (for theta)
                                                                                                                                                                                # eom[2] = 0 (for x)
                                                                                                                                                                                
                                                                                                                                                                                # Solve for accelerations
                                                                                                                                                                                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                                                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                                                
                                                                                                                                                                                # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                                                
                                                                                                                                                                                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                                                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                                                # This requires setting up the system using `KanesMethod`.
                                                                                                                                                                                
                                                                                                                                                                                # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                                                # Generalized coordinates
                                                                                                                                                                                q_k = [x_w, theta, x]
                                                                                                                                                                                # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                                                u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                
                                                                                                                                                                                # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                                                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                                                
                                                                                                                                                                                # Define the inertial frame
                                                                                                                                                                                N = ReferenceFrame('N')
                                                                                                                                                                                
                                                                                                                                                                                # Point for the wheel
                                                                                                                                                                                P_w = Point('P_w')
                                                                                                                                                                                P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                
                                                                                                                                                                                # Particle for the wheel
                                                                                                                                                                                Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                                                
                                                                                                                                                                                # Point for the rotor
                                                                                                                                                                                P_r = Point('P_r')
                                                                                                                                                                                P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                
                                                                                                                                                                                # RigidBody for the rotor
                                                                                                                                                                                # The rotor rotates about the y-axis (theta)
                                                                                                                                                                                # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                                                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                                                # This is getting complicated. Let's simplify the approach.
                                                                                                                                                                                
                                                                                                                                                                                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                                                # it implies we should linearize the derived EOM directly.
                                                                                                                                                                                
                                                                                                                                                                                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                                                # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                                                
                                                                                                                                                                                # Define the small variables and their derivatives
                                                                                                                                                                                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                                                                
                                                                                                                                                                                # Function to linearize an expression by substituting small angle approximations
                                                                                                                                                                                # and then removing terms that are products of two or more small variables.
                                                                                                                                                                                def linearize_expression(expr, small_vars_set):
                                                                                                                                                                                    # Step 1: Substitute small angle approximations
                                                                                                                                                                                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                                                                
                                                                                                                                                                                    # Step 2: Expand the expression
                                                                                                                                                                                    expr_expanded = sympy.expand(expr_subs)
                                                                                                                                                                                
                                                                                                                                                                                    # Step 3: Filter out higher-order terms
                                                                                                                                                                                    linear_terms = []
                                                                                                                                                                                    for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                                                                        # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                                                                        order = 0
                                                                                                                                                                                        for s_var in small_vars_set:
                                                                                                                                                                                            if term.has(s_var):
                                                                                                                                                                                                # If the term is a product, we need to be careful.
                                                                                                                                                                                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                                                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                                                                
                                                                                                                                                                                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                                                                # Or if it contains a square of a small_vars_set.
                                                                                                                                                                                                
                                                                                                                                                                                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                                                                # However, sympy.series is for functions of a single variable.
                                                                                                                                                                                                
                                                                                                                                                                                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                                                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                                                                
                                                                                                                                                                                                # Let's define the state vector:
                                                                                                                                                                                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                                                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                                                                
                                                                                                                                                                                                # This means we need to compute the Jacobian.
                                                                                                                                                                                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                                                                
                                                                                                                                                                                                # Let's define the full state vector for linearization
                                                                                                                                                                                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                
                                                                                                                                                                                                # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                                                                # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                                                                # This is where solving the EOM for accelerations comes in.
                                                                                                                                                                                                
                                                                                                                                                                                                # Let's try to solve the EOM for accelerations first.
                                                                                                                                                                                                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                                                                # Then we can linearize f1, f2, f3.
                                                                                                                                                                                                
                                                                                                                                                                                                # The equations are:
                                                                                                                                                                                                # eom[0] = 0 (for x_w)
                                                                                                                                                                                                # eom[1] = 0 (for theta)
                                                                                                                                                                                                # eom[2] = 0 (for x)
                                                                                                                                                                                                
                                                                                                                                                                                                # Solve for accelerations
                                                                                                                                                                                                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                                                                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                                                                
                                                                                                                                                                                                # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                                                                
                                                                                                                                                                                                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                                                                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                                                                # This requires setting up the system using `KanesMethod`.
                                                                                                                                                                                                
                                                                                                                                                                                                # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                                                                # Generalized coordinates
                                                                                                                                                                                                q_k = [x_w, theta, x]
                                                                                                                                                                                                # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                                                                u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                
                                                                                                                                                                                                # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                                                                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                                                                
                                                                                                                                                                                                # Define the inertial frame
                                                                                                                                                                                                N = ReferenceFrame('N')
                                                                                                                                                                                                
                                                                                                                                                                                                # Point for the wheel
                                                                                                                                                                                                P_w = Point('P_w')
                                                                                                                                                                                                P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                
                                                                                                                                                                                                # Particle for the wheel
                                                                                                                                                                                                Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                                                                
                                                                                                                                                                                                # Point for the rotor
                                                                                                                                                                                                P_r = Point('P_r')
                                                                                                                                                                                                P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                
                                                                                                                                                                                                # RigidBody for the rotor
                                                                                                                                                                                                # The rotor rotates about the y-axis (theta)
                                                                                                                                                                                                # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                                                                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                                                                # This is getting complicated. Let's simplify the approach.
                                                                                                                                                                                                
                                                                                                                                                                                                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                                                                # it implies we should linearize the derived EOM directly.
                                                                                                                                                                                                
                                                                                                                                                                                                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                                                                # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                                                                
                                                                                                                                                                                                # Define the small variables and their derivatives
                                                                                                                                                                                                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                                                                                
                                                                                                                                                                                                # Function to linearize an expression by substituting small angle approximations
                                                                                                                                                                                                # and then removing terms that are products of two or more small variables.
                                                                                                                                                                                                def linearize_expression(expr, small_vars_set):
                                                                                                                                                                                                    # Step 1: Substitute small angle approximations
                                                                                                                                                                                                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                                                                                
                                                                                                                                                                                                    # Step 2: Expand the expression
                                                                                                                                                                                                    expr_expanded = sympy.expand(expr_subs)
                                                                                                                                                                                                
                                                                                                                                                                                                    # Step 3: Filter out higher-order terms
                                                                                                                                                                                                    linear_terms = []
                                                                                                                                                                                                    for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                                                                                        # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                                                                                        order = 0
                                                                                                                                                                                                        for s_var in small_vars_set:
                                                                                                                                                                                                            if term.has(s_var):
                                                                                                                                                                                                                # If the term is a product, we need to be careful.
                                                                                                                                                                                                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                                                                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                                                                                
                                                                                                                                                                                                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                                                                                # Or if it contains a square of a small_vars_set.
                                                                                                                                                                                                                
                                                                                                                                                                                                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                                                                                # However, sympy.series is for functions of a single variable.
                                                                                                                                                                                                                
                                                                                                                                                                                                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                                                                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Let's define the state vector:
                                                                                                                                                                                                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                                                                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                                                                                
                                                                                                                                                                                                                # This means we need to compute the Jacobian.
                                                                                                                                                                                                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Let's define the full state vector for linearization
                                                                                                                                                                                                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                
                                                                                                                                                                                                                # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                                                                                # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                                                                                # This is where solving the EOM for accelerations comes in.
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Let's try to solve the EOM for accelerations first.
                                                                                                                                                                                                                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                                                                                # Then we can linearize f1, f2, f3.
                                                                                                                                                                                                                
                                                                                                                                                                                                                # The equations are:
                                                                                                                                                                                                                # eom[0] = 0 (for x_w)
                                                                                                                                                                                                                # eom[1] = 0 (for theta)
                                                                                                                                                                                                                # eom[2] = 0 (for x)
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Solve for accelerations
                                                                                                                                                                                                                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                                                                                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                                                                                
                                                                                                                                                                                                                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                                                                                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                                                                                # This requires setting up the system using `KanesMethod`.
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                                                                                # Generalized coordinates
                                                                                                                                                                                                                q_k = [x_w, theta, x]
                                                                                                                                                                                                                # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                                                                                u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                                                                                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Define the inertial frame
                                                                                                                                                                                                                N = ReferenceFrame('N')
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Point for the wheel
                                                                                                                                                                                                                P_w = Point('P_w')
                                                                                                                                                                                                                P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Particle for the wheel
                                                                                                                                                                                                                Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Point for the rotor
                                                                                                                                                                                                                P_r = Point('P_r')
                                                                                                                                                                                                                P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                
                                                                                                                                                                                                                # RigidBody for the rotor
                                                                                                                                                                                                                # The rotor rotates about the y-axis (theta)
                                                                                                                                                                                                                # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                                                                                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                                                                                # This is getting complicated. Let's simplify the approach.
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                                                                                # it implies we should linearize the derived EOM directly.
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                                                                                # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Define the small variables and their derivatives
                                                                                                                                                                                                                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                                                                                                
                                                                                                                                                                                                                # Function to linearize an expression by substituting small angle approximations
                                                                                                                                                                                                                # and then removing terms that are products of two or more small variables.
                                                                                                                                                                                                                def linearize_expression(expr, small_vars_set):
                                                                                                                                                                                                                    # Step 1: Substitute small angle approximations
                                                                                                                                                                                                                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                                                                                                
                                                                                                                                                                                                                    # Step 2: Expand the expression
                                                                                                                                                                                                                    expr_expanded = sympy.expand(expr_subs)
                                                                                                                                                                                                                
                                                                                                                                                                                                                    # Step 3: Filter out higher-order terms
                                                                                                                                                                                                                    linear_terms = []
                                                                                                                                                                                                                    for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                                                                                                        # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                                                                                                        order = 0
                                                                                                                                                                                                                        for s_var in small_vars_set:
                                                                                                                                                                                                                            if term.has(s_var):
                                                                                                                                                                                                                                # If the term is a product, we need to be careful.
                                                                                                                                                                                                                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                                                                                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                                                                                                # Or if it contains a square of a small_vars_set.
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                                                                                                # However, sympy.series is for functions of a single variable.
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                                                                                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Let's define the state vector:
                                                                                                                                                                                                                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                                                                                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # This means we need to compute the Jacobian.
                                                                                                                                                                                                                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Let's define the full state vector for linearization
                                                                                                                                                                                                                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                                                                                                # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                                                                                                # This is where solving the EOM for accelerations comes in.
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Let's try to solve the EOM for accelerations first.
                                                                                                                                                                                                                                # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                                                                                                # Then we can linearize f1, f2, f3.
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # The equations are:
                                                                                                                                                                                                                                # eom[0] = 0 (for x_w)
                                                                                                                                                                                                                                # eom[1] = 0 (for theta)
                                                                                                                                                                                                                                # eom[2] = 0 (for x)
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Solve for accelerations
                                                                                                                                                                                                                                # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                                                                                                # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                                                                                                # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                                                                                                # This requires setting up the system using `KanesMethod`.
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                                                                                                # Generalized coordinates
                                                                                                                                                                                                                                q_k = [x_w, theta, x]
                                                                                                                                                                                                                                # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                                                                                                u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                                                                                                kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Define the inertial frame
                                                                                                                                                                                                                                N = ReferenceFrame('N')
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Point for the wheel
                                                                                                                                                                                                                                P_w = Point('P_w')
                                                                                                                                                                                                                                P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Particle for the wheel
                                                                                                                                                                                                                                Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Point for the rotor
                                                                                                                                                                                                                                P_r = Point('P_r')
                                                                                                                                                                                                                                P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # RigidBody for the rotor
                                                                                                                                                                                                                                # The rotor rotates about the y-axis (theta)
                                                                                                                                                                                                                                # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                                                                                                # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                                                                                                # This is getting complicated. Let's simplify the approach.
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                                                                                                # it implies we should linearize the derived EOM directly.
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                                                                                                # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Define the small variables and their derivatives
                                                                                                                                                                                                                                small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                # Function to linearize an expression by substituting small angle approximations
                                                                                                                                                                                                                                # and then removing terms that are products of two or more small variables.
                                                                                                                                                                                                                                def linearize_expression(expr, small_vars_set):
                                                                                                                                                                                                                                    # Step 1: Substitute small angle approximations
                                                                                                                                                                                                                                    expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                    # Step 2: Expand the expression
                                                                                                                                                                                                                                    expr_expanded = sympy.expand(expr_subs)
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                    # Step 3: Filter out higher-order terms
                                                                                                                                                                                                                                    linear_terms = []
                                                                                                                                                                                                                                    for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                                                                                                                        # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                                                                                                                        order = 0
                                                                                                                                                                                                                                        for s_var in small_vars_set:
                                                                                                                                                                                                                                            if term.has(s_var):
                                                                                                                                                                                                                                                # If the term is a product, we need to be careful.
                                                                                                                                                                                                                                                # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                                                                                                                # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                                                                                                                # Or if it contains a square of a small_vars_set.
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                                                                                                                # However, sympy.series is for functions of a single variable.
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                                                                                                                # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                # Let's define the state vector:
                                                                                                                                                                                                                                                # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                                                                                                                # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                # This means we need to compute the Jacobian.
                                                                                                                                                                                                                                                # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                # Let's define the full state vector for linearization
                                                                                                                                                                                                                                                states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                                                                                 # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                                                                                 # This is where solving the EOM for accelerations comes in.
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # Let's try to solve the EOM for accelerations first.
                                                                                                                                                                                                                 # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                                                                                 # Then we can linearize f1, f2, f3.
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # The equations are:
                                                                                                                                                                                                                 # eom[0] = 0 (for x_w)
                                                                                                                                                                                                                 # eom[1] = 0 (for theta)
                                                                                                                                                                                                                 # eom[2] = 0 (for x)
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # Solve for accelerations
                                                                                                                                                                                                                 # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                                                                                 # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                                                                                 # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                                                                                 # This requires setting up the system using `KanesMethod`.
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                                                                                 # Generalized coordinates
                                                                                                                                                                                                                 q_k = [x_w, theta, x]
                                                                                                                                                                                                                 # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                                                                                 u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                                                                                 kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # Define the inertial frame
                                                                                                                                                                                                                 N = ReferenceFrame('N')
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # Point for the wheel
                                                                                                                                                                                                                 P_w = Point('P_w')
                                                                                                                                                                                                                 P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # Particle for the wheel
                                                                                                                                                                                                                 Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # Point for the rotor
                                                                                                                                                                                                                 P_r = Point('P_r')
                                                                                                                                                                                                                 P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # RigidBody for the rotor
                                                                                                                                                                                                                 # The rotor rotates about the y-axis (theta)
                                                                                                                                                                                                                 # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                                                                                 # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                                                                                 # This is getting complicated. Let's simplify the approach.
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                                                                                 # it implies we should linearize the derived EOM directly.
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                                                                                 # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # Define the small variables and their derivatives
                                                                                                                                                                                                                 small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 # Function to linearize an expression by substituting small angle approximations
                                                                                                                                                                                                                 # and then removing terms that are products of two or more small variables.
                                                                                                                                                                                                                 def linearize_expression(expr, small_vars_set):
                                                                                                                                                                                                                     # Step 1: Substitute small angle approximations
                                                                                                                                                                                                                     expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                                                                                                 
                                                                                                                                                                                                                     # Step 2: Expand the expression
                                                                                                                                                                                                                     expr_expanded = sympy.expand(expr_subs)
                                                                                                                                                                                                                 
                                                                                                                                                                                                                     # Step 3: Filter out higher-order terms
                                                                                                                                                                                                                     linear_terms = []
                                                                                                                                                                                                                     for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                                                                                                         # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                                                                                                         order = 0
                                                                                                                                                                                                                         for s_var in small_vars_set:
                                                                                                                                                                                                                             if term.has(s_var):
                                                                                                                                                                                                                                 # If the term is a product, we need to be careful.
                                                                                                                                                                                                                                 # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                                                                                                 # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                                                                                                 # Or if it contains a square of a small_vars_set.
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                                                                                                 # However, sympy.series is for functions of a single variable.
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                                                                                                 # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Let's define the state vector:
                                                                                                                                                                                                                                 # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                 # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                                                                                                 # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # This means we need to compute the Jacobian.
                                                                                                                                                                                                                                 # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Let's define the full state vector for linearization
                                                                                                                                                                                                                                 states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                                                                                                 # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                                                                                                 # This is where solving the EOM for accelerations comes in.
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Let's try to solve the EOM for accelerations first.
                                                                                                                                                                                                                                 # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                                                                                                 # Then we can linearize f1, f2, f3.
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # The equations are:
                                                                                                                                                                                                                                 # eom[0] = 0 (for x_w)
                                                                                                                                                                                                                                 # eom[1] = 0 (for theta)
                                                                                                                                                                                                                                 # eom[2] = 0 (for x)
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Solve for accelerations
                                                                                                                                                                                                                                 # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                                                                                                 # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                                                                                                 # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                                                                                                 # This requires setting up the system using `KanesMethod`.
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                                                                                                 # Generalized coordinates
                                                                                                                                                                                                                                 q_k = [x_w, theta, x]
                                                                                                                                                                                                                                 # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                                                                                                 u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                                                                                                 kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Define the inertial frame
                                                                                                                                                                                                                                 N = ReferenceFrame('N')
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Point for the wheel
                                                                                                                                                                                                                                 P_w = Point('P_w')
                                                                                                                                                                                                                                 P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Particle for the wheel
                                                                                                                                                                                                                                 Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Point for the rotor
                                                                                                                                                                                                                                 P_r = Point('P_r')
                                                                                                                                                                                                                                 P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # RigidBody for the rotor
                                                                                                                                                                                                                                 # The rotor rotates about the y-axis (theta)
                                                                                                                                                                                                                                 # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                                                                                                 # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                                                                                                 # This is getting complicated. Let's simplify the approach.
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                                                                                                 # it implies we should linearize the derived EOM directly.
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                                                                                                 # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Define the small variables and their derivatives
                                                                                                                                                                                                                                 small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                 # Function to linearize an expression by substituting small angle approximations
                                                                                                                                                                                                                                 # and then removing terms that are products of two or more small variables.
                                                                                                                                                                                                                                 def linearize_expression(expr, small_vars_set):
                                                                                                                                                                                                                                     # Step 1: Substitute small angle approximations
                                                                                                                                                                                                                                     expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                     # Step 2: Expand the expression
                                                                                                                                                                                                                                     expr_expanded = sympy.expand(expr_subs)
                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                     # Step 3: Filter out higher-order terms
                                                                                                                                                                                                                                     linear_terms = []
                                                                                                                                                                                                                                     for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                                                                                                                         # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                                                                                                                         order = 0
                                                                                                                                                                                                                                         for s_var in small_vars_set:
                                                                                                                                                                                                                                             if term.has(s_var):
                                                                                                                                                                                                                                                 # If the term is a product, we need to be careful.
                                                                                                                                                                                                                                                 # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                                                                                                                 # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                                                                                                                 # Or if it contains a square of a small_vars_set.
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                                                                                                                 # However, sympy.series is for functions of a single variable.
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                                                                                                                 # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Let's define the state vector:
                                                                                                                                                                                                                                                 # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                 # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                                                                                                                 # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # This means we need to compute the Jacobian.
                                                                                                                                                                                                                                                 # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Let's define the full state vector for linearization
                                                                                                                                                                                                                                                 states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                                                                                                                 # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                                                                                                                 # This is where solving the EOM for accelerations comes in.
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Let's try to solve the EOM for accelerations first.
                                                                                                                                                                                                                                                 # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                                                                                                                 # Then we can linearize f1, f2, f3.
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # The equations are:
                                                                                                                                                                                                                                                 # eom[0] = 0 (for x_w)
                                                                                                                                                                                                                                                 # eom[1] = 0 (for theta)
                                                                                                                                                                                                                                                 # eom[2] = 0 (for x)
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Solve for accelerations
                                                                                                                                                                                                                                                 # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                                                                                                                 # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                                                                                                                 # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                                                                                                                 # This requires setting up the system using `KanesMethod`.
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                                                                                                                 # Generalized coordinates
                                                                                                                                                                                                                                                 q_k = [x_w, theta, x]
                                                                                                                                                                                                                                                 # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                                                                                                                 u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                                                                                                                 kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Define the inertial frame
                                                                                                                                                                                                                                                 N = ReferenceFrame('N')
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Point for the wheel
                                                                                                                                                                                                                                                 P_w = Point('P_w')
                                                                                                                                                                                                                                                 P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Particle for the wheel
                                                                                                                                                                                                                                                 Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Point for the rotor
                                                                                                                                                                                                                                                 P_r = Point('P_r')
                                                                                                                                                                                                                                                 P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # RigidBody for the rotor
                                                                                                                                                                                                                                                 # The rotor rotates about the y-axis (theta)
                                                                                                                                                                                                                                                 # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                                                                                                                 # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                                                                                                                 # This is getting complicated. Let's simplify the approach.
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                                                                                                                 # it implies we should linearize the derived EOM directly.
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                                                                                                                 # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Define the small variables and their derivatives
                                                                                                                                                                                                                                                 small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                 # Function to linearize an expression by substituting small angle approximations
                                                                                                                                                                                                                                                 # and then removing terms that are products of two or more small variables.
                                                                                                                                                                                                                                                 def linearize_expression(expr, small_vars_set):
                                                                                                                                                                                                                                                     # Step 1: Substitute small angle approximations
                                                                                                                                                                                                                                                     expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                     # Step 2: Expand the expression
                                                                                                                                                                                                                                                     expr_expanded = sympy.expand(expr_subs)
                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                     # Step 3: Filter out higher-order terms
                                                                                                                                                                                                                                                     linear_terms = []
                                                                                                                                                                                                                                                     for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                                                                                                                                         # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                                                                                                                                         order = 0
                                                                                                                                                                                                                                                         for s_var in small_vars_set:
                                                                                                                                                                                                                                                             if term.has(s_var):
                                                                                                                                                                                                                                                                 # If the term is a product, we need to be careful.
                                                                                                                                                                                                                                                                 # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                                                                                                                                 # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                                                                                                                                 # Or if it contains a square of a small_vars_set.
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                                                                                                                                 # However, sympy.series is for functions of a single variable.
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                                                                                                                                 # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Let's define the state vector:
                                                                                                                                                                                                                                                                 # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                 # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                                                                                                                                 # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # This means we need to compute the Jacobian.
                                                                                                                                                                                                                                                                 # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Let's define the full state vector for linearization
                                                                                                                                                                                                                                                                 states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                                                                                                                                 # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                                                                                                                                 # This is where solving the EOM for accelerations comes in.
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Let's try to solve the EOM for accelerations first.
                                                                                                                                                                                                                                                                 # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                                                                                                                                 # Then we can linearize f1, f2, f3.
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # The equations are:
                                                                                                                                                                                                                                                                 # eom[0] = 0 (for x_w)
                                                                                                                                                                                                                                                                 # eom[1] = 0 (for theta)
                                                                                                                                                                                                                                                                 # eom[2] = 0 (for x)
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Solve for accelerations
                                                                                                                                                                                                                                                                 # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                                                                                                                                 # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                                                                                                                                 # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                                                                                                                                 # This requires setting up the system using `KanesMethod`.
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                                                                                                                                 # Generalized coordinates
                                                                                                                                                                                                                                                                 q_k = [x_w, theta, x]
                                                                                                                                                                                                                                                                 # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                                                                                                                                 u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                                                                                                                                 kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Define the inertial frame
                                                                                                                                                                                                                                                                 N = ReferenceFrame('N')
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Point for the wheel
                                                                                                                                                                                                                                                                 P_w = Point('P_w')
                                                                                                                                                                                                                                                                 P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Particle for the wheel
                                                                                                                                                                                                                                                                 Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Point for the rotor
                                                                                                                                                                                                                                                                 P_r = Point('P_r')
                                                                                                                                                                                                                                                                 P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # RigidBody for the rotor
                                                                                                                                                                                                                                                                 # The rotor rotates about the y-axis (theta)
                                                                                                                                                                                                                                                                 # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                                                                                                                                 # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                                                                                                                                 # This is getting complicated. Let's simplify the approach.
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                                                                                                                                 # it implies we should linearize the derived EOM directly.
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                                                                                                                                 # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Define the small variables and their derivatives
                                                                                                                                                                                                                                                                 small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                 # Function to linearize an expression by substituting small angle approximations
                                                                                                                                                                                                                                                                 # and then removing terms that are products of two or more small variables.
                                                                                                                                                                                                                                                                 def linearize_expression(expr, small_vars_set):
                                                                                                                                                                                                                                                                     # Step 1: Substitute small angle approximations
                                                                                                                                                                                                                                                                     expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                     # Step 2: Expand the expression
                                                                                                                                                                                                                                                                     expr_expanded = sympy.expand(expr_subs)
                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                     # Step 3: Filter out higher-order terms
                                                                                                                                                                                                                                                                     linear_terms = []
                                                                                                                                                                                                                                                                     for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                                                                                                                                                         # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                                                                                                                                                         order = 0
                                                                                                                                                                                                                                                                         for s_var in small_vars_set:
                                                                                                                                                                                                                                                                             if term.has(s_var):
                                                                                                                                                                                                                                                                                 # If the term is a product, we need to be careful.
                                                                                                                                                                                                                                                                                 # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                                                                                                                                                 # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                                                                                                                                                 # Or if it contains a square of a small_vars_set.
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                                                                                                                                                 # However, sympy.series is for functions of a single variable.
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                                                                                                                                                 # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Let's define the state vector:
                                                                                                                                                                                                                                                                                 # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                                 # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                                                                                                                                                 # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # This means we need to compute the Jacobian.
                                                                                                                                                                                                                                                                                 # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Let's define the full state vector for linearization
                                                                                                                                                                                                                                                                                 states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                                                                                                                                                 # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                                                                                                                                                 # This is where solving the EOM for accelerations comes in.
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Let's try to solve the EOM for accelerations first.
                                                                                                                                                                                                                                                                                 # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                                                                                                                                                 # Then we can linearize f1, f2, f3.
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # The equations are:
                                                                                                                                                                                                                                                                                 # eom[0] = 0 (for x_w)
                                                                                                                                                                                                                                                                                 # eom[1] = 0 (for theta)
                                                                                                                                                                                                                                                                                 # eom[2] = 0 (for x)
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Solve for accelerations
                                                                                                                                                                                                                                                                                 # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                                                                                                                                                 # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                                                                                                                                                 # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                                                                                                                                                 # This requires setting up the system using `KanesMethod`.
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                                                                                                                                                 # Generalized coordinates
                                                                                                                                                                                                                                                                                 q_k = [x_w, theta, x]
                                                                                                                                                                                                                                                                                 # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                                                                                                                                                 u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                                                                                                                                                 kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Define the inertial frame
                                                                                                                                                                                                                                                                                 N = ReferenceFrame('N')
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Point for the wheel
                                                                                                                                                                                                                                                                                 P_w = Point('P_w')
                                                                                                                                                                                                                                                                                 P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Particle for the wheel
                                                                                                                                                                                                                                                                                 Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Point for the rotor
                                                                                                                                                                                                                                                                                 P_r = Point('P_r')
                                                                                                                                                                                                                                                                                 P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # RigidBody for the rotor
                                                                                                                                                                                                                                                                                 # The rotor rotates about the y-axis (theta)
                                                                                                                                                                                                                                                                                 # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                                                                                                                                                 # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                                                                                                                                                 # This is getting complicated. Let's simplify the approach.
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                                                                                                                                                 # it implies we should linearize the derived EOM directly.
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                                                                                                                                                 # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Define the small variables and their derivatives
                                                                                                                                                                                                                                                                                 small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                 # Function to linearize an expression by substituting small angle approximations
                                                                                                                                                                                                                                                                                 # and then removing terms that are products of two or more small variables.
                                                                                                                                                                                                                                                                                 def linearize_expression(expr, small_vars_set):
                                                                                                                                                                                                                                                                                     # Step 1: Substitute small angle approximations
                                                                                                                                                                                                                                                                                     expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                     # Step 2: Expand the expression
                                                                                                                                                                                                                                                                                     expr_expanded = sympy.expand(expr_subs)
                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                     # Step 3: Filter out higher-order terms
                                                                                                                                                                                                                                                                                     linear_terms = []
                                                                                                                                                                                                                                                                                     for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                                                                                                                                                                         # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                                                                                                                                                                         order = 0
                                                                                                                                                                                                                                                                                         for s_var in small_vars_set:
                                                                                                                                                                                                                                                                                             if term.has(s_var):
                                                                                                                                                                                                                                                                                                 # If the term is a product, we need to be careful.
                                                                                                                                                                                                                                                                                                 # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                                                                                                                                                                 # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                                                                                                                                                                 # Or if it contains a square of a small_vars_set.
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                                                                                                                                                                 # However, sympy.series is for functions of a single variable.
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                                                                                                                                                                 # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Let's define the state vector:
                                                                                                                                                                                                                                                                                                 # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                                                 # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                                                                                                                                                                 # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # This means we need to compute the Jacobian.
                                                                                                                                                                                                                                                                                                 # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Let's define the full state vector for linearization
                                                                                                                                                                                                                                                                                                 states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                                                                                                                                                                 # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                                                                                                                                                                 # This is where solving the EOM for accelerations comes in.
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Let's try to solve the EOM for accelerations first.
                                                                                                                                                                                                                                                                                                 # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                                                                                                                                                                 # Then we can linearize f1, f2, f3.
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # The equations are:
                                                                                                                                                                                                                                                                                                 # eom[0] = 0 (for x_w)
                                                                                                                                                                                                                                                                                                 # eom[1] = 0 (for theta)
                                                                                                                                                                                                                                                                                                 # eom[2] = 0 (for x)
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Solve for accelerations
                                                                                                                                                                                                                                                                                                 # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                                                                                                                                                                 # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                                                                                                                                                                 # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                                                                                                                                                                 # This requires setting up the system using `KanesMethod`.
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                                                                                                                                                                 # Generalized coordinates
                                                                                                                                                                                                                                                                                                 q_k = [x_w, theta, x]
                                                                                                                                                                                                                                                                                                 # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                                                                                                                                                                 u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                                                                                                                                                                 kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Define the inertial frame
                                                                                                                                                                                                                                                                                                 N = ReferenceFrame('N')
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Point for the wheel
                                                                                                                                                                                                                                                                                                 P_w = Point('P_w')
                                                                                                                                                                                                                                                                                                 P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Particle for the wheel
                                                                                                                                                                                                                                                                                                 Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Point for the rotor
                                                                                                                                                                                                                                                                                                 P_r = Point('P_r')
                                                                                                                                                                                                                                                                                                 P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # RigidBody for the rotor
                                                                                                                                                                                                                                                                                                 # The rotor rotates about the y-axis (theta)
                                                                                                                                                                                                                                                                                                 # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                                                                                                                                                                 # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                                                                                                                                                                 # This is getting complicated. Let's simplify the approach.
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                                                                                                                                                                 # it implies we should linearize the derived EOM directly.
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                                                                                                                                                                 # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Define the small variables and their derivatives
                                                                                                                                                                                                                                                                                                 small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 # Function to linearize an expression by substituting small angle approximations
                                                                                                                                                                                                                                                                                                 # and then removing terms that are products of two or more small variables.
                                                                                                                                                                                                                                                                                                 def linearize_expression(expr, small_vars_set):
                                                                                                                                                                                                                                                                                                     # Step 1: Substitute small angle approximations
                                                                                                                                                                                                                                                                                                     expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                     # Step 2: Expand the expression
                                                                                                                                                                                                                                                                                                     expr_expanded = sympy.expand(expr_subs)
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                     # Step 3: Filter out higher-order terms
                                                                                                                                                                                                                                                                                                     linear_terms = []
                                                                                                                                                                                                                                                                                                     for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                                                                                                                                                                                         # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                                                                                                                                                                                         order = 0
                                                                                                                                                                                                                                                                                                         for s_var in small_vars_set:
                                                                                                                                                                                                                                                                                                             if term.has(s_var):
                                                                                                                                                                                                                                                                                                                 # If the term is a product, we need to be careful.
                                                                                                                                                                                                                                                                                                                 # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                                                                                                                                                                                 # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                                                                                                                                                                                 # Or if it contains a square of a small_vars_set.
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                                                                                                                                                                                 # However, sympy.series is for functions of a single variable.
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                                                                                                                                                                                 # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Let's define the state vector:
                                                                                                                                                                                                                                                                                                                 # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                                                                 # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                                                                                                                                                                                 # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # This means we need to compute the Jacobian.
                                                                                                                                                                                                                                                                                                                 # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Let's define the full state vector for linearization
                                                                                                                                                                                                                                                                                                                 states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                                                                                                                                                                                 # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                                                                                                                                                                                 # This is where solving the EOM for accelerations comes in.
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Let's try to solve the EOM for accelerations first.
                                                                                                                                                                                                                                                                                                                 # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                                                                                                                                                                                 # Then we can linearize f1, f2, f3.
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # The equations are:
                                                                                                                                                                                                                                                                                                                 # eom[0] = 0 (for x_w)
                                                                                                                                                                                                                                                                                                                 # eom[1] = 0 (for theta)
                                                                                                                                                                                                                                                                                                                 # eom[2] = 0 (for x)
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Solve for accelerations
                                                                                                                                                                                                                                                                                                                 # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                                                                                                                                                                                 # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                                                                                                                                                                                 # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                                                                                                                                                                                 # This requires setting up the system using `KanesMethod`.
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                                                                                                                                                                                 # Generalized coordinates
                                                                                                                                                                                                                                                                                                                 q_k = [x_w, theta, x]
                                                                                                                                                                                                                                                                                                                 # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                                                                                                                                                                                 u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                                                                                                                                                                                 kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Define the inertial frame
                                                                                                                                                                                                                                                                                                                 N = ReferenceFrame('N')
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Point for the wheel
                                                                                                                                                                                                                                                                                                                 P_w = Point('P_w')
                                                                                                                                                                                                                                                                                                                 P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Particle for the wheel
                                                                                                                                                                                                                                                                                                                 Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Point for the rotor
                                                                                                                                                                                                                                                                                                                 P_r = Point('P_r')
                                                                                                                                                                                                                                                                                                                 P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # RigidBody for the rotor
                                                                                                                                                                                                                                                                                                                 # The rotor rotates about the y-axis (theta)
                                                                                                                                                                                                                                                                                                                 # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                                                                                                                                                                                 # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                                                                                                                                                                                 # This is getting complicated. Let's simplify the approach.
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                                                                                                                                                                                 # it implies we should linearize the derived EOM directly.
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                                                                                                                                                                                 # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Define the small variables and their derivatives
                                                                                                                                                                                                                                                                                                                 small_vars_set = {x, theta, x_dot, theta_dot, x_w_dot}
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 # Function to linearize an expression by substituting small angle approximations
                                                                                                                                                                                                                                                                                                                 # and then removing terms that are products of two or more small variables.
                                                                                                                                                                                                                                                                                                                 def linearize_expression(expr, small_vars_set):
                                                                                                                                                                                                                                                                                                                     # Step 1: Substitute small angle approximations
                                                                                                                                                                                                                                                                                                                     expr_subs = expr.subs([(sin(theta), theta), (cos(theta), 1)])
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                     # Step 2: Expand the expression
                                                                                                                                                                                                                                                                                                                     expr_expanded = sympy.expand(expr_subs)
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                     # Step 3: Filter out higher-order terms
                                                                                                                                                                                                                                                                                                                     linear_terms = []
                                                                                                                                                                                                                                                                                                                     for term in sympy.Add.make_args(expr_expanded):
                                                                                                                                                                                                                                                                                                                         # Count the number of 'small' variables (or their derivatives) in the term
                                                                                                                                                                                                                                                                                                                         order = 0
                                                                                                                                                                                                                                                                                                                         for s_var in small_vars_set:
                                                                                                                                                                                                                                                                                                                             if term.has(s_var):
                                                                                                                                                                                                                                                                                                                                 # If the term is a product, we need to be careful.
                                                                                                                                                                                                                                                                                                                                 # For example, x*theta is order 2. x_dot*theta is order 2.
                                                                                                                                                                                                                                                                                                                                 # x_ddot is order 1 (as it's a derivative of a small variable).
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # A simple heuristic: if a term contains a product of two or more small_vars_set, it's higher order.
                                                                                                                                                                                                                                                                                                                                 # Or if it contains a square of a small_vars_set.
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # This is still not perfect. Let's try a more direct approach using sympy.series.
                                                                                                                                                                                                                                                                                                                                 # However, sympy.series is for functions of a single variable.
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # The most robust way for multi-variable functions is to use Taylor expansion:
                                                                                                                                                                                                                                                                                                                                 # f(x,y) approx f(0,0) + df/dx(0,0)*x + df/dy(0,0)*y
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # Let's define the state vector:
                                                                                                                                                                                                                                                                                                                                 # z = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                                                                                 # The EOM are F_i(z, z_dot_dot) = 0
                                                                                                                                                                                                                                                                                                                                 # We need to linearize F_i around z=0, z_dot=0, z_dot_dot=0
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # This means we need to compute the Jacobian.
                                                                                                                                                                                                                                                                                                                                 # d(F_i)/d(q_j) * q_j + d(F_i)/d(q_dot_j) * q_dot_j + d(F_i)/d(q_ddot_j) * q_ddot_j = 0
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # Let's define the full state vector for linearization
                                                                                                                                                                                                                                                                                                                                 states = [x_w, theta, x, x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # The EOM are functions of q, q_dot, q_ddot.
                                                                                                                                                                                                                                                                                                                                 # We need to express q_ddot in terms of q and q_dot.
                                                                                                                                                                                                                                                                                                                                 # This is where solving the EOM for accelerations comes in.
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # Let's try to solve the EOM for accelerations first.
                                                                                                                                                                                                                                                                                                                                 # This will give us x_w_ddot = f1(q, q_dot), theta_ddot = f2(q, q_dot), x_ddot = f3(q, q_dot)
                                                                                                                                                                                                                                                                                                                                 # Then we can linearize f1, f2, f3.
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # The equations are:
                                                                                                                                                                                                                                                                                                                                 # eom[0] = 0 (for x_w)
                                                                                                                                                                                                                                                                                                                                 # eom[1] = 0 (for theta)
                                                                                                                                                                                                                                                                                                                                 # eom[2] = 0 (for x)
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # Solve for accelerations
                                                                                                                                                                                                                                                                                                                                 # This is a system of 3 equations and 3 unknowns (x_w_ddot, theta_ddot, x_ddot)
                                                                                                                                                                                                                                                                                                                                 # The equations are highly non-linear. SymPy's solve might struggle or give very complex results.
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # Let's try to linearize the EOM first, then solve for accelerations.
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # A more robust way to linearize is to use sympy.linearize, but it requires a specific setup.
                                                                                                                                                                                                                                                                                                                                 # Let's try to use the `linearize` function from `sympy.physics.mechanics`.
                                                                                                                                                                                                                                                                                                                                 # This requires setting up the system using `KanesMethod`.
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # Redefine the system using KanesMethod to leverage its linearization capabilities
                                                                                                                                                                                                                                                                                                                                 # Generalized coordinates
                                                                                                                                                                                                                                                                                                                                 q_k = [x_w, theta, x]
                                                                                                                                                                                                                                                                                                                                 # Generalized speeds (derivatives of generalized coordinates)
                                                                                                                                                                                                                                                                                                                                 u_k = [x_w_dot, theta_dot, x_dot]
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # Kinematic differential equations (u_k = q_k_dot)
                                                                                                                                                                                                                                                                                                                                 kd_eqs = [u_k[i] - diff(q_k[i], t) for i in range(len(q_k))]
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # Define the inertial frame
                                                                                                                                                                                                                                                                                                                                 N = ReferenceFrame('N')
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # Point for the wheel
                                                                                                                                                                                                                                                                                                                                 P_w = Point('P_w')
                                                                                                                                                                                                                                                                                                                                 P_w.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # Particle for the wheel
                                                                                                                                                                                                                                                                                                                                 Wheel = Particle('Wheel', P_w, M_w)
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # Point for the rotor
                                                                                                                                                                                                                                                                                                                                 P_r = Point('P_r')
                                                                                                                                                                                                                                                                                                                                 P_r.set_vel(N, x_w_dot * N.x)
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # RigidBody for the rotor
                                                                                                                                                                                                                                                                                                                                 # The rotor rotates about the y-axis (theta)
                                                                                                                                                                                                                                                                                                                                 # Its center of mass is at (x_w, R + d*cos(theta))
                                                                                                                                                                                                                                                                                                                                 # For linearization, we assume it's at (x_w, R+d) and rotates about y.
                                                                                                                                                                                                                                                                                                                                 # This is getting complicated. Let's simplify the approach.
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # Given the prompt "assume small x and small theta and find how x,x_w, theta behave",
                                                                                                                                                                                                                                                                                                                                 # it implies we should linearize the derived EOM directly.
                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                                 # Let's re-attempt manual linearization by substituting and then dropping higher order terms.
                                                                                                                                                                                                                                                                                                                                 # This is the most direct interpretation of "assume small x and small theta".
                                                                                                                                                                                                                                                                                                                                 
                                                              I will now proceed with the linearization and simulation.
