[![View RKN1210 - A 12th/10th order Runge-Kutta-Nyström integrator on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/25291-rkn1210-a-12th-10th-order-runge-kutta-nystrom-integrator)

[![Donate to Rody](https://i.stack.imgur.com/bneea.png)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=4M7RMVNMKAXXQ&source=url)

# FEX-RKN1210
High-order Runge Kutta Nyström integrator

RKN1210 12th/10th order Runge-Kutta-Nyström integrator
RKN1210() is a 12th/10th order variable-step numerical integrator for second-order ordinary differential equations of the form
y'' = f(t, y) (1)
with initial conditions

y(t0) = y0
y'(t0) = yp0 (2)

This second-order differential equation is integrated with a Runge-Kutta-Nyström method using 17 function evaluations per step. RKN12(10) is a very high-order method, to be used in problems with *extremely* stringent error tolerances.

The RKN-class of integrators is especially suited for problems of type (1). Compared to a classic Runge-Kutta integration scheme, the same accuracy can be obtained with fewer function evaluations. Also, it has been shown in various studies that this particular integration method is overall more efficient than (symplectic) multi-step or extrapolation methods that give the same accuracy.

RKN1210's behavior is very similar MATLAB's ODE-integrator suite; you can set options via ODESET, and input and output values are also practically the same.

Both output functions and event functions are fully supported.

The construction of RKN12(10) is described in
High-Order Embedded Runge-Kutta-Nyström Formulae
J. R. DORMAND, M. E. A. EL-MIKKAWY, AND P. J. PRINCE
IMA Journal of Numerical Analysis (1987) 7, 423-430

Coefficients obtained from
http://www.tampa.phys.ucl.ac.uk/rmat/test/rknint.f
These are also available in any format on request to these authors.

If you like this work, please consider [a donation](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=4M7RMVNMKAXXQ&source=url).
