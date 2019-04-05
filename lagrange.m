syms x(t) y(t) a(t) phi1(t) phi2(t) OC e Ic I0 m mw len hgt kx ky bx by...
    g x_dot a_dot y_dot phi1_dot phi2_dot xc1 xc2 yc1 yc2 xc1_dot xc2_dot yc1_dot yc2_dot %hgt - od srodka do dolnej krawedzi(w osi y) len - od srodka do osi sprezyny(w osi x)

x_dot = diff(x, t)
y_dot = diff(y, t)
a_dot = diff(a, t)
phi1_dot = diff(phi1, t)
phi2_dot = diff(phi2, t)
xc1 = -OC*cos((pi/3)-a)-e*cos(phi1)+x
yc1 = OC*sin((pi/3)-a)+e*sin(phi1)+y
xc1_dot = diff(xc1, t)
yc1_dot = diff(yc1, t)
xc2 = OC*cos((pi/3)-a)+e*cos(phi2)+x
yc2 = OC*sin((pi/3)-a)+e*sin(phi2)+y
xc2_dot = diff(xc2, t)
yc2_dot = diff(yc2, t)

T = 0.5*m*(x_dot^2+y_dot^2) + 0.5*Ic*a_dot^2 + 0.5*mw*xc1_dot^2 + 0.5*mw*yc1_dot^2 + 0.5*I0*phi1_dot^2+...
    0.5*mw*xc2_dot^2 + 0.5*mw*yc2_dot^2 + 0.5*I0*phi2_dot^2

V = m*g*y + mw*g*yc1 + mw*g*yc2 + 0.5*kx*(-x*hgt*a*x)^2 + 0.5*ky*(y+len*a)^2
L = T-V
N = 0.5*bx*(-x_dot-hgt*a_dot)^2+0.5*by*(y_dot+len*a_dot)^2

dldx = functionalDerivative(L, x)
dldy = functionalDerivative(L, y)
dlda = functionalDerivative(L, a)
dldphi1 = functionalDerivative(L, phi1)
dldphi2 = functionalDerivative(L, phi2)
syms xvar yvar avar phi1var phi2var
syms x_dotvar y_dotvar a_dotvar phi1_dotvar phi2_dotvar
Lvar = subs(L, x, xvar)
Lvar = subs(Lvar, y, yvar)
Lvar = subs(Lvar, a, avar)
Lvar = subs(Lvar, phi1, phi1var)
Lvar = subs(Lvar, phi2, phi2var)
Lvar = subs(Lvar, x_dot, x_dotvar)
Lvar = subs(Lvar, y_dot, y_dotvar)
Lvar = subs(Lvar, phi1_dot, phi1_dotvar)
Lvar = subs(Lvar, phi2_dot, phi2_dotvar)
Nvar = subs(N, x, xvar)
Nvar = subs(Nvar, y, yvar)
Nvar = subs(Nvar, a, avar)
Nvar = subs(Nvar, phi1, phi1var)
Nvar = subs(Nvar, phi2, phi2var)
Nvar = subs(Nvar, x_dot, x_dotvar)
Nvar = subs(Nvar, y_dot, y_dotvar)
Nvar = subs(Nvar, phi1_dot, phi1_dotvar)
Nvar = subs(Nvar, phi2_dot, phi2_dotvar)
%
dldx_dotvar = diff(Lvar, x_dotvar)
dldy_dotvar = diff(Lvar, y_dotvar)
dlda_dotvar = diff(Lvar, a_dotvar)
dldphi1_dotvar = diff(Lvar, phi1_dotvar)
dldphi2_dotvar = diff(Lvar, phi2_dotvar)
dndx_dotvar = diff(Nvar, x_dotvar)
dndy_dotvar = diff(Nvar, y_dotvar)
dnda_dotvar = diff(Nvar, a_dotvar)
dndphi1_dotvar = diff(Nvar, phi1_dotvar)
dndphi2_dotvar = diff(Nvar, phi2_dotvar)
%zmian dldx_dotvar na dldx_dot
dldx_dot = subs(dldx_dotvar, xvar, x)
dldx_dot = subs(dldx_dot, yvar, y)
dldx_dot = subs(dldx_dot, avar, a)
dldx_dot = subs(dldx_dot, phi1var, phi1)
dldx_dot = subs(dldx_dot, phi2var, phi2)
dldx_dot = subs(dldx_dot, x_dotvar, x_dot)
dldx_dot = subs(dldx_dot, y_dotvar, y_dot)
dldx_dot = subs(dldx_dot, phi1_dotvar, phi1_dot)
dldx_dot = subs(dldx_dot, phi2_dotvar, phi2_dot)

%zmiana dldy_dotvar na dldy_dot
dldy_dot = subs(dldy_dotvar, xvar, x)
dldy_dot = subs(dldy_dot, yvar, y)
dldy_dot = subs(dldy_dot, avar, a)
dldy_dot = subs(dldy_dot, phi1var, phi1)
dldy_dot = subs(dldy_dot, phi2var, phi2)
dldy_dot = subs(dldy_dot, x_dotvar, x_dot)
dldy_dot = subs(dldy_dot, y_dotvar, y_dot)
dldy_dot = subs(dldy_dot, phi1_dotvar, phi1_dot)
dldy_dot = subs(dldy_dot, phi2_dotvar, phi2_dot)

%zmiana dlda_dotvar na dlda_dot
dlda_dot = subs(dlda_dotvar, xvar, x)
dlda_dot = subs(dlda_dot, yvar, y)
dlda_dot = subs(dlda_dot, avar, a)
dlda_dot = subs(dlda_dot, phi1var, phi1)
dlda_dot = subs(dlda_dot, phi2var, phi2)
dlda_dot = subs(dlda_dot, x_dotvar, x_dot)
dlda_dot = subs(dlda_dot, y_dotvar, y_dot)
dlda_dot = subs(dlda_dot, phi1_dotvar, phi1_dot)
dlda_dot = subs(dlda_dot, phi2_dotvar, phi2_dot)

%zmiana dldphi1_dotvar na dldphi1_dot
dldphi1_dot = subs(dldphi1_dotvar, xvar, x)
dldphi1_dot = subs(dldphi1_dot, yvar, y)
dldphi1_dot = subs(dldphi1_dot, avar, a)
dldphi1_dot = subs(dldphi1_dot, phi1var, phi1)
dldphi1_dot = subs(dldphi1_dot, phi2var, phi2)
dldphi1_dot = subs(dldphi1_dot, x_dotvar, x_dot)
dldphi1_dot = subs(dldphi1_dot, y_dotvar, y_dot)
dldphi1_dot = subs(dldphi1_dot, phi1_dotvar, phi1_dot)
dldphi1_dot = subs(dldphi1_dot, phi2_dotvar, phi2_dot)


%zmiana dldphi2_dotvar na dldphi2_dot
dldphi2_dot = subs(dldphi2_dotvar, xvar, x)
dldphi2_dot = subs(dldphi2_dot, yvar, y)
dldphi2_dot = subs(dldphi2_dot, avar, a)
dldphi2_dot = subs(dldphi2_dot, phi1var, phi1)
dldphi2_dot = subs(dldphi2_dot, phi2var, phi2)
dldphi2_dot = subs(dldphi2_dot, x_dotvar, x_dot)
dldphi2_dot = subs(dldphi2_dot, y_dotvar, y_dot)
dldphi2_dot = subs(dldphi2_dot, phi1_dotvar, phi1_dot)
dldphi2_dot = subs(dldphi2_dot, phi2_dotvar, phi2_dot)
%zmian dndx_dotvar na dndx_dot
dndx_dot = subs(dndx_dotvar, xvar, x)
dndx_dot = subs(dndx_dot, yvar, y)
dndx_dot = subs(dndx_dot, avar, a)
dndx_dot = subs(dndx_dot, phi1var, phi1)
dndx_dot = subs(dndx_dot, phi2var, phi2)
dndx_dot = subs(dndx_dot, x_dotvar, x_dot)
dndx_dot = subs(dndx_dot, y_dotvar, y_dot)
dndx_dot = subs(dndx_dot, phi1_dotvar, phi1_dot)
dndx_dot = subs(dndx_dot, phi2_dotvar, phi2_dot)

%zmiana dndy_dotvar na dndy_dot
dndy_dot = subs(dndy_dotvar, xvar, x)
dndy_dot = subs(dndy_dot, yvar, y)
dndy_dot = subs(dndy_dot, avar, a)
dndy_dot = subs(dndy_dot, phi1var, phi1)
dndy_dot = subs(dndy_dot, phi2var, phi2)
dndy_dot = subs(dndy_dot, x_dotvar, x_dot)
dndy_dot = subs(dndy_dot, y_dotvar, y_dot)
dndy_dot = subs(dndy_dot, phi1_dotvar, phi1_dot)
dndy_dot = subs(dndy_dot, phi2_dotvar, phi2_dot)

%zmiana dnda_dotvar na dnda_dot
dnda_dot = subs(dnda_dotvar, xvar, x)
dnda_dot = subs(dnda_dot, yvar, y)
dnda_dot = subs(dnda_dot, avar, a)
dnda_dot = subs(dnda_dot, phi1var, phi1)
dnda_dot = subs(dnda_dot, phi2var, phi2)
dnda_dot = subs(dnda_dot, x_dotvar, x_dot)
dnda_dot = subs(dnda_dot, y_dotvar, y_dot)
dnda_dot = subs(dnda_dot, phi1_dotvar, phi1_dot)
dnda_dot = subs(dnda_dot, phi2_dotvar, phi2_dot)

%zmiana dndphi1_dotvar na dndphi1_dot
dndphi1_dot = subs(dndphi1_dotvar, xvar, x)
dndphi1_dot = subs(dndphi1_dot, yvar, y)
dndphi1_dot = subs(dndphi1_dot, avar, a)
dndphi1_dot = subs(dndphi1_dot, phi1var, phi1)
dndphi1_dot = subs(dndphi1_dot, phi2var, phi2)
dndphi1_dot = subs(dndphi1_dot, x_dotvar, x_dot)
dndphi1_dot = subs(dndphi1_dot, y_dotvar, y_dot)
dndphi1_dot = subs(dndphi1_dot, phi1_dotvar, phi1_dot)
dndphi1_dot = subs(dndphi1_dot, phi2_dotvar, phi2_dot)


%zmiana dndphi2_dotvar na dndphi2_dot
dndphi2_dot = subs(dndphi2_dotvar, xvar, x)
dndphi2_dot = subs(dndphi2_dot, yvar, y)
dndphi2_dot = subs(dndphi2_dot, avar, a)
dndphi2_dot = subs(dndphi2_dot, phi1var, phi1)
dndphi2_dot = subs(dndphi2_dot, phi2var, phi2)
dndphi2_dot = subs(dndphi2_dot, x_dotvar, x_dot)
dndphi2_dot = subs(dndphi2_dot, y_dotvar, y_dot)
dndphi2_dot = subs(dndphi2_dot, phi1_dotvar, phi1_dot)
dndphi2_dot = subs(dndphi2_dot, phi2_dotvar, phi2_dot)
ddt_dlddqi = [diff(dldx_dot, t); diff(dldy_dot, t); diff(dlda_dot, t); diff(dldphi1_dot, t); diff(dldphi2_dot, t)]
dldqi = [dldx; dldy; dlda; dldphi1; dldphi2]
dnddqi = [diff(dndx_dot, t); diff(dndy_dot, t); diff(dnda_dot, t); diff(dndphi1_dot, t); diff(dndphi2_dot, t)]
m_eqns = -ddt_dlddqi+dldqi-dnddqi

dsolve(m_eqns ==0)
