function [dA] = modelv1(t, J)


m = 11.53; %kg - masa rynny
mw = 3.86; %kg - masa wibratora
len = 0.5; %m - 1/2 odleglosc miedzy sprezynami
hgt = 0.135; %m - odleglosc miedzy srodkiem masy a dolem rynny
kx = 628.4; %N/m - stala sprezystosci w x
ky = 3706; %N/m - stala sprezystosci w y
bx = 0.0239; %N/(m/s)tlumienie w x
by = 0.0027; %N/(m/s) tlumienie w y
I = 5; %kgm^2 moment bezwladnosci rynny
Iw = 3; %kgm^2 moment bezwlasnosci wibratora
e = 0.02; %m - odleglosc srodka ciezkosci od osi obrotu wibratora
g = 9.81; %ms^-2 - stala grawitacyjna
OC = 0.22; %odleglosc od srodka masy rynny do osi obrotu wibratora


x_dot = J(1);
x = J(2);
y_dot = J(3);
y = J(4);
a_dot = J(5);
a = J(6);
f1_dot = J(7);
f1 = J(8);
f2_dot = J(9);
f2 = J(10);




Mm = sparse(10, 10);
Mm(1, 1) = m;
Mm(2, 2) = 1;
Mm(3, 3) = m;
Mm(4, 4) = 1;
Mm(5, 5) = I;
Mm(6, 6) = 1;
Mm(7, 7) = Iw;
Mm(8, 8) = 1;
Mm(9, 9) = Iw;
Mm(10, 10) = 1;


overload = 2.5;
f_main = 50; %Hz
n_sync = 1400;
w_sync = 2*pi*n_sync/60;
n_nom = n_sync * 0.99; %RPM
w_nom = 2*pi*n_nom/60; %s^-1
P_nom = 120; %W
M_nom = P_nom/w_nom;
S_nom = (n_sync-n_nom)/n_sync;
S_stall = S_nom*(overload+sqrt((overload^2)-1));
n_stall = n_sync-S_stall*n_sync;
w_stall = 2*pi*n_stall/60;
M_stall = overload*M_nom;

m_el =@(w) (2.*M_stall.*(w_sync-w_stall).*(w_sync-w))/((w_sync-w_stall).^2+(w_sync-w).^2);


Mw = [ -kx*x-kx*hgt*a+bx*x_dot+bx*hgt*a_dot; ...
       x_dot;...
       -ky*y+2*mw*g+g*m-ky*len*a-by*y_dot-by*len*a_dot;...
       y_dot;...
       -0.5*mw*OC^2*a_dot^2*sin(pi/3+2*a)-ky*len*y-kx*hgt^2*a-ky*len^2*a+kx*hgt*x-0.5*mw*OC^2*a_dot^2*sin(pi/3+2*a)+2*mw*g*OC*sin(pi/6+a)-0.5*mw*OC*a_dot*f1_dot*e*cos(pi/6+a+f1)+mw*OC*a_dot*f2_dot*e*cos(pi/6+a+f2)+0.5*mw*OC*a_dot*f1_dot*e*cos(pi/6+a+f1)+0.5*mw*OC*a_dot*f1_dot*cos(pi/6*a*f1*e)-0.5*mw*OC*a_dot*f1_dot*cos(pi/6+a+f1)*e-bx*a_dot*hgt^2-bx*x_dot*hgt-by*a_dot*len^2-by*y_dot*len+m_el(f1)-m_el(f2);...
       a_dot;...
       -mw*g*cos(f1)*e-m_el(f1);...
       f1_dot;...
       mw*g*cos(f2)*e+m_el(f2);...
       f2_dot ];

     dA = Mm\Mw;

end

