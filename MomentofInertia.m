function [Ilw,Irw,rho] = MomentofInertia(c,filename, Tr, Tt, ct, cr, b, m, lam,nacan)

% c: chord length
% Tr: thickness at root of wing
% Tt: thickness at tip of wing
% ct: chord length at tip of wing
% cr: chord length at root of wing
% b: span of ONE wing on aircraft
% m: mass of wing
% lam: wing sweep angle [deg]
% nacan: NACA number of airfoil

% Must adjust range of excel file, dependent on NACA # and number of points

xu = readmatrix(filename,'Range','A1:A66'); %Importing airfoil coordinates 
yu = readmatrix(filename,'Range','B1:B66'); %Importing airfoil coordinates

xcu = xu/c;

mdl = fittype({'sqrt(x)','x','x^2','x^3','x^4'}, 'coefficients', {'a0','a1','a2','a3','a4'}); %form of line

fittedmdl = fit(xcu,yu,mdl); %Finding coefficients
coeff = coeffvalues(fittedmdl);

a0u = coeff(1);
a1u = coeff(2);
a2u = coeff(3);
a3u = coeff(4);
a4u = coeff(5);

y0u = (a0u * sqrt(xcu)) + (a1u*xcu) + (a2u*xcu.^2) + (a3u * xcu.^3) + (a4u * xcu.^4);

%Ensuring line of best fit is appropriate 

figure
plot(xu,yu)
hold on 
plot(xu,y0u)
title(['NACA', nacan, 'Upper Half'])
xlabel('Chord Length')
ylabel('Camber')
legend('NACA Coordinates','Line of Best Fit')


% Fitting bottom half of the airfoil

xb = readmatrix(filename,'Range','A67:A132'); %Importing airfoil coordinates
yb = readmatrix(filename,'Range','B67:B132'); %Importing airfoil coordinates

xcb = xb/c;

mdl = fittype({'sqrt(x)','x','x^2','x^3','x^4'}, 'coefficients', {'a0','a1','a2','a3','a4'}); %form of line

fittedmdl = fit(xcb,yb,mdl); %Finding coefficients
coeff = coeffvalues(fittedmdl);


% Absolute value is added for a symmetric airfoil, remove if your airfoil
% is not symmetric
a0b = abs(coeff(1));
a1b = abs(coeff(2));
a2b = abs(coeff(3));
a3b = abs(coeff(4));
a4b = abs(coeff(5));

y0b = (a0b * sqrt(xcb)) + (a1b*xcb) + (a2b*xcb.^2) + (a3b * xcb.^3) + (a4b * xcb.^4);


%Ensuring line of best fit is appropriate 

figure
plot(xb,yb)
hold on 
plot(xb,y0b)
title(['NACA', nacan, 'Bottom Half'])
xlabel('Chord Length')
ylabel('Camber')
legend('NACA Coordinates','Line of Best Fit')

%Now averaging the bottom and the top to find the overall line of best fit
a0 = (a0u + a0b) / 2;
a1 = (a1u + a1b) / 2;
a2 = (a2u + a2b) / 2;
a3 = (a3u + a3b) / 2;
a4 = (a4u + a4b) / 2;

x = [xu; xb];
y = [yu; yb];
xc = x / c;
y0 = (a0b * sqrt(xc)) + (a1b*xc) + (a2b*xc.^2) + (a3b * xc.^3) + (a4b * xc.^4);

figure
plot(xc,y)
hold on
plot(xc,y0)
title('NACA', nacan)
xlabel('Chord Length')
ylabel('Camber')
legend('NACA Coordinates','Line of Best Fit')

% Calculating volume of wing 

v0 = (1/60) * ((40*a0) + (30*a1) + (20*a2)+ (15*a3) + (12*a4));

Ka = (Tr * ((3*(cr^2)) + (2*cr*ct) + (ct^2))) + (Tt * ((cr^2) + (2*cr*ct) + (3*(ct^2)))); %Wing planform coefficient

V = (b/12) * Ka * v0; %Volume of Wing [m^3]
rho = m / V; %Calculating the density of the wing

% Calculating Moments

deltar = 1; %Right wing
deltal = -1; %Left wing

Kb = (Tr * ((4*(cr^3))+ (3*(cr^2)*ct) + (2*cr*(ct^2)) + (ct^3))) + (Tt * ((cr^3) + (2*(cr^2)*ct)...
    + (3*cr*(ct^2)) + (4*(ct^3))));
Kc = (Tr * ((3*(cr^2)) + (4*cr*ct) + (3*(ct^2)))) + ((2*Tt) * ((cr^2) + (3*cr*ct) + (6*(ct^2))));

v1 = (1/60) * ((56*a0) + (50*a1) + (40*a2) + (33*a3) + (28*a4));

Myz = -rho * (b/240) * ((3*Kb*v1)+ (4*b*Kc*v0*tand(lam))); %[N m]
Mxz = (rho * deltar * (b^2) * Kc * v0) / 60;
Mxy = 0;

% Finding center of gravity with respect to wing location

x_bar = -1*(((3*Kb*v1) + (4*b*Kc*v0*tand(lam))) / (20*Ka*v0));
y_bar = (deltar * b * Kc) / (5*Ka);
z_bar = 0;

% Finding moment of inertia matrix for right wing

Kd = (Tr * (cr+ct) * ((2*(cr^2)) + (cr*ct) + (2*(ct^2)))) + (Tt*((cr^3) + (3*(cr^2)*ct) + (6*cr*(ct^2))...
    + (10*(ct^3))));
Ke = (Tr * ((5*(cr^4)) + (4*(cr^3)*ct)) + (3*(cr^2)*(ct^2)) + (2*cr*(ct^3)) + (ct^4))...
    + (Tt * ((cr^4) + (2*(cr^3)*ct) + (3*(cr^2)*(ct^2)) + (4*cr*(ct^3)) + (5*(ct^4))));


Kf = (Tr * ((cr^2) + (2*cr*ct) + (2*(ct^2)))) + (Tt * ((cr^2) + (4*cr*ct) + (10*(ct^2))));


Kg = ((Tr^3) * ((35*(cr^4)) + (20*(cr^3)*ct) + (10*(cr^2)*(ct^2)) + (4*cr*(ct^3)) + (ct^4)))...
    + (((Tr^2)*Tt) * ((15*(cr^4))...
    + (20*(cr^3)*ct) + (18*(cr^2)*(ct^2)) + (12*cr*(ct^3)) + (5*(ct^4)))) + ((Tr*(Tt^2)) * ((5*(cr^4))...
    + (12*(cr^3)*ct)...
    + (18*(cr^2)*(ct^2)) + (20*cr*(ct^3)) + (15*(ct^4)))) + ((Tt^3) * ((cr^4) + (4*(cr^3)*ct)...
    + (10*(cr^2)*(ct^2))...
    + (20*cr*(ct^3)) + (35*(ct^4))));

v2 = (1/980) * ((856*a0) + (770*a1) + (644*a2) + (553*a3) + (484*a4));
v3 = (2/5 * (a0^3)) + ((a0^2)*a1) + ((3*(a0^2)*a2)/4) + ((3*(a0^2)*a3)/5) + (((a0^2)*a4)/2) + ((6*a0*(a1^2))/7)...
    + ((4*a0*a1*a2)/3) + ((12*a0*a1*a3)/11) + ((12*a0*a1*a4)/13) + ((6*a0*(a2^2))/11) + ((12*a0*a2*a3)/13)...
    + ((4*a0*a2*a4)/5) + ((2*a0*(a3^2))/5) + ((12*a0*a3*a4)/17) + ((6*a0*(a4^2))/19) + ((a1^3)/4)...
    + ((3*(a1^2)*a2)/5) + (((a1^2)*a3)/2) + ((3*(a1^2)*a4)/7) + ((a1*(a2^2))/2) + ((6*a1*a2*a3)/7)...
    + ((3*a1*a2*a4)/4) + ((3*a1*(a3^2))/8) + ((2*a1*a3*a4)/3) + ((3*a1*(a4^2))/10) + ((a2^3)/7)...
    + ((3*(a2^2)*a3)/8) + (((a2^2)*a4)/3) + ((a2*(a3^2))/3) + ((3*a2*a3*a4)/5) + ((3*a2*(a4^2))/11)...
    + ((a3^3)/10) + ((3*(a3^2)*a4)/11) + ((a3*(a4^2))/4) + ((a4^3)/13);



Ixx_rw = ((rho * b)/3360) * ((56*(b^2)*Kf*v0) + (Kg*v3));
Iyy_rw = ((rho * b)/10080) * ((84*b*v1*(2*b*Kf*(tand(lam))^2 + (Kd*tand(lam)))) + (49*Ke*v2) + (3*Kg*v3));
Izz_rw = ((rho*b)/1440) * ((12*b) * (2*b*Kf*v0*((tand(lam))^2 + 1)) + (Kd*v1*tand(lam)) + (7*Ke*v2));
Ixy_rw = ((-rho*deltar*(b^2))/240) * ((4*b*Kf*v0*tand(lam)) + (Kd*v1));
Iyx_rw = Ixy_rw;
Ixz = 0;
Izx = 0;
Iyz = 0;
Izy = 0;

I0rw = [Ixx_rw, -Ixy_rw, -Ixz; -Iyx_rw, Iyy_rw, -Iyz; -Izx, -Izy, Izz_rw];

Irw = I0rw + m* [(y_bar^2 + z_bar^2), (-x_bar*y_bar), (-x_bar*z_bar); ...
    (-x_bar*y_bar), (x_bar^2 + z_bar^2), (-y_bar*z_bar);...
    (-x_bar*z_bar), (-y_bar*z_bar), (x_bar^2 + y_bar^2)];

% Calculating Moment of Inertia for Left Wing

Mxz = (rho * deltal * (b^2) * Kc * v0) / 60;
y_bar = (deltal * b * Kc) / (5*Ka);
Ixy_lw = ((-rho*deltal*(b^2))/240) * ((4*b*Kf*v0*tand(lam)) + (Kd*v1));
Ixx_lw = Ixx_rw;
Iyx_lw = Ixy_lw;
Iyy_lw = Iyy_rw;
Izz_lw = Izz_rw;


I0lw = [Ixx_lw, -Ixy_lw, -Ixz; -Iyx_rw, Iyy_rw, -Iyz; -Izx, -Izy, Izz_rw];
Ilw = I0lw + m* [(y_bar^2 + z_bar^2), (-x_bar*y_bar), (-x_bar*z_bar); ...
    (-x_bar*y_bar), (x_bar^2 + z_bar^2), (-y_bar*z_bar);...
    (-x_bar*z_bar), (-y_bar*z_bar), (x_bar^2 + y_bar^2)];

end