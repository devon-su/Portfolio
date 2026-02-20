function [r, v] = gauss_IOD(L1, L2, L3, R1, R2, R3, tau1, tau3, mu)
% Gauss method with iterative improvement to calculate the state vector of an orbiting body from angles-only observations at three closely spaced times

% Difference between last and first pass
tau  = tau3 - tau1;

% Cross products among the direction cosine vectors
p1 = cross(L2,L3);
p2 = cross(L1,L3);
p3 = cross(L1,L2);

Do = dot(L1,p1);

D  = [[dot(R1,p1) dot(R1,p2) dot(R1,p3)]
      [dot(R2,p1) dot(R2,p2) dot(R2,p3)]
      [dot(R3,p1) dot(R3,p2) dot(R3,p3)]];

E = dot(R2,L2);

A = 1/Do*(-D(1,2)*tau3/tau + D(2,2) + D(3,2)*tau1/tau);
B = 1/6/Do*(D(1,2)*(tau3^2 - tau^2)*tau3/tau + D(3,2)*(tau^2 - tau1^2)*tau1/tau);

a = -(A^2 + 2*A*E + norm(R2)^2);
b = -2*mu*B*(A + E);
c = -(mu*B)^2;

% Calculate the roots of Equation using built-in MATLAB solver
Roots = roots([1 0 a 0 0 b 0 0 c]);

% Only want positive real root
x = posroot(Roots);

f1 =    1 - 1/2*mu*tau1^2/x^3;
f3 =    1 - 1/2*mu*tau3^2/x^3;

g1 = tau1 - 1/6*mu*(tau1/x)^3;
g3 = tau3 - 1/6*mu*(tau3/x)^3;

rho2 = A + mu*B/x^3;

rho1 = 1/Do*((6*(D(3,1)*tau1/tau3 + D(2,1)*tau/tau3)*x^3 + mu*D(3,1)*(tau^2 - tau1^2)*tau1/tau3)/(6*x^3 + mu*(tau^2 - tau3^2)) - D(1,1));

rho3 = 1/Do*((6*(D(1,3)*tau3/tau1 - D(2,3)*tau/tau1)*x^3 + mu*D(1,3)*(tau^2 - tau3^2)*tau3/tau1)/(6*x^3 + mu*(tau^2 - tau1^2)) - D(3,3)); 

r1 = R1 + rho1*L1;
r2 = R2 + rho2*L2;
r3 = R3 + rho3*L3;

v2 = (-f3*r1 + f1*r3)/(f1*g3 - f3*g1);

% Save r2 and v2 as initial estimates
r = r2;
v = v2;

function x = posroot(Roots)
% Positive root solver

    posroots = Roots(find(Roots>0 & ~imag(Roots)));
    npositive = length(posroots);
    if npositive == 0
        fprintf('\n\n ** There are no positive roots. \n\n')
        return
    end
    if npositive == 1
        x = posroots;
    else
        fprintf('\n\n ** There are two or more positive roots.\n')
    
    for i = 1:npositive
        fprintf('\n root #%g = %g',i,posroots(i))
    end
        fprintf('\n\n Make a choice:\n')
        nchoice = 0;
    while nchoice < 1 || nchoice > npositive
        nchoice = input(' Use root #? ');
    end
        x = posroots(nchoice);
        fprintf('\n We will use %g .\n', x)
    end
    end
end
