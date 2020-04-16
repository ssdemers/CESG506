% set initial conditions
EA = 2100;
nDof = 2;
I = eye(nDof);
u = [0,0];
Lavec = [5.5,0.5];
Lbvec = [-4,0.5];
La = sqrt(dot(Lavec,Lavec));
Lb = sqrt(dot(Lbvec,Lbvec));
P = [0,0];
Ans = zeros();
iter = zeros();

% find Pcr -- refine iterations to get more accurate Pcr
for i = 1:22
    P(2) = 0.971+.0005*i;
    for it = 1:1000
        lavec = Lavec+u;
        lbvec = Lbvec+u;
        la = sqrt(dot(lavec,lavec));
        lb = sqrt(dot(lbvec,lbvec));
        na = lavec/la;
        nb = lbvec/lb;
        straina = log(la/La);
        strainb = log(lb/Lb);
        fa = EA*straina;
        fb = EA*strainb;
        R = P + fb*nb + fa*na;
        ka = EA/la*(na'*na) + fa/la*(I - na'*na);
        kb = EA/lb*(nb'*nb) + fb/lb*(I - nb'*nb);
        kt = ka+kb;
        u = -(kt\R')' + u;
        if norm(R) < 10*10^-13
            Ans(i) = [u(2)];
            iter(i) = it;
            break
        end
    end
end

% plot load-displacement 
p = 0.9715:.0005:0.982;
figure
plot(Ans,p,'b')
title('Pcr')
xlabel('displacement (m)')
ylabel('Py (kN)')