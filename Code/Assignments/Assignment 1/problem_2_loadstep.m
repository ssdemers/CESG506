% set initial conditions
EA = 2100;
nDof = 2;
I = eye(nDof);
u = [0,0];
Lavec = [5.5,0.5];
Lbvec = [-4,0.5];
La = sqrt(dot(Lavec,Lavec));
Lb = sqrt(dot(Lbvec,Lbvec));
Pcr = [0,0.9815];
lambda = [0.25, 0.5, 0.75, 0.99, 0.999];
% placeholders
Ans = zeros();
iter = zeros();

% find disp
for i = 1:5
    P = Pcr*lambda(i)
    for it = 1:10
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
        Ans(it,i) = norm(R);
        if norm(R) < 10*10^-13
            iter(i) = it;
            break
        end
    end
end


