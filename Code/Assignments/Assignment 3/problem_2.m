% set initial conditions
EA = 2100;
nDof = 2;
I = eye(nDof);
u = [0;0];
Lavec = [5.5;0.5];
Lbvec = [-4.0;0.5];
La = sqrt(dot(Lavec,Lavec));
Lb = sqrt(dot(Lbvec,Lbvec));
Na = Lavec/La;
Nb = Lbvec/Lb;
P = [0;-0.999];
alpha = 0.00;

% solution placeholders
lambdaT = [];
uT = {};
Residuals = [];
arclength = [];

% initialize
R = [0;0];
g = 0;
du0 = kt\R;
u1 = kt\P;
un = [0;0];
lambda = 0;
lambdan = 0;
s = 0;

lambdaT(1) = lambdan;
uT{1} = un;

iter = 1;
while lambdan < 2
    iter = iter+1;
    if iter == 2
        lambda = lambdan + 0.2;
        u = un + u1*lambda;
        ds = sqrt(dot(u-un,u-un) + alpha*(lambda - lambdan)^2);
    else
        u = 2*uT{end} - uT{end-1};
        lambda = 2*lambdaT(end) - lambdaT(end-1);
    end

    for i = 1:10
        lavec = Lavec + u;
        lbvec = Lbvec + u;
        la = sqrt(dot(lavec,lavec));
        lb = sqrt(dot(lbvec,lbvec));
        na = lavec/la;
        nb = lbvec/lb;
        straina = log(la/La);
        strainb = log(lb/Lb);
        fa = EA*straina;
        fb = EA*strainb;
        Ft = fa*na + fb*nb;
        udiff = u - un;
        lambdadiff = lambda - lambdan;
        g = -ds^2 + dot(udiff,udiff) + alpha*lambdadiff^2;
        R = P*lambda - Ft;
        Residuals(i,iter) = norm(R+g);
        if norm(R+g) < 10*10^-10
            uT{iter} = u;
            lambdaT(iter) = lambda;
            arclength(iter) = s + ds;
            un = u;
            lambdan = lambda;
            s = s + ds;
            break
        end
        ka = EA/la*(na*na') + fa/la*(I - na*na');
        kb = EA/lb*(nb*nb') + fb/lb*(I - nb*nb');
        kt = ka+kb;
        du0 = kt\R;
        u1 = kt\P;
        dlambda = (-g - 2*dot(udiff,du0))/(2*dot(udiff,u1)+2*alpha*lambdadiff);
        du = du0 + u1*dlambda;
        lambda = lambda + dlambda;
        u = u + du;
    end
    if iter > 200   % stop loop if lambda doesnt exceed specified value 
        break
    end
end

disp = cell2mat(uT);

% run reference plot
reference_plot_p2;

% plot load-displacement 
figure
scatter(disp(1,:),lambdaT)
hold on
scatter(disp(2,:),lambdaT)
hold on
plot(v,Ans(:,1),'r',Ans(:,2),Ans(:,1),'b')
title('Load-Displacement')
xlabel('displacement (m)')
ylabel('Lambda')
legend('arc length ux','arc length uy','disp uy','disp ux')

