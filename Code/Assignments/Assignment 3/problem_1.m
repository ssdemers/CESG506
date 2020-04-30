% set initial conditions
EA = 2100;
nDof = 1;
I = eye(nDof);
u = 0;
Lavec = [5.5;0.5];
La = sqrt(dot(Lavec,Lavec));
Na = Lavec/La;
P = -0.3;
alpha = 0;

% solution placeholders
lambdaT = [];
uT = [];
Residuals = [];
arclength = [];
Res = [];

% initialize
R = 0;
kt = EA/La*(Na(2)*Na(2));
g = 0;
du0 = kt\R;
u1 = kt\P;
un = 0;
lambda = 0;
lambdan = 0;
s = 0;

lambdaT(1) = lambdan;
uT(1) = un;

iter = 1;
count = 0;
while lambdan < 2
    iter = iter+1;
    if iter == 2
        lambda = lambdan + 0.25;
        u = un + u1*lambda;
        ds = sqrt(dot(u-un,u-un) + alpha*(lambda - lambdan)^2);
    else
        u = 2*uT(end) - uT(end-1);
        lambda = 2*lambdaT(end) - lambdaT(end-1);
    end

    for i = 1:10
        count = count + 1;
        lavec = Lavec + [0;u];
        la = sqrt(dot(lavec,lavec));
        na = lavec/la;
        straina = log(la/La);
        fa = EA*straina;
        Ft = fa*na(2);
        udiff = u - un;
        lambdadiff = lambda - lambdan;
        g = -ds^2 + dot(udiff,udiff) + alpha*lambdadiff^2;
        R = P*lambda - Ft;
        Residuals(count) = norm(R+g);
        if norm(R+g) < 10*10^-10
            uT(iter) = u;
            lambdaT(iter) = lambda;
            arclength(iter) = s + ds;
            Res(iter) = norm(R+g);
            un = u;
            lambdan = lambda;
            s = s + ds;
            break
        end
        kt = EA/la*(na(2)*na(2)) + fa/la*(I - na(2)*na(2))
        du0 = kt\R;
        u1 = kt\P;
        dlambda = (-g - 2*dot(udiff,du0))/(2*dot(udiff,u1)+2*alpha*lambdadiff);
        du = du0 + u1*dlambda;
        lambda = lambda + dlambda;
        u = u + du;
    end
    if iter > 200  
        break
    end
end

% run script for reference plot
reference_plot_p1;
% plot load-displacement 
figure
scatter(uT,lambdaT)
hold on
plot(-x,Pd4,'r')
title('Load-Displacement')
xlabel('displacement (m)')
ylabel('Lambda')
legend('arc length','Hw 1 solution')

% plot error and step count
figure
ct = linspace(1,count,count);
plot(ct,Residuals,'r')
set(gca, 'YScale', 'log')
title('Norm of R')
xlabel('Step count')
ylabel('R')

% plot error and arc length
figure
plot(arclength,Res,'r')
set(gca, 'YScale', 'log')
title('Norm of R')
xlabel('Arc length')
ylabel('R')