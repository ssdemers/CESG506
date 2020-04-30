% set initial conditions
EA = 2100;
nDof = 2;
I = eye(nDof);
u = [0,0];
Lavec = [5.5,0.5];
Lbvec = [-4,0.5];
La = sqrt(dot(Lavec,Lavec));
Lb = sqrt(dot(Lbvec,Lbvec));
Na = Lavec/La;
Nb = Lbvec/Lb;
P = [0,-0.999];

% solution placeholders
Ans = zeros();
iter = zeros();

% initialize
ev = [0,1];
R = [0;0];
ka = EA/La*(Na'*Na); 
kb = EA/Lb*(Nb'*Nb);
kt = ka+kb;
gu = dot(ev,u)+0.01;
du0 = -(kt\R);
u1 = -(kt\P');
dlambda = -(gu+dot(ev,du0))/dot(ev,u1');
lambda = dlambda;
u = du0'+(u1*dlambda)';


% solve
for i = 1:120
    u(2) = 0-.01*i;
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
        ka = EA/la*(na'*na) + fa/la*(I - na'*na);
        kb = EA/lb*(nb'*nb) + fb/lb*(I - nb'*nb);
        kt = ka+kb;
        R = P'*lambda + fb*nb' + fa*na';
        if norm(R+gu) < 10*10^-10
            Ans(i,1) = [P(2)*lambda];
            Ans(i,2) = u(1);
            iter(i) = it;
            break
        else
        gu = dot(ev,u)-(0-.01*i);
        du0 = -(kt\R);
        u1 = -(kt\P');
        dlambda = -(gu+dot(ev,du0))/dot(ev,u1');
        lambda = lambda + dlambda;
        u = u+ du0'+(u1*dlambda)';
        end
    end
end

% plot load-displacement 
v = -0.01:-.01:-1.2;
% figure
% plot(v,Ans(:,1),'b',Ans(:,2),Ans(:,1),'r')
% title('Load-Displacement')
% xlabel('displacement (m)')
% ylabel('Py (kN)')
% legend('uy','ux')
% 
% % plot displacement of node 2
% figure
% plot(Ans(:,2),v,'b')
% title('Tracing node 2')
% xlabel('ux (m)')
% ylabel('uy (m)')

