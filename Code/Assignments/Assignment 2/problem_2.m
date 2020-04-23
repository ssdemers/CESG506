% set initial conditions
EA = 2100;
nDof = 3;
I = eye(nDof);
u = [0,0,0,0,0,0];
Lavec = [5.5,0.5,1.25];
Lbvec = [5.5,0.5,3.75];
Lcvec = [5.5,0.5,-1.25];
Ldvec = [-4,0.5,1.25];
Levec = [-4,0.5,-3.75];
Lfvec = [-4,0.5,-1.25];
Lgvec = [0,0,-2.5];
La = sqrt(dot(Lavec,Lavec));
Lb = sqrt(dot(Lbvec,Lbvec));
Lc = sqrt(dot(Lcvec,Lcvec));
Ld = sqrt(dot(Ldvec,Ldvec));
Le = sqrt(dot(Levec,Levec));
Lf = sqrt(dot(Lfvec,Lfvec));
Lg = sqrt(dot(Lgvec,Lgvec));
Na = Lavec/La;
Nb = Lbvec/Lb;
Nc = Lcvec/Lc;
Nd = Ldvec/Ld;
Ne = Levec/Le;
Nf = Lfvec/Lf;
Ng = Lgvec/Lg;
P = [0,-.999,0,0,0,0];

% solution placeholders
Ans = zeros();
iter = zeros();

% initialize
ev = [0,1,0,0,0,0];
R = [0;0;0;0;0;0];
ka = EA/La*(Na'*Na);
kb = EA/Lb*(Nb'*Nb);
kc = EA/Lc*(Nc'*Nc);
kd = EA/Ld*(Nd'*Nd);
ke = EA/Le*(Ne'*Ne);
kf = EA/Lf*(Nf'*Nf);
kg = EA/Lg*(Ng'*Ng);
kt = [ka+kd+ke+kg,-kg;-kg,kb+kc+kf+kg];
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
        lavec = Lavec+u(1:3);
        lbvec = Lbvec+u(4:6);
        lcvec = Lcvec+u(4:6);
        ldvec = Ldvec+u(1:3);
        levec = Levec+u(1:3);
        lfvec = Lfvec+u(4:6);
        lgvec = Lgvec+u(1:3)-u(4:6);
        la = sqrt(dot(lavec,lavec));
        lb = sqrt(dot(lbvec,lbvec));
        lc = sqrt(dot(lcvec,lcvec));
        ld = sqrt(dot(ldvec,ldvec));
        le = sqrt(dot(levec,levec));
        lf = sqrt(dot(lfvec,lfvec));
        lg = sqrt(dot(lgvec,lgvec));
        na = lavec/la;
        nb = lbvec/lb;
        nc = lcvec/lc;
        nd = ldvec/ld;
        ne = levec/le;
        nf = lfvec/lf;
        ng = lgvec/lg;
        straina = log(la/La);
        strainb = log(lb/Lb);
        strainc = log(lc/Lc);
        straind = log(ld/Ld);
        straine = log(le/Le);
        strainf = log(lf/Lf);
        straing = log(lg/Lg);
        fa = EA*straina;
        fb = EA*strainb;
        fc = EA*strainc;
        fd = EA*straind;
        fe = EA*straine;
        ff = EA*strainf;
        fg = EA*straing;
        ka = EA/La*(na'*na) + fa/la*(I - na'*na);
        kb = EA/Lb*(nb'*nb) + fb/lb*(I - nb'*nb);
        kc = EA/Lc*(nc'*nc) + fc/lc*(I - nc'*nc);
        kd = EA/Ld*(nd'*nd) + fd/ld*(I - nd'*nd);
        ke = EA/Le*(ne'*ne) + fe/le*(I - ne'*ne);
        kf = EA/Lf*(nf'*nf) + ff/lf*(I - nf'*nf);
        kg = EA/Lg*(ng'*ng) + fg/lg*(I - ng'*ng);
        kt = [ka+kd+ke+kg,-kg;-kg,kb+kc+kf+kg];
        K = kt*u';
        R(1:3) = P(1:3)'*lambda + fa*na'+fd*nd'+fe*ne'+fg*ng';
        R(4:6) = fb*nb'+fc*nc'+ff*nf'-fg*ng';
        if norm(R+gu) < 10*10^-12
            Ans(i,1) = [P(2)*lambda];
            Ans(i,2) = u(1);  
            Ans(i,3) = u(3);
            Ans(i,4) = u(4);
            Ans(i,5) = u(5);
            Ans(i,6) = u(6);
            iter(i) = it;
            break
        else
        gu = dot(ev,u)-(.0-.01*i);
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
figure
plot(v,Ans(:,1),'b',Ans(:,2),Ans(:,1),'r',Ans(:,3),Ans(:,1),'g',Ans(:,5),Ans(:,1),'m')
title('Load-Displacement')
xlabel('displacement (m)')
ylabel('Py (kN)')
legend('Node 5 uy','Node 5 ux','Node 6 ux','Node 6 uy')

% plot node 5/6 path (z-x view)
figure
plot(Ans(:,3),v,'b',Ans(:,6)+2.5,Ans(:,5),'r')
title('Tracing node 5/6 (z-x view)')
xlabel('uz (m)')
ylabel('ux (m)')
legend('Node 5','Node 6')



