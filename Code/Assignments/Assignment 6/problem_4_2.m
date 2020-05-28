% set initial conditions
EA = 10000;
EI = 10;
w0 = 1;
Nnode = 5;  % plug in 1 more than desired number of elements in mesh
nDOF = 3;
alpha = 0.01;
idxFixed = [1:2,Nnode*3-1];
nNds = Nnode;
idxFree = [3:1:Nnode*3-3,Nnode*3-2,Nnode*3]; 
Lel = (5/(Nnode-1));
U = zeros(3,Nnode);

% %     nodeID     X       Y
nodes=[  ];
x = 0;
for i = 1:Nnode;
   nodes(i,1) = i;
   nodes(i,2) = x;
   nodes(i,3) = -0.008*x^2+.04*x;
   x = x + Lel;
end
  
% %      Member       i       j        EA       EI
mesh=[  ];
for i = 1:(Nnode-1);
   mesh(i,1) = i;
   mesh(i,2) = i;
   mesh(i,3) = i+1;
   mesh(i,4) = EA;
   mesh(i,5) = EI;
end
       
% %      loads
P = zeros(3,Nnode);
P(1,Nnode) = -3.96;

[Pcr, Fsys, Ksys] = BeamAssemble(nodes, mesh, support, P, U, idxFixed, idxFree, nDOF, nNds);       

% solution placeholders
lambdaT = [];
uT = {};
midx = [];
midy = [];
q1x = [];
q1y = [];
q3x = [];
q3y = [];
dk = [];
Residuals = [];

% initialize
R = zeros(Nnode*3-3,1);
g = 0;
du0 = Ksys\R;
u1 = Ksys\Pcr;
un = zeros(1,Nnode*3-3)';
lambda = 0;
lambdan = 0;
s = 0;

lambdaT(1) = lambdan;
uT{1} = un;

iter = 1;
for i = 1:50                      
    iter = iter+1;
    if iter == 2
        lambda = lambdan + 0.5;
        u = un + u1*lambda;
        ds = sqrt(dot(u-un,u-un) + alpha*(lambda - lambdan)^2);
    else
        u = 2*uT{end} - uT{end-1};
        lambda = 2*lambdaT(end) - lambdaT(end-1);
    end

    for i = 1:10
        U(1,2:Nnode) = u(2:3:end);
        U(2,2:Nnode-1) = u(3:3:end-1);
        U(3,1:Nnode) = [u(1:3:end);u(end)];
        [Pcr, Fsys, Ksys] = BeamAssemble(nodes, mesh, support, P, U, idxFixed, idxFree, nDOF, nNds);
        udiff = u - un;
        lambdadiff = lambda - lambdan;
        g = -(ds^2) + dot(udiff,udiff) + alpha*lambdadiff^2;
        R = Pcr*lambda - Fsys;
        Residuals(i,iter) = norm(R+g);
        if norm(R+g) < 10*10^-10
            uT{iter} = u;
            lambdaT(iter) = lambda;
            un = u;
            lambdan = lambda;
            s = s + ds;
            midx(iter) = U(1,Nnode/2+.5);
            midy(iter) = U(2,Nnode/2+.5);
            q1x(iter) = U(1,Nnode/4+.75);
            q1y(iter) = U(2,Nnode/4+.75);
            q3x(iter) = U(1,Nnode*.75+.25);
            q3y(iter) = U(2,Nnode*.75+.25);
            dk(iter) = min(eig(Ksys));
            break
        end
        du0 = Ksys\R;
        u1 = Ksys\Pcr;
        dlambda = (-g - 2*dot(udiff,du0))/(2*dot(udiff,u1)+2*alpha*lambdadiff);
        du = du0 + u1*dlambda;
        lambda = lambda + dlambda;
        u = u + du;
    end
      
end

% plot final disp
figure
plot(nodes(:,2)'+U(1,:),nodes(:,3)'+U(2,:),'b',nodes(:,2),nodes(:,3),'r')
title('Beam Displacement')
xlabel('position (m)')
ylabel('position (m)')
legend('Deformed','Undeformed')
% plot disp vs load factor
figure
plot(midx,lambdaT,'b',midy,lambdaT,'r',q1x,lambdaT,'-.b',q1y,lambdaT,'-.r',q3x,lambdaT,'-db',q3y,lambdaT,'-dr','MarkerSize',5)
title('Beam Displacement')
xlabel('position (m)')
ylabel('load factor')
legend('1/2 u','1/2 v','1/4 u','1/4 v','3/4 u','3/4 v')
% plot disp vs load factor
figure
plot(lambdaT(2:end),dk(2:end),'b')
title('Det(K)')
xlabel('lambda')
ylabel('det(K)')
P_cr = -(dk(end-1))*((lambdaT(end)-lambdaT(end-1))/(dk(end)-dk(end-1)))+lambdaT(end-1)
