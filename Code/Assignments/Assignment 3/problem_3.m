% set initial conditions
EA1 = 5000;
EA2 = 2000;
U = zeros(2,24);
u = zeros(1,44);
P = zeros(2,24);
P(2,24) = -3.1;
P(2,23) = -3.1;
alpha = .1;
h = 5/11;
idxFixed = [1;2;3;4];
idxFree = [5:1:48]; 

% %     nodeID     X       Y
nodes=[  1         0.0     0.0;
         2         0.25    0.0;
         3         0.0     h;
         4         0.25    h;
         5         0.0     2*h;
         6         0.25    2*h;
         7         0.0     3*h;
         8         0.25    3*h;
         9         0.0     4*h;
         10        0.25    4*h;
         11        0.0     5*h;
         12        0.25    5*h;
         13        0.0     6*h;
         14        0.25    6*h;
         15        0.0     7*h;
         16        0.25    7*h;
         17        0.0     8*h;
         18        0.25    8*h;
         19        0.0     9*h;
         20        0.25    9*h;
         21        0.0     10*h;
         22        0.25    10*h;
         23        0.0     11*h;
         24        0.25    11*h];
       
  
% %      Member       i       j        EA
mesh=[  1            1        3        EA2;
        2            1        4        EA1;
        3            2        4        EA2;
        4            3        5        EA2;
        5            3        6        EA1;
        6            3        4        EA1;
        7            4        6        EA2;
        8            5        7        EA2;
        9            5        8        EA1;
        10           5        6        EA1;
        11           6        8        EA2;
        12           7        9        EA2;
        13           7        10       EA1;
        14           7        8        EA1;
        15           8        10       EA2;
        16           9        11       EA2;
        17           9        12       EA1;
        18           9        10       EA1;
        19           10       12       EA2;
        20           11       13       EA2;
        21           11       14       EA1;
        22           11       12       EA1;
        23           12       14       EA2;
        24           13       15       EA2;
        25           13       16       EA1;
        26           13       14       EA1;
        27           14       16       EA2;
        28           15       17       EA2;
        29           15       18       EA1;
        30           15       16       EA1;
        31           16       18       EA2;
        32           17       19       EA2;
        33           17       20       EA1;
        34           17       18       EA1;
        35           18       20       EA2;
        36           19       21       EA2;
        37           19       22       EA1;
        38           19       20       EA1;
        39           20       22       EA2;
        40           21       23       EA2;
        41           21       24       EA1;
        42           21       22       EA1;
        43           22       24       EA2;
        44           23       24       EA1];

%          Node    X       Y
support=[  1       1       1;
           2       1       1];
       
[Pcr, Fsys, Ksys] = Assemble(nodes, mesh, support, P, U, idxFixed, idxFree);       
     
% solution placeholders
lambdaT = [];
uT = {};
dispx = [];
dispy = [];
trackx = [];
tracky = [];
trackx1 = [];
tracky1 = [];
Residuals = [];

% initialize
R = zeros(44,1);
g = 0;
du0 = Ksys\R;
u1 = Ksys\Pcr;
un = zeros(44,1);
lambda = 0;
lambdan = 0;
s = 0;

lambdaT(1) = lambdan;
uT{1} = un;

iter = 1;
while u(44) > -5
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
        U(1,3:end) = u(1:2:43,1);
        U(2,3:end) = u(2:2:44,1);
        [Pcr, Fsys, Ksys] = Assemble(nodes, mesh, support, P, U, idxFixed, idxFree);
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
            dispx(iter) = U(1,24);
            dispy(iter) = U(2,24);
            dispx1(iter) = U(1,23);
            dispy1(iter) = U(2,23);
            trackx(iter-1) = U(1,24) + nodes(24,2);
            tracky(iter-1) = U(2,24) + nodes(24,3);
            trackx1(iter-1) = U(1,24) + nodes(23,2);
            tracky1(iter-1) = U(2,24) + nodes(23,3);
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

%disp = cell2mat(uT);

% plot node tracking 
figure
plot(trackx,tracky,'b',trackx1,tracky1,'r')
title('Top node tracking')
xlabel('displacement (m)')
ylabel('displacement (m)')
legend('top right node','top left node')
% plot load displacement
figure
plot(dispx,lambdaT,'b',dispy,lambdaT,'r',dispx1,lambdaT,'g',dispy1,lambdaT,'y')
title('Load-Displacement')
xlabel('displacement (m)')
ylabel('Lambda')
legend('top right: ux','top right: uy','top left: ux','top left: uy')

