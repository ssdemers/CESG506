function [ F , K ] = CurvedBeamElement( Ui , Uj , dispi , dispj , EA, EI )

% initialize
F = zeros(6,1);
K = zeros(6,6);

% elevation
Lel = Uj(1) - Ui(1);
hi = Ui(2);
hj = Uj(2);
dh = (hj - hi)/Lel;

qu = [dispi(1); dispj(1)];
qv = [dispi(2);dispi(3);dispj(2);dispj(3)];

% Gauss points for int.
gp = 2;
xk = [0.2113, 0.7887] ;
wk = [0.5, 0.5] ;

% solve for K and F
for it = 1:gp
    x = xk(it);
    w = wk(it)*Lel;
    dN_u = [-1,1]/Lel ;   % linear
    dN_v  = [-6*x+6*x^2, Lel*(1-4*x+3*x^2), 6*x-6*x^2, Lel*(-2*x+3*x^2)]/Lel; % cubic
    ddN_v = [-6+12*x, Lel*(-4+6*x), 6-12*x, Lel*(-2+6*x)]/(Lel^2); % cubic
                         
    du0 = dN_u*qu ;
    dv0 = dN_v*qv ; 
    phi = ddN_v*qv ; 
    ep0 = du0 + dh*dv0 + 0.5*dv0*dv0;
    f = EA*ep0;
    M = EI*phi;
    
    dNv = [0,dN_v(1),dN_v(2),0,dN_v(3),dN_v(4)];
    Bphi = [0,ddN_v(1),ddN_v(2),0,ddN_v(3),ddN_v(4)];
    Bep = [dN_u(1),(dh+dv0)*dN_v(1),(dh+dv0)*dN_v(2),dN_u(2),(dh+dv0)*dN_v(3),(dh+dv0)*dN_v(4)];
    
    K = K + (Bep'*EA*Bep+Bphi'*EI*Bphi+f*dNv'*dNv)*w;
    F = F + (Bep'*f+Bphi'*M)*w;
    
end 
end