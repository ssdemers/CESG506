%% problem 1
I = eye(2);

% part 1/2
A = [1,1,2;1,4,1;1,4,4];
B = [3,4;5,1;6,3];
f = A\B;
F = f(2:3,:)';

% part 3
W = [3;-1];
w = [2'-3];
C = F'*F;
[v,d1] = eig(C);
N1 = v(:,1);
N2 = v(:,2);
U = sqrtm(C);
R = F*U^-1;
V = R*U*R';
[v,d2] = eig(V);
n1 = v(:,1);
n2 = v(:,2);
Wh = U*W;
wcheck1 = R*Wh;
wh = R*W;
wcheck2 = V*wh;

% part 4
E = 0.5*(C-I);
e = 0.5*(I-inv(F')*inv(F));

% part 5
[vE,dE] = eig(E);
[ve,de] = eig(e);
EcheckN1 = (E - dE(1,1)*I)*N1;
EcheckN2 = (E - dE(2,2)*I)*N2;
echeckn1 = (e - de(1,1)*I)*n1;
echeckn2 = (e - de(2,2)*I)*n2;

% part 6
n1 = R*N1;
n2 = R*N2;

% part 7
push = inv(F')*E*inv(F);

% part 8
push = F'*e*F