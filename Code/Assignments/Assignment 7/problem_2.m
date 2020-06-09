%% problem 2

% part 1 on paper

% part 2
J = [2,3;-3,-1];
J0 = [3,3;-1,2];

% part 3
F = J*inv(J0);

% part 4
Wh = [1;0];
W = J*Wh;
w = J0*Wh;

% part 5
Eh = 0.5*(J'*J-J0'*J0);

% part 6
push = inv(J')*Eh*inv(J);

% part 6
pull = inv(J0')*Eh*inv(J0);