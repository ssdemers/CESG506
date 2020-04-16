% set placeholders for later plotting
P1 = zeros();
P2 = zeros();
P3 = zeros();
P4 = zeros();
Pd1 = zeros();
Pd2 = zeros();
Pd3 = zeros();
Pd4 = zeros();
it = 1;
% iterate through displacements solving Py
for u = 0:0.05:1.25
    L = sqrt(5.5^2+0.5^2);
    l = sqrt(5.5^2 + (0.5-u)^2);
    la = l/L;
    e = [la-1
        0.5*(la^2-1)
        0.5*(1-1/la^2)
        0.5*log(la^2)];
    EA = 1300;
    fk = zeros();
    for i = 1:length(e)
        f(1,i) = EA*e(i,1);
    end
    % fill in undeformed forces
    P1(it) = -f(1)*(0.5)/L;
    P2(it) = -f(2)*(0.5)/L;
    P3(it) = -f(3)*(0.5)/L;
    P4(it) = -f(4)*(0.5)/L;
    % fill in deformed forces
    Pd1(it) = -f(1)*(0.5-u)/L;
    Pd2(it) = -f(2)*(0.5-u)/L;
    Pd3(it) = -f(3)*(0.5-u)/L;
    Pd4(it) = -f(4)*(0.5-u)/L;
    it = it+1;
end

% plot undeformed
x = 0:0.05:1.25;
figure
plot(x,P1,'b',x,P2,'g',x,P3,'r',x,P4,'g')
legend('1a','1b','1c','1d')
title('Undeformed')
xlabel('displacement (m)')
ylabel('Py (kN)')

% plot deformed
figure
plot(x,Pd1,'k',x,Pd2,'c',x,Pd3,'y',x,Pd4,'m')
legend('1a','1b','1c','1d')
title('Deformed')
xlabel('displacement (m)')
ylabel('Py (kN)')