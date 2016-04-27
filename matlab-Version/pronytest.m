clc;
close all;
clear;

%% init setting
fs = 40;
NX = 50;
t = 0:1/fs:NX/fs;
x = 2*exp(-3*t).*cos(8*pi*t) + exp(-4*t).*cos(6*pi*t);

figure;
% subplot(211)
% plot(x,'g')
hold on;

p = 4; % ??

c = zeros(p,p);
for i = 1:p
    for j = 1:p
        c(i,j) = x(p-i+1:NX-i)*x(p-j+1:NX-j)' + x(i+1:NX-p+i)*x(j+1:NX-p+j)';
    end
end

b = zeros(p,1);
for i = 1:p
    b(i,1) = x(p-i+1:NX-i)*x(p+1:NX)' + x(i+1:NX-p+i)*x(1:NX-p)';
end

a = [1;c\(-b)]';

z = roots(a);

zmat = zeros(p,p);
for i = 1:p
    for j = 1:p
        zmat(i,j) = real(z(j)^i);
    end
end

% z = sort(z,'descend');

A = zmat\x(1:p)';

% A = sort(A,'descend');

rx = zeros(1,NX);
for i = 1:p
    rx(i) = x(i);
end
%%
for i = p+1:NX
    rx(i) = -sum(x(i-p:i-1).*fliplr(a(2:p+1)));
end

% subplot(212)
% plot(rx,'r:');

%% calc ValidLength
maxer = abs(A(1));
for i = 1:p
    if sum(abs(A(i:p))) < maxer*0.005
        ValidLength = i;
        break;
    end
    ValidLength = i;
end
ValidLength = p;

%% re-construction by prony
rxPRONY = zeros(1, NX);
error = zeros(1,NX);
for i = 1:NX
    rxTmp = 0;
    for j = 1:ValidLength
        rxTmp = rxTmp + A(j)*z(j)^i;
        
    end
   
    error(i) = rxTmp - x(i);
    rxPRONY(i) = real(rxTmp);
end

% subplot(211)
plot(x);
hold on;
% subplot(212)
plot(rxPRONY,'r');hold on;
plot(rxPRONY - x(1:NX), ':');





