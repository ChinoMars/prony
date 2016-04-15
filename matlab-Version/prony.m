function [pronyParam, ValidLength] = prony(x, p, fs)
%
%
%
%

%% get C matrix by Marple Algorithm
NX = length(x);

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

%% get z expression and zmat as a matrix
z = roots(a);

zmat = zeros(p,p);
for i = 1:p
    for j = 1:p
        zmat(i,j) = z(j)^i;
    end
end
z = sort(z, 'descend');

%% get A, the Amplitude of prony expression
A = zmat\x(1:p)';
A = sort(A, 'descend');

%% re-construction by AR module
rxAR = zeros(1, NX);
for i = 1:p
    rxAR(i) = x(i);
end
for i = p+1:NX
    rxAR(i) = -sum(x(i-p:i-1).*fliplr(a(2:p+1)));
end

%% get Valid Length of prony parameters, too small parameters will be abandoned
maxer = abs(A(1));
for i = 1:p
    if sum(abs(A(i:p))) < maxer*0.005
        ValidLength = i;
        break;
    end
    ValidLength = i;
end

%% re-construction by prony
ztmp = zeros(1, NX);
rxPRONY = zeros(1, NX);
for i = 1:ValidLength
    for j = 1:NX
        ztmp(j) = z(j)^j;
    end
    rxPRONY = rxPRONY + ztmp;
end


%% get alpha_i and omega_i
ai = zeros(1,p);
oi = zeros(1,p);
for i = 1:p
    ai(i) = real(log(z(i)))*fs;
    oi(i) = image(log(z(i)))*fs;
end

%% return a struct
pronyParam.a = a;
pronyParam.z = z;
pronyParam.A = A;
pronyParam.alpha = ai;
pronyParam.omega = oi;


