function [A] = RouseModelA(N,beta)
% N is the number of particles
% beta is the coefficient, notice that final MSD ~ t^(1-1/beta)
% the interaction matrix is A(n,m) in unit of k
A=zeros(N);
p=1:N-1;
for n=1:N
    for m=1:N
        s=((sin(p*pi/2/N)).^beta).*cos((n-0.5)*p*pi/N).*cos((m-0.5)*p*pi/N);
        A(n,m)=sum(s)*4/N; % k is omitted
    end
end
end