function [ MSD, MSDerr, Num ]=RMSDspecify(t,varargin)
% This function calculate the relative MSD of given data
% input should be x,y,z, ... coordinates in matrix form, with same size
% each input is one dimension
% each input data contain 2n columns, which are the coordinates of two dots
% ------------ How to use ------------
% [MSD] = RMSDspecify(time, xt,yt,... );
[m,n]=size(varargin{1});
if isempty(t)
    t=1:m-1;
end
n=2*floor(n/2);
MSD=zeros(length(t),1);
MSDerr=MSD;
Num=MSD;
for q=1:length(t)
    k=t(q);
    xk=zeros(m-k,n/2);
    for p=1:nargin-1
        xk=xk+(varargin{p}(k+1:end,1:2:n-1)-varargin{p}(k+1:end,2:2:n)-varargin{p}(1:end-k,1:2:n-1)+varargin{p}(1:end-k,2:2:n)).^2;
    end
    MSD(q)=nanmean(xk(:));
    MSDerr(q)=nanstd(xk(:));
    Num(q)=sum(~isnan(xk(:)));
    if mod(q,30)==0
        disp([num2str(q),' intervals calculated!'])
    end
end
end