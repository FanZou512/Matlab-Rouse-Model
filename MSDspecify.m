function [ MSD ]=MSDspecify(t,varargin)
% This function calculate the MSD of given data
% input should be x,y,z, ... coordinates in matrix form, with same size
% each input is one dimension
% ------------ How to use ------------
% [MSD] = MSDspecify(time, xt,yt,... );
[m,n]=size(varargin{1});
MSD=zeros(length(t),1);
for q=1:length(t)
    k=t(q);
    xk=zeros(m-k,n);
    for p=1:nargin-1
        xk=xk+(varargin{p}(k+1:end,:)-varargin{p}(1:end-k,:)).^2;
    end
    %threshold=0.3*k/0.0036;
    %xk(xk>threshold)=NaN;
    MSD(q)=mean(xk(:));
    if mod(q,20)==0
        disp([num2str(q),' intervals calculated!'])
    end
end
end