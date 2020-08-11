clc;clear;close all

global N H1 H2 Y const1 const2 lamda1 lamda2

for N=2:2:10
    
    lamda1=0;
    lamda2=-1;
    Y=inline('s*exp(1)+1');
    const1=inline('0*(si+sk)');
    const2=inline('si+((sk+1)/2)');
    H1=inline('(C*T).^2');
    H2=inline('exp(C*T)');
    
    c=ones(1,N+1);
    cx=fsolve(@FV_1D,c);
    
    t=0:0.1:1;
    
    
    for i=1:length(t)
        z=t(i);
        for k=1:N+1
            if k==1
                to(k,1)=1;
            elseif k==2
                to(k,1)=z;
            else
                to(k,1)=(2*z*to(k-1))-to(k-2);
            end
        end
        est(i,1)=cx*to;
    end
    true=t(:);
    err=abs(true-est);
    table(true,est,err)
    figure()
    plot(t,true,'b');hold on
    plot(t(2:end),est(2:end),'*r')
    legend('Exact','Estimated')
    title(['N=',num2str(N)])
end