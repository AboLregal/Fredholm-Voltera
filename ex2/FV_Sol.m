clc;clear;close all

% global N H1 H2 Y const1 const2 lamda1 lamda2
global N
for N=2:2:10
    
%     lamda1=1;
%     lamda2=1;
%     Y=inline('((s^6)/(-30))+((s^4)/3)-(s^2)+(5*s/3)-(5/4)');
%     const1=inline('(si)-(0.5*(si)*(sk+1))');
%     const2=inline('si+((sk+1)/2)');
%     H1=inline('(C*T)^2');
%     H2=inline('(C*T)');
    
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
        est(i)=cx*to;
    end
    true=cos(t(:));
    est=est(:);
    err=abs(est-true);
    fprintf('\n ********** \n N=%d',N)
    table(true,est,err)
    figure()
    plot(t,true,'b');hold on
    plot(t(2:end),est(2:end),'*r')
    legend('Exact','Estimated')
    title(['N=',num2str(N)])
end