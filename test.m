clc
clear
Uacc = [zeros(1,240)];
Uomega = [zeros(1,240)];
%state vector: x,y 2D position, velocity, heading
%simulate no noise system
[Xnonoise,n,t,dt,nt]=simulate_2Dcar(Uacc,Uomega);
nk=nt;
xtrue=[4;9]
[xe]=Elocation(xtrue);
%% Simulation Function calls
function [Xnonoise,n,t,dt,nt]=simulate_2Dcar(Uacc,Uomega);
dt=0.1;
nt=length(Uacc);
t=[0:dt:dt*(nt-1)];
n=4;
Xnonoise=zeros(n,nt);
Xnonoise(:,1)=[1;1;norm([7-1,5-1])/24;0.588];
for k=1:(nt-1),
    Vk=Xnonoise(3,k);
    Tk=Xnonoise(4,k);    
    Xnonoise(:,k+1) = Xnonoise(:,k) +...
        dt*[Vk*cos(Tk);Vk*sin(Tk);Uacc(k);Uomega(k)];   
end
end

%% 
function [xe]=Elocation(xtrue)
rng(100)
bA = [2;4];
bB = [6;6];
bC = [4;2];
%True location of the robot
%xtrue = [Xnonoise(1,:);Xnonoise(2,:)];
%Perfect measurements
RA = sqrt((bA(1) - xtrue(1)).^2 + (bA(2) - xtrue(2)).^2);
RB = sqrt((bB(1) - xtrue(1)).^2 + (bB(2) - xtrue(2)).^2);
RC = sqrt((bC(1) - xtrue(1)).^2 + (bC(2) - xtrue(2)).^2);
htrue=[RA RB RC];
na=100;
Rpart = diag([10 10 10]);
va = sqrtm(Rpart)*randn(3,na);
za = [RA;RB;RC]+ va;
xhat=[0;0;0];
NLS_pass=1;J=[];Jold=1;iter=0;
while NLS_pass,
    iter=iter+1;
    xj=xhat(:,iter);
    RAhat=sqrt((bA(1) - xj(1))^2 + (bA(2) - xj(2))^2);
    RBhat=sqrt((bB(1) - xj(1))^2 + (bB(2) - xj(2))^2);
    RChat=sqrt((bC(1) - xj(1))^2 + (bC(2) - xj(2))^2);
    H=[(-bA(1)+xj(1))/RA  (-bA(2)+xj(2))/RA;
        (-bB(1)+xj(1))/RB (-bB(2)+xj(2))/RB;
        (-bC(1)+xj(1))/RC (-bC(2)+xj(2))/RC];
    M1=0;M2=0;J(iter)=0;
    for k=1:na,
        M1=M1+H'*inv(Rpart)*H;
        ek=[za(:,k)-[RAhat;RBhat;RChat]];
        M2=M2+H'*inv(Rpart)*ek;
        J(iter)=J(iter)+0.5*ek'*inv(Rpart)*ek;
    end
    Sigma_x = inv(M1);
    xhat(1:2,iter+1)=xj(1:2,:)+inv(M1)*M2;
    if abs(Jold-J(iter))<1E-3,
        NLS_pass=0;
    else,
        Jold=J(iter);
    end
end
xe=xj(1:2);
end