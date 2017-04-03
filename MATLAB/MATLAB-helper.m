clc
clear all
close all
N = 4; % number of nodes
%%%mg hpl ep RA master
% workload information 
C1 = [0.0006    0.0001    0.0006    0.0013];
C2 = [-0.2389    -0.0783    -0.2846   -0.4834];
C3 = [20.2984   8.5808  27.3119  39.6226];
vl = 130*ones(N,1); %% lower limit
vh = 200*ones(N,1); %% higher limit



A=ones(1,N); %% constraint zero: total power budger
A2 =  repmat([1,0], 1,2); %% Circut breaker 0 
A3 =  repmat([0,1], 1,2); %% Circuit breaker 1
H=2*diag(C1);
D = 650; %% total budget
D2 = 300; %% circuit breaker 0 limit
D3 = 300; %% circuit breaker 1 limit
prs = [1 0 0 1]; %% priorities
rng(1);
Nh = [0 1 0 1;1 0 1 0; 0 1 0 1; 1 0 1 0]; %% neighborhood matrix  
mu1=.05; %% mu
T=6000; %% maximum number of iteration to simulate
epsilon0=2; %% epsilon
conv_thr = .01; %% conversion threshhol

iter_mat = [];
error_mat = [];
prs_mat = [];

v_0 = vl; %% starting point.
v=[]; %% has the power cap for each iteration in its columns. 
v(:,1)=v_0; %% starting from v_0
%% cost function 1/2*x'*H*x+C2*x+sum(C3)
b=D;
e_0=(A*v_0-b)/N; %% initial condition
e_02=(A2*v_0 - D2)/(N/2);
e_03=(A3*v_0 - D3)/(N/2);
et=[];
et2 = [];
et3 = []; 
et(:,1)=e_0*ones(N,1);
et2(:,1)=e_02 .* A2; 
et3(:,1)=e_03 .* A3; 
[OptSol,OptCost]=quadprog(H.*(2.^repmat(prs,N,1)),C2.*(2.^prs),[A;A2;A3],[b;D2;D3],[],[],vl,vh); %% optimal solution
vhat = zeros(N, T); %% changes during each iteration
OptCost=(1/2*OptSol'*H*OptSol+C2*OptSol+sum(C3));
flag = 1; 
for t=1:T
    cost(:,t)=1/2*v(1:N,t)'*H*v(1:N,t)+C2*v(1:N,t)+sum(C3(1:N));
    eij=zeros(N,N);
    eij2=zeros(N,N);
    eij3=zeros(N,N);
    for i=1:N
        %epsilonT(i,t)=epsilon0;
        vhat(i,t)=-epsilon0*((2*C1(i)*v(i,t)+C2(i)).*(2.^prs(i))+(2*mu1*max(0,et(i,t))*A(i)+2*mu1*max(0,et2(i,t))*A2(i)+ 2*mu1*max(0,et3(i,t))*A3(i)));
        for j=1:N
            if Nh(i,j) %% if they are neighbors
                eij(i,j)=epsilon0*2*mu1*(et(i,t)-et(j,t));
                eij2(i,j)=epsilon0*2*mu1*(et2(i,t)-et2(j,t));
                eij3(i,j)=epsilon0*2*mu1*(et3(i,t)-et3(j,t));
            end
        end
    end
    v(:,t+1)=v(:,t)+vhat(:,t);
    v(:,t+1)=max(min(v(:,t+1),vh(1:N)),vl(1:N));
    vhat(:,t)=v(:,t+1)-v(:,t);
    et(:,t+1)=et(:,t)+sum(eij)'-sum(eij')'+A'.*vhat(:,t);
    et2(:,t+1)=et2(:,t)+sum(eij2)'-sum(eij2')'+A2'.*vhat(:,t);
    et3(:,t+1)=et3(:,t)+sum(eij3)'-sum(eij3')'+A3'.*vhat(:,t);
    %%CCgrid convergence condition
    %if abs(cost(end)-OptCost)<0.01*abs(OptCost) 
    num_neigh = 0;
    num_conv = 0; 
    for i =1:N 
        for j= 1:N 
            if (Nh(i,j) ==1) 
                num_neigh = num_neigh+1;
                if (abs(et(i,t+1)- et(j,t+1)) < conv_thr && abs(et2(i,t+1)- et2(j,t+1)) < conv_thr && abs(et3(i,t+1)- et3(j,t+1)) < conv_thr && vhat(i,t) <.05)
                    num_conv = num_conv +1;
                end
            end
        end
    end
    if (num_neigh == num_conv)
        [OptSol,OptCost]=quadprog(H.*(2.^repmat(prs,N,1)),C2.*(2.^prs),[A;A2;A3],[b;D2;D3],[],[],vl,vh);
        error = (OptSol - v(1:N,t));
        break %% converged!
    end
end
% 
figure;
h=plot(1:t,v(:,1:t),1:t,sum(v(:,1:t))-12,'m-.');
axis([0,t,100,250])
%legend('server 1','server 2', 'server 3', 'server 4', 'Total power - capacity');
ylabel('power, v','FontSize',14);
xlabel('Iterations','FontSize',14);
set(gca,'FontSize',14);
set(h,'linewidth',2);
figure;
subplot(1,3,1);
h=plot(1:t,et(:,1:t));
ylabel('error, et','FontSize',14);
xlabel('Iterations','FontSize',14);
set(gca,'FontSize',14);
set(h,'linewidth',2);
subplot(1,3,2);
h=plot(1:t,et2(:,1:t));
ylabel('error, et2','FontSize',14);
xlabel('Iterations','FontSize',14);
set(gca,'FontSize',14);
set(h,'linewidth',2);
subplot(1,3,3);
h=plot(1:t,et3(:,1:t));
ylabel('error, et3','FontSize',14);
xlabel('Iterations','FontSize',14);
set(gca,'FontSize',14);
set(h,'linewidth',2);
FINAL_POWER = v(:,end)'
SUM_POWER = sum(v(:,end))



