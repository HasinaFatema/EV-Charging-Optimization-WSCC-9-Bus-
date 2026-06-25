clc;
clear;
close all;

j = sqrt(-1);

%% ================= LINE DATA =================
% from to R X Bc tap
linedata = [
    1 4 0     0.0576 0     1
    2 7 0     0.0625 0     1
    3 9 0     0.0586 0     1
    4 5 0.01  0.085  0.088 1
    4 6 0.017 0.092  0.079 1
    5 7 0.032 0.161  0.153 1
    6 9 0.039 0.17   0.179 1
    7 8 0.0085 0.072 0.0745 1
    8 9 0.0119 0.1008 0.1045 1
];

%% ================= BUS DATA =================
%busno pg qg pl ql v del bustype qmax qmin
busdata = [
    1 0    0 0    0    1.04  0 1 1   -1
    2 1.63 0 0    0    1.025 0 3 0.06 -1
    3 0.85 0 0    0    1.025 0 3 0.5  -0.1
    4 0    0 0    0    1     0 2 0    0
    5 0    0 1.25 0.5  1     0 2 0    0
    6 0    0 0.9  0.3  1     0 2 0    0
    7 0    0 0    0    1     0 2 0    0
    8 0    0 1    0.35 1     0 2 0    0
    9 0    0 0    0    1     0 2 0    0
];

%% ======================================================
% PSO SETTINGS
%% ======================================================
nPop = 50;
MaxIter = 40;

% Decision variables:
% x(1)=Q5 , x(2)=Q6
% Q8 = 0.3 - Q5 - Q6

lb = 0;
ub = 0.3;

%% ================= INITIAL POPULATION =================
x = zeros(2,nPop);

for i = 1:nPop
    q5 = rand*0.3;
    q6 = rand*(0.3-q5);
    x(:,i) = [q5;q6];
end

v = zeros(size(x));

pbest = x;
gbest = x(:,1);

BestLoss = zeros(MaxIter,1);

%% ================= INITIAL GLOBAL BEST =================
vals = zeros(1,nPop);
for i = 1:nPop
    vals(i) = fx(x(:,i),busdata,linedata);
end
[~,idx] = min(vals);
gbest = x(:,idx);

%% ======================================================
% PSO MAIN LOOP
%% ======================================================
for it = 1:MaxIter

    w = 0.9 - (0.5*it/MaxIter);
    c1 = 2;
    c2 = 2.5;

    for i = 1:nPop

        v(:,i) = w*v(:,i) ...
               + c1*rand*(pbest(:,i)-x(:,i)) ...
               + c2*rand*(gbest-x(:,i));

        x(:,i) = x(:,i) + v(:,i);

        % limits
        x(x<lb)=lb;
        x(x>ub)=ub;

        % enforce Q5+Q6 <= 0.3
        if sum(x(:,i)) > 0.3
            x(:,i)=0.3*x(:,i)/sum(x(:,i));
        end

        % update pbest
        if fx(x(:,i),busdata,linedata) < fx(pbest(:,i),busdata,linedata)
            pbest(:,i)=x(:,i);
        end
    end

    % update gbest
    vals = zeros(1,nPop);
    for i = 1:nPop
        vals(i)=fx(x(:,i),busdata,linedata);
    end

    [BestLoss(it),idx]=min(vals);
    gbest=x(:,idx);

end

%% ================= BASE CASE =================
BaseLoss = NR_loadflow(busdata, linedata, [0 0 0]);

%% ======================================================
% FINAL RESULTS
%% ======================================================
Q5 = gbest(1);
Q6 = gbest(2);
Q8 = 0.3 - Q5 - Q6;

MinLoss = NR_loadflow(busdata,linedata,[Q5 Q6 Q8]);

disp('=========== RESULTS ===========')

fprintf('Base Case Loss (Qev=0 at Bus 5,6,8) = %.4f MW\n', BaseLoss)

disp('=========== OPTIMAL EV REACTIVE POWER ===========')
fprintf('Bus 5 Qev = %.4f pu\n',Q5)
fprintf('Bus 6 Qev = %.4f pu\n',Q6)
fprintf('Bus 8 Qev = %.4f pu\n',Q8)

fprintf('\nMinimum Total Loss = %.4f MW\n',MinLoss)

% individual display
fprintf('\nLoss contribution cases checked simultaneously at buses 5,6,8\n')

%% GRAPH
figure
plot(BestLoss,'LineWidth',2)
grid on
xlabel('Iteration')
ylabel('Loss (MW)')
title('PSO Convergence')

%% ======================================================
% OBJECTIVE FUNCTION
%% ======================================================
function loss = fx(x,busdata,linedata)

Q5 = x(1);
Q6 = x(2);
Q8 = 0.3 - Q5 - Q6;

if Q8 < 0
    loss = 1e6;
else
    loss = NR_loadflow(busdata,linedata,[Q5 Q6 Q8]);
end

end

%% ======================================================
% NR LOAD FLOW FUNCTION
%% ======================================================
function Total_Loss = NR_loadflow(busdata, linedata, Qval)

j = sqrt(-1);

% Add EV reactive loads simultaneously
busdata(5,5) = busdata(5,5) + Qval(1);
busdata(6,5) = busdata(6,5) + Qval(2);
busdata(8,5) = busdata(8,5) + Qval(3);

[ybus,nbus] = ybus1(linedata);

vm = busdata(:,6);
delta = busdata(:,7);

dpq=zeros(1,nbus);
dpv=dpq;
dns=dpq;

dpq(busdata(:,8)==2)=1;
dpv(busdata(:,8)==3)=1;
dns(busdata(:,8)>1)=1;

psh = busdata(:,2)-busdata(:,4);
qsh = busdata(:,3)-busdata(:,5);

iter=1;
maxerror=100;

while maxerror>1e-4 && iter<200

    v = vm.*cos(delta)+j*vm.*sin(delta);
    s = v.*conj(ybus*v);

    pcalc = real(s);
    qcalc = imag(s);

    delp = psh-pcalc;
    delq = qsh-qcalc;

    dc = [delp(dns==1);delq(dpq==1)];

    Jx = ybus.*repmat(v.',nbus,1).*repmat(conj(v),1,nbus);

    J11 = -imag(Jx);
    J11 = J11-diag(sum(J11,2));
    J11 = J11(dns==1,dns==1);

    J12 = real(Jx);
    J12 = (J12+diag(sum(J12,2)))./repmat(vm',nbus,1);
    J12 = J12(dns==1,dpq==1);

    J21 = -real(Jx);
    J21 = J21-diag(sum(J21,2));
    J21 = J21(dpq==1,dns==1);

    J22 = -imag(Jx);
    J22 = (J22+diag(sum(J22,2)))./repmat(vm',nbus,1);
    J22 = J22(dpq==1,dpq==1);

    J = [J11 J12;J21 J22];

    delv = J\dc;

    edit = [dns';dpq'];
    edit(edit==1)=delv;

    delta = delta + edit(1:nbus);
    vm = vm + edit(nbus+1:end);

    iter = iter + 1;
    maxerror = max(abs(delv));

end

%% LOSS
y = 1./(linedata(:,3)+1j*linedata(:,4));
Vn = vm.*(cos(delta)+1j*sin(delta));

Iline = (Vn(linedata(:,1))-Vn(linedata(:,2))).*y;
Sl = (abs(Iline).^2)./y;

Total_Loss = sum(real(Sl))*100;

end

%% ======================================================
% YBUS
%% ======================================================
function [yb,nbus] = ybus1(linedata)

j = sqrt(-1);

nbus=max(max(linedata(:,1)),max(linedata(:,2)));
nbr=size(linedata,1);

z=linedata(:,3)+j*linedata(:,4);
y=1./z;

yb=zeros(nbus,nbus);

for i=1:nbr

    from=linedata(i,1);
    to=linedata(i,2);
    tap=linedata(i,6);
    b=linedata(i,5);

    yb(from,from)=yb(from,from)+y(i)/(tap^2)+j*b;
    yb(to,to)=yb(to,to)+y(i)/(tap^2)+j*b;

    yb(from,to)=yb(from,to)-y(i)/tap;
    yb(to,from)=yb(from,to);

end

end