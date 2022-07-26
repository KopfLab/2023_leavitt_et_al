%% Prediction of hydrogen alphas for archaeal lipid biosynthesis
% N. maritimus autotrophy only
% Non-linear simulated annealing 
% Leavitt, Kopf, Weber, Pearson - July, 2022
clc; clear all

% "Scenario 2", Fd is the reductant for GGR: Calculate BP stoichiometry; nothing in lines 7-29 should be modified!
% First, decide if AcCo-A enol water exchange is allowed?
Exch =      0;              %set from 0 to 1 to reflect tautomerization of AcCoA with water; also set alpha_x = 1 if Exch is 0
% Now calculate stoichiometry of H sources for biphytanes
n_TOT =     [80; 78; 76; 74];   %total H inventories
nRings =    [0; 1; 2; 3];       %BPs 0-3
A_w =       2/3;                %fraction of AcCoA from water
A_NADPH =   1-A_w;              %fraction of AcCoA from NAD(P)H
n_NADPH =   16 + A_NADPH*(47+1/3); %number of H from NAD(P)H in DGGGP
n_w =       (2/3) + (47+1/3)*A_w;  %number of H total from water
n_a =       (2/3) + A_w*(47+1/3)*(1-Exch); %number of H from unspecified water, including NAD/FAD water
n_x =       n_w - n_a;          %number of H specifically from enol tautomerization; may have same or smaller KIE as n_a?
% Next: Adding a ring means Hs are not added by Fd, i.e., for every ring, two fewer H: one less from the Fd hydride pool, 
% 2/3 less from the general water pool, and 1/3 less from the NAD(P)H pool; don't touch the enol water pool.
n_a_GGR =   8;                  %initial max number of H from GGR-water
n_NADPH_GGR = 0;                %initial max number of H from GGR-NAD(P)H
n_Fd_GGR  = 8;                  %initial max number of H from GGR-Fd (Fd-hydride)
f_W =       (n_a + n_a_GGR - nRings.*(2/3))./n_TOT; F_W=[f_W; f_W; f_W]; %vector, 12 entries
f_x =       n_x./n_TOT; F_x=[f_x; f_x; f_x]; %vector, 12 entries
f_Fd =      (n_Fd_GGR - nRings.*(1))./n_TOT; F_Fd=[f_Fd; f_Fd; f_Fd]; %vector, 12 entries
f_NADPH =   (n_NADPH - nRings.*(1/3))./n_TOT; F_NADPH=[f_NADPH; f_NADPH; f_NADPH]; %vector, 12 entries
%These are all the relative fluxes to yield BPs 0-3: fW, fx, fFd, fNADPH
%There is no gamma (substrate-derived H) in N. maritimus

% Set the upper and lower bounds for alphas and lambda
alphaW=      0.9;    %Based on Wijker 2019 and simple extrapolation of our BP data to intercept; is also predicted >0.9 if let this parameter float 
lb=[0.01, 0.01, 0.01, 0.01, 0.01];  %lower bounds of x1,..., x5
ub=[1, 1, 1, 0.5, 1];               %upper bounds of x1,..., x5
x0=[0.5, 0.5, 0.5, 0.25, 0.5];      %starting values of x1,..., x5
% alpha_Fd, alpha_NADPH, alpha_E, lambda, alpha_TH = x1,..., x5

%% Input data and run for N. maritimus
Td=[30.8  46.2  92.5];  % Td (avg doubling times, in hours) 
Tdmin=20; Tdmax=120;    % set limits of organism growth rate
fTd=(Td-Tdmin)./(Tdmax-Tdmin); FTd=repelem(fTd,4)'; %vector, 12 entries, fractional limitation of growth rate, XTd
SMOW=1.5576e-4;         % d2H SMOW
rW=1.4901e-4; RW=repelem(rW,12)'; % Isotope ratio of water for Td 30.8, 46.2, 92.5; 4 BPs each

%data, d2H values averaged for 3 growth rates, 4 BPs
%all values of Y for N. mari (solutions for solver)
davg_Nmari= [-310.55; -307.07; -299.36; -297.05; -306.73; -299.05; -292.61; -288.57; -303.60; -299.42; -285.14; -283.98]; %put in a single vector      
Y=(davg_Nmari./1000+1).*SMOW./RW - F_W.*alphaW; 

%set up coefficients, A 
A=[F_Fd F_NADPH (1-FTd)]; 

%equation is RBP/RW-fW*aW = fFd*aFd+fNADPH*aNADPH*aE*aTH/(aTH + lambda(1-Xtd)(1-aTH))
%calculate Y=a1*x1 + a2*x2*x3*x5/(x5 + a3*x4*(1-x5))

%Simulated annealing -- less biased to boundary conditions than non-linear least squares  
numit=10000;        %10,000 trials; Monte Carlo estimation of free parameters
value=(100);        %placeholder  
results=zeros(numit*0.2, 5);
f=zeros(numit*0.2, 1);
i=0;
 
for i=1:numit
funscalar = @(x) sum(abs((x(1).*A(:,1) + A(:,2).*x(2).*x(3).*x(5)./(x(5)+A(:,3).*x(4).*(1-x(5)))) - Y)); 
[x, fval, exitFlag, output] = simulannealbnd(funscalar, x0, lb, ub)
value=x;
    results(i,:)=value;
    f(i,:)=fval;
    clear value 
end

combine=[(f), (results)];
sort=sortrows(combine,1);
best_Nmari=sort(1:1000,:);  %Keep only the best 1000 results (10% of results)
%disp('[order in best_Nmari(:, [2:6]) is alphaFd  alphaNADPH  alphaE  lambda  alphaTH]')

disp('This approach establishes that alpha_Fd is narrowly constrained')
NmariALPHAF_mean=mean(best_Nmari(:,2))  
NmariSigF=std(best_Nmari(:,2)) 

disp('The product of alphaNADPH x alphaE = 0.669 (always, i.e., +/- 0.003)')
mean_aNxaE=mean(best_Nmari(:, 3).*best_Nmari(:,4))
sig_anxaE=std(best_Nmari(:, 3).*best_Nmari(:,4))

disp('Keep only the rows where alphaTH is >= 0.22 and <= 0.57; limits of THs from Wijker 2019')
disp('Find mean and std of lambda within these constraints')
disp('The full range of allowed lambda is 0.01-0.06 or possibly 0.07')
iTH=find(best_Nmari(:, 6)>=0.22 & best_Nmari(:, 6)<=0.57);  %index alphaTH between 0.22 & 0.57
L_best_Nmari=best_Nmari(iTH,:); 
NmariL_mean=mean(L_best_Nmari(:,5)) 
NmariSigL=std(L_best_Nmari(:,5))
minLambda=min(L_best_Nmari(:,5))
maxLambda=max(L_best_Nmari(:,5))

%% Plot histogram of the best fit options

figure(1)
subplot(1,3,1) %alpha Fd, column 2
bin_spacing = [0.0 0.0:0.01:0.5 0.5]; %Set bin edge values
hold on
histogram(best_Nmari(:,2),bin_spacing, 'FaceColor',[0.5 0.5 0.5], 'FaceAlpha', .6); %Create histogram plot  
[m_n2,s_n2]=normfit(best_Nmari(:,2));
title('Model Output Histogram')
xlabel('alphaFd')

subplot(1,3,2) %alphaE vs. alphaNADPH
hold on
plot(best_Nmari(:,3),best_Nmari(:,4), 'ko')
line=0.669./(0.6:0.001:1);  %best fit, i.e., the product of alphaE x alphaNADPH = 0.669
plot(0.6:0.001:1, line, 'LineWidth', 2)
axis([0.6 1 0.6 1]);
xlabel('alphaNADPH')
ylabel('alphaE')

subplot(1,3,3) %lambda, column 5
bin_spacing = [0.0 0.0:0.01:0.5 0.5]; %Set bin edge values
hold on
histogram(L_best_Nmari(:,5),bin_spacing, 'FaceColor',[0.5 0.5 0.5], 'FaceAlpha', .6); %Create histogram plot  
[m_n5,s_n5]=normfit(best_Nmari(:,5));
title('Model Output Histogram')
xlabel('lambda')

 