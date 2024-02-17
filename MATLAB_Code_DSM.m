clear all
close all
clc

tic

%% settings for text, legends and labels
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Import Data
N = 100;
energy_need = (readmatrix('energy_need.csv'))';
p1 = (readmatrix('p1.csv'))';
p2 = (readmatrix('p2.csv'))';
xbar = (readmatrix('xbar.csv'))';

%% Exercise 1.e 

% Build Jacobian
h=24;
N = 100;
JF = eye(h);
for i = 1:7
    JF = [JF,JF;JF,JF];
end
JF = JF(1:N*h,1:N*h)+eye(N*h);
for i = 1:N
    JF(((i-1)*h+1):i*h,:) = p1(i)*JF(((i-1)*h+1):i*h,:);
end

% Compute μ, L and γ_max
% eigvals_JFT = eig(JF');
% eigvals_JF = eig(JF);
% mu = min(eigvals_JFT);
% L = max(eigvals_JF);
% gamma_max = 2*mu/L^2;
% tau_opt = 1-mu^2/L^2;
% gamma_max = 0.012995075352174;
% gamma = 0.5*gamma_max;
gamma = 1.5604;

% Iteration Parameters
max_it = 200;
A = eye(N*h)-gamma*JF;
% rho = max(abs(eig(A)));
b = gamma*kron(p2,ones(h,1));

%% Initialize optimmization problem framework
ub = reshape(xbar', [], 1);

xProj = sdpvar(N*h,1);
x = sdpvar(N*h, 1);  % parameters containing information before projection
objective = (xProj-x)'*(xProj-x);
constraints = [xProj >= 0, xProj <= ub];
for i = 1:N
    constraints = [constraints, sum(xProj(((i-1)*h+1):i*h)) == sum(energy_need(i,:))];
end

opts = sdpsettings('verbose', 0);
projector = optimizer( ...
    constraints, objective, opts, ...
    x, ...  % input parameter of the optimization
    xProj ...  % result
);

%%
% Random Initialisation: x0 ∈ K
x0 = 20*rand(N*h,1);
[x_it, errorcode] = projector(x0);

% Best-response algorithm 
% norm_b = norm(b);
for it = 1:max_it
    x_prev = x_it;
    % Update
    x_it = A*x_it-b;
    % Projection
    [x_it, errorcode] = projector(x_it);
end



%% Print of NE and file generation
xstar = [];
for i = 1:N
    xstar(:,i) = x_it(((i-1)*h+1):i*h);
end
disp('    ');
disp('The Nash Equilibrium of the game is: ');
disp(xstar);

time_spent = toc;

% file generation
% filename = 'Graziano_NE.csv';
% writematrix(xstar,filename);
writematrix(xstar, sprintf("x_star_%03.0f.csv", time_spent));


%% 1.f - Plot
figure(1)
plot(1:h,sum(energy_need))
hold on
plot(1:h,sum(xstar'))
plot(1:h,sum(xbar))
legend("$e(h)\ =\ \sum\limits_{i = 1}^N{e_i^h}$", ...
    '$x^*(h)\ =\ \sum\limits_{i = 1}^N{x_i^{*,h}}$', ...
    '$\bar{x}(h)\ =\ \sum\limits_{i = 1}^N{\bar{x}_i^{h}}$', ...
    "Location",'west','FontSize',23);
xlabel('$h$','FontSize',25)
ylabel('$e(h)\ \ \land\ \ x^*(h)\ \ \land\ \ \bar{x}(h)$','FontSize',25)
title('$Aggregate\ Nominal\ and\ NE\ energy\ consumption\ vs\ time$','FontSize',26)
grid on
xlim([1 24])
set(figure(1),'Position', get(figure(1), 'Position')+[-200, -170, 400, 400]);
hold off


%% The code runs up to this point. The code for the checks is below.

return

%% Exercise 1.c - Assess Existence and Uniqueness of NE
% non-empty strategy spaces
bool = 1;
empty_spaces = [];
for i = 1:N
    di = sum(energy_need(i,:));
    max_capability = sum(xbar(i,:));
    if (di > max_capability | di < 0) 
        bool = 0;
        empty_spaces = [empty_spaces;i];
    end
end
disp(' ');
if bool==1
    disp('All strategy spaces are non-empty!');
else
    disp('The following players have empty strategy spaces:');
    disp(empty_spaces');
end
disp(' ');

% convexity of cost functions for fixed x_{-i}
warning = find(p1<0);
if warning ~= []
    disp('The following players have a personalised aggregate energy consumption parameter p_{i,1} < 0:');
    disp(warning);
else 
    disp('All players have a personalised aggregate energy consumption parameter p_{i,1} ≥ 0. '); 
end
disp(' ');

% Uniqueness of the NE
warning = find(p1<0);
if warning ~= []
    disp('The following players have a personalised aggregate energy consumption parameter p_{i,1} ≤ 0:');
    disp(warning);
else 
    disp('All players have a personalised aggregate energy consumption parameter p_{i,1} > 0. '); 
end
disp(' ');

%% Exercise 1.e - Condition on the eigenvalues of JF and JF^T
eigvals_JFT = eig(JF');
eigvals_JF = eig(JF);
disp("    ");
if all(eigvals_JFT > 0)
    disp('The eigenvalues of JF^T are all strictly positve!');
else
    disp('The eigenvalues of JF^T are NOT all strictly positve!');
end
disp("    ");
if all(eigvals_JF > 0)
    disp('The eigenvalues of JF are all strictly positve!');
else
    disp('The eigenvalues of JF are NOT all strictly positve!');
end
disp("    ");




