 function main
clc
close all
%*******************main code for the 1D finite element solver*************

%**************************build mesh**************************************

SMin=0;
SMax=50;
nEls=20;
[nN,nodes,connect,nB,bEls,bPts]=mesh(SMin,SMax,nEls);


%************************** Time vars **************************************

t0 = 0;                                 % Begin
T = 1;                                  % End
Nt = 200;                                 % Time steps
Nt = Nt - 1;
dt = (T-t0) / Nt ;                       % Step size
% time = linespace(t0, T, Nt+1);          % List of times

%*******************degrees of freedom*************************************

pDeg=zeros(nEls,1);
pDeg(:,1)=1;                            %set polynomial degree for each element
%pDeg=[1,2,2,1]';                       %could set each element individually
pType=zeros(nEls,1);
pType(:,1)=1;                           %set element type: 1=Lagrangian, 2=hierarchical

[elDof,dFreedom]=dof(nEls,pDeg,connect);

%*****************build element matrices k and vectors f*******************

%get quadrature points and weights
[xiQ,wQ]=gQuad;

%build k and f for each element; sum to find K and F
% [K,F]=element(nEls,nodes,connect,xiQ,wQ,pDeg,pType,elDof,dFreedom);

% build m, a and f for each; sum to find M, A and F
[M, A, F] = time_element(nEls, nodes, connect, xiQ,wQ,pDeg, pType, elDof, dFreedom, t0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD POINT LOAD f0*delta(x - x0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x0 = 0;        % [CHANGE ME] location of point load
% f0 = 0;        % [CHANGE ME] magnitude of point load

% --- Step 1: Find the element that contains x0
% elem = -1;
% 
% for e = 1:nEls
%     x1 = nodes(connect(e,1));
%     x2 = nodes(connect(e,2));
%     if x0 >= x1 && x0 <= x2
%         elem = e;
%         break;
%     end
% end
% 
% if elem == -1
%     error('Point load applied outside of the domain.');
% end

% --- Step 2: Map x0 to reference coordinate xi0
% x1 = nodes(connect(elem,1));
% x2 = nodes(connect(elem,2));
% h  = x2 - x1;
% xi0 = (2*x0 - (x1 + x2)) / h;   % inverse mapping

% --- Step 3: Evaluate shape functions at xi0
% [N, ~] = shape(xi0, elem, pDeg, pType);

% --- Step 4: Add point load to the global RHS
% nLoc = elDof(elem);
% loc  = dFreedom(elem, 1:nLoc);  % global DOF indices
% 
% for i = 1:nLoc
%     F(loc(i)) = F(loc(i)) - f0 * N(i);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%*****************add boundary conditions and solve************************

%set boundary condition type (set values within boundaryC.m)
% bType=zeros(2,1);   %boundary condition types; use 1 for Dirichlet,
% bType(1)=1;         %2 for Neumann, 3 for Robin
% bType(2)=2;
% 
% [K,F]=boundaryC(nB,bEls,bPts,bType,dFreedom,K,F);
% u=K\F;


%***************** add time boundary conditions ************************
% bType=zeros(2,1);   %boundary condition types; use 1 for Dirichlet,
% bType(1)=1;         %2 for Neumann, 3 for Robin
% % bType(2)=2;
% 
% [M, F]=boundaryC(nB,bEls,bPts,bType,dFreedom,M,F);

%****************************** solve **********************************

K = 0.5 * SMax;            % strike
r = 1;                     % risk

m = max(max(dFreedom));
U = zeros(m, Nt+1);

U0 = zeros(m, 1);

% Compute U0
for e = 1:nEls
    x1 = nodes(connect(e,1));
    x2 = nodes(connect(e,2));

    nLoc = elDof(e);

    % Reference sample points for this element
    xiSamples = linspace(-1, 1, nLoc);

    Nmat = zeros(nLoc, nLoc);   % local shape matrix
    fvec = zeros(nLoc, 1);      % local target (payoff) values

    for k = 1:nLoc
        xi = xiSamples(k);

        % shape functions at this reference point
        [N, ~] = shape(xi, e, pDeg, pType);

        % fill row of Nmat
        Nmat(k, :) = N(:).';

        % map xi -> physical S
        Sx = (x1*(1 - xi) + x2*(1 + xi))/2;

        % payoff at this point
        fvec(k) = max(Sx - K, 0);
    end

    % Solve for local DOFs that best fit the payoff at sample points
    u_loc = Nmat \ fvec;

    % Scatter to global U0
    for a = 1:nLoc
        gi = dFreedom(e, a);
        U0(gi) = u_loc(a);
    end
end

%====================================%
%======= Manufactured Solution ======%
%====================================%
% Time stepping (Crank–Nicolson in tau)
% u0_fun = @(S) (1-t0) .* S .* (SMax - S);   % MMS at tau=t0 (usually 0)
% U0 = buildU0(nEls, nodes, connect, elDof, dFreedom, pDeg, pType, u0_fun);

% Initialize the global solution vector U with the initial condition U0
U(:,1) = U0;

%====================================%
%====== Black-Scholes Solution ======%
%====================================%

F = zeros(m,1);          % assume no source term for Black–Scholes

L = M + 0.5*dt*A;
R = M - 0.5*dt*A;

for n = 1:Nt
    tau_n   = n*dt/2;        % or time(n) if you stored it
    bn_raw  = R*U(:,n) + dt/2*F;

    % Apply Dirichlet BCs for this time:
   [L_bc, bn] = applyBC_time(L, bn_raw, dFreedom, SMax, K, r, tau_n);


   U(:,n+1) = L_bc \ bn;
end

%====================================%
%======= Manufactured Solution ======%
%====================================%

% idxL = leftLocalIdx(1, pDeg, pType);
% idxR = rightLocalIdx(nEls, pDeg, pType);
% 
% cL = dFreedom(1, idxL);
% cR = dFreedom(nEls, idxR);
% 
% for n = 1:(Nt)
%     tau_n   = (n-1)*dt;
%     tau_np1 = (n)*dt;
% 
%     % enforce BCs on current state BEFORE R*U(:,n)
%     U(cL,n) = 0;
%     U(cR,n) = 0;
% 
%     [M, A, Fn]   = time_element(nEls, nodes, connect, xiQ, wQ, pDeg, pType, elDof, dFreedom, tau_n);
%     [~, ~, Fnp1] = time_element(nEls, nodes, connect, xiQ, wQ, pDeg, pType, elDof, dFreedom, tau_np1);
% 
%     L = M + 0.5*dt*A;
%     R = M - 0.5*dt*A;
% 
%     bn_raw = R*U(:,n) + 0.5*dt*(Fn + Fnp1);
% 
%     [L_bc, bn] = applyBC_MMS(L, bn_raw, dFreedom, pDeg, pType);
% 
%     U(:,n+1) = L_bc \ bn;
% 
%     % enforce again after solve
%     U(cL,n+1) = 0;
%     U(cR,n+1) = 0;
% end


u = U(:,end-1);   % this is U(T,·) in tau, i.e. price at t=0


%**************post-processing*********************************************
figure(1)
postProc(nEls,nodes,connect,elDof,dFreedom,pDeg,pType,u);
% figure(2)
% postProc(nEls,nodes,connect,elDof,dFreedom,pDeg,pType,U0);


time = linspace(0, T, Nt+1);

plot3D(nEls, nodes, connect, elDof, dFreedom, pDeg, pType, U, time);


% compute and display the errors:
% Calculate the exact solution for comparison
% [uL2, uH1] = errors(nEls, nodes, connect, elDof, dFreedom, pDeg, pType, u);
% [uL2, uH1] = numericalErrors(5, 200, 2, 2);   % compare nEls=10 to a 200-element reference

% fprintf('L2 error = %.6e\n', uL2);
% fprintf('H1 error = %.6e\n', uH1);

[L2_t, H1_t, L2agg, H1agg] = errors_time( ...
    nEls, nodes, connect, elDof, dFreedom, pDeg, pType, U, time);

fprintf('L2 error = %.6e\n', L2agg);
fprintf('H1 error = %.6e\n', H1agg);
fprintf('Min time: %.6e\n', min(time));
fprintf('First time: %6e\n', time(2));

% Build a time vector that matches L2_t/H1_t length (main-only fix)
t_err = t0 + (1:length(L2_t)) * dt;   % assumes first error corresponds to t0+dt

figure
hold on
plot(t_err, L2_t, 'b-',  'LineWidth', 2)
plot(t_err, H1_t, 'r--', 'LineWidth', 2)
grid on
xlabel('Time')
ylabel('Error norm')
title('Error norms vs time')
legend('L^2 error', 'H^1 error', 'Location', 'best')
hold off


% Nelem = [10 20 40 80 160 320];
% results = convergenceStudy_MMS(Nelem, 200, 2, 0, 50, 0, 1);

%% ================================
% QoI convergence study (h + p)
% ================================

Sstar = 1.30*K;

nEls_list = [10, 20, 40, 80, 160, 320]; % 6 meshes
p_list    = [1 2];                      % p-enrichment

Nt_test = 400;                          % keep time error small and fixed

% Reference solution (fine)
nEls_ref = 1000;
Nt_ref   = 200;
p_ref    = 2;

fprintf("\n=== Computing Qref (reference) ===\n");
Qref = solve_and_QoI( ...
    SMin,SMax,nEls_ref,p_ref,t0,T,Nt_ref, ...
    K,r,Sstar);
fprintf("Qref = %.12e\n\n", Qref);

% Sweep and print tables
for p = p_list
    Qh  = zeros(numel(nEls_list),1);
    Err = zeros(numel(nEls_list),1);
    hh  = zeros(numel(nEls_list),1);

    fprintf("=== QoI convergence for p=%d ===\n", p);

    for i = 1:numel(nEls_list)
        nEls_i = nEls_list(i);
        hh(i) = (SMax - SMin)/nEls_i;

        Qh(i) = solve_and_QoI( ...
            SMin,SMax,nEls_i,p,t0,T,Nt_test, ...
            K,r,Sstar);

        Err(i) = abs(Qh(i) - Qref);

        fprintf("nEls=%4d  h=%.4e  Qh=%.12e  |Qh-Qref|=%.3e\n", ...
            nEls_i, hh(i), Qh(i), Err(i));
    end

    % Plot error vs h (log-log)
    figure(100); hold on; grid on;
    loglog(hh, Err, '-o', 'LineWidth', 2);

    TBL = table(nEls_list(:), hh, Qh, Err, ...
        'VariableNames', {'nEls','h','Qh','AbsError'});
    disp(TBL);
end

fprintf('\n\n');

figure(100);
set(gca,'XDir','reverse');
xlabel('h = (S_{max}-S_{min})/nEls');
ylabel('|Q_h - Q_{ref}|');
title(sprintf('QoI convergence at S* = %.4g, T = %.4g', Sstar, T));
legend('p=1','p=2','Location','best');