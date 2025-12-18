function results = convergenceStudy_MMS(nEls_list, Nt, p, SMin, SMax, t0, T)
% convergenceStudy_MMS
% Runs an h-convergence study for your MMS time-dependent 1D FEM code.
%
% Requires your existing functions:
%   mesh, dof, gQuad, time_element, buildU0, applyBC_MMS, errors_time,
%   leftLocalIdx, rightLocalIdx
%
% Output:
%   results: table with nEls, h, L2agg, H1agg, rateL2, rateH1

    if nargin < 1 || isempty(nEls_list), nEls_list = [5, 10 20 40 80 160 320]; end
    if nargin < 2 || isempty(Nt),        Nt = 200; end
    if nargin < 3 || isempty(p),         p  = 2;   end
    if nargin < 4 || isempty(SMin),      SMin = 0; end
    if nargin < 5 || isempty(SMax),      SMax = 50; end
    if nargin < 6 || isempty(t0),        t0 = 0; end
    if nargin < 7 || isempty(T),         T  = 1; end

    Nt = Nt - 1;                 % keep your convention
    dt = (T - t0) / Nt;
    time = linspace(t0, T, Nt+1);

    nRuns = numel(nEls_list);

    hh    = zeros(nRuns,1);
    L2agg = zeros(nRuns,1);
    H1agg = zeros(nRuns,1);

    % For error-vs-time plots
    L2t_cell = cell(nRuns,1);
    H1t_cell = cell(nRuns,1);

    fprintf("\n==============================\n");
    fprintf(" MMS h-Convergence Study (p=%d)\n", p);
    fprintf(" Domain: [%.3g, %.3g], T=%.3g, Nt=%d, dt=%.3e\n", SMin, SMax, T, Nt+1, dt);
    fprintf("==============================\n\n");

    for i = 1:nRuns
        nEls = nEls_list(i);
        hh(i) = (SMax - SMin) / nEls;

        fprintf("---- Run %d/%d: nEls=%d (h=%.4e) ----\n", i, nRuns, nEls, hh(i));

        % Solve + errors
        [L2_t, H1_t, L2_i, H1_i] = solve_one_MMS_and_errors(nEls, p, SMin, SMax, t0, T, Nt);

        L2agg(i) = L2_i;
        H1agg(i) = H1_i;

        L2t_cell{i} = L2_t(:);
        H1t_cell{i} = H1_t(:);

        fprintf("   L2agg = %.6e\n", L2agg(i));
        fprintf("   H1agg = %.6e\n\n", H1agg(i));
    end

    % Compute convergence rates (slope between consecutive meshes)
    rateL2 = nan(nRuns,1);
    rateH1 = nan(nRuns,1);

    for i = 2:nRuns
        rateL2(i) = log(L2agg(i-1)/L2agg(i)) / log(hh(i-1)/hh(i));
        rateH1(i) = log(H1agg(i-1)/H1agg(i)) / log(hh(i-1)/hh(i));
    end

    % Build neat table
    results = table(nEls_list(:), hh, L2agg, rateL2, H1agg, rateH1, ...
        'VariableNames', {'nEls','h','L2_error','L2_rate','H1_error','H1_rate'});

    fprintf("\n==============================\n");
    fprintf(" Final Convergence Table\n");
    fprintf("==============================\n");
    disp(results);

% -----------------------------
% Plot error vs time (all meshes) - robust to varying lengths
% -----------------------------

figure; hold on; grid on;
for i = 1:nRuns
    e = L2t_cell{i};
    % build matching time vector of same length
    t_plot = linspace(t0 + dt, T, numel(e));
    plot(t_plot, e(:), 'LineWidth', 1.8);
end
xlabel('Time');
ylabel('L^2 error');
title('L^2 error vs time (all meshes)');
legend(compose('nEls=%d', nEls_list), 'Location','best');
hold off;

figure; hold on; grid on;
for i = 1:nRuns
    e = H1t_cell{i};
    t_plot = linspace(t0 + dt, T, numel(e));
    plot(t_plot, e(:), 'LineWidth', 1.8);
end
xlabel('Time');
ylabel('H^1 error');
title('H^1 error vs time (all meshes)');
legend(compose('nEls=%d', nEls_list), 'Location','best');
hold off;


    % -----------------------------
    % Plot aggregate error vs h (log-log)
    % -----------------------------
    figure; grid on; hold on;
    loglog(hh, L2agg, '-o', 'LineWidth', 2);
    loglog(hh, H1agg, '-s', 'LineWidth', 2);
    set(gca,'XDir','reverse');
    xlabel('h = (S_{max}-S_{min})/nEls');
    ylabel('Aggregate error');
    title(sprintf('Aggregate errors vs h (p=%d)', p));
    legend('L^2 aggregate','H^1 aggregate','Location','best');
    hold off;
end


function [L2_t, H1_t, L2agg, H1agg] = solve_one_MMS_and_errors(nEls, p, SMin, SMax, t0, T, Nt)
% One full solve on a given mesh, then compute time-dependent + aggregate errors.

    % Mesh
    [~, nodes, connect, nB, bEls, bPts] = mesh(SMin, SMax, nEls);

    % FEM setup
    pDeg  = p * ones(nEls,1);
    pType = ones(nEls,1);     % 1 = Lagrangian

    [elDof, dFreedom] = dof(nEls, pDeg, connect);

    [xiQ, wQ] = gQuad;

    % Time
    dt = (T - t0) / Nt;
    time = linspace(t0, T, Nt+1);

    m = max(max(dFreedom));
    U = zeros(m, Nt+1);

    % MMS initial condition
    u0_fun = @(S) (1 - t0) .* S .* (SMax - S);
    U(:,1) = buildU0(nEls, nodes, connect, elDof, dFreedom, pDeg, pType, u0_fun);

    % Boundary DOFs (use your robust local-index helpers)
    idxL = leftLocalIdx(1, pDeg, pType);
    idxR = rightLocalIdx(nEls, pDeg, pType);
    cL = dFreedom(1, idxL);
    cR = dFreedom(nEls, idxR);

    % March in time (your MMS CN scheme)
    for n = 1:Nt
        tau_n   = (n-1)*dt;
        tau_np1 = (n)*dt;

        % enforce BCs on current state
        U(cL,n) = 0;
        U(cR,n) = 0;

        [M, A, Fn]   = time_element(nEls, nodes, connect, xiQ, wQ, pDeg, pType, elDof, dFreedom, tau_n);
        [~, ~, Fnp1] = time_element(nEls, nodes, connect, xiQ, wQ, pDeg, pType, elDof, dFreedom, tau_np1);

        L = M + 0.5*dt*A;
        R = M - 0.5*dt*A;

        bn_raw = R*U(:,n) + 0.5*dt*(Fn + Fnp1);

        [L_bc, bn] = applyBC_MMS(L, bn_raw, dFreedom, pDeg, pType);

        U(:,n+1) = L_bc \ bn;

        % enforce again after solve
        U(cL,n+1) = 0;
        U(cR,n+1) = 0;
    end

    % Errors (your function)
    [L2_t, H1_t, L2agg, H1agg] = errors_time( ...
        nEls, nodes, connect, elDof, dFreedom, pDeg, pType, U, time);
end
