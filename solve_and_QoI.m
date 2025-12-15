function Qh = solve_and_QoI(SMin,SMax,nEls,pVal,t0,T,Nt_in,K,r,Sstar)

% -----------------------------
% Mesh
% -----------------------------
[~,nodes,connect,~,~,~] = mesh(SMin,SMax,nEls);

% -----------------------------
% FEM metadata
% -----------------------------
pDeg  = pVal*ones(nEls,1);
pType = ones(nEls,1);                 % 1 = Lagrange
[elDof,dFreedom] = dof(nEls,pDeg,connect);
m = max(max(dFreedom));

% -----------------------------
% Quadrature
% -----------------------------
[xiQ,wQ] = gQuad;

% -----------------------------
% Assemble M, A (as in main)
% -----------------------------
[M,A,~] = time_element(nEls, nodes, connect, xiQ, wQ, pDeg, pType, elDof, dFreedom, t0);

% -----------------------------
% Time vars (as in main)
% -----------------------------
Nt = Nt_in;
Nt = Nt - 1;                      % <-- matches your main
dt = (T - t0) / Nt;

% -----------------------------
% Initial condition U0 (as in main: local fit via Nmat\fvec)
% -----------------------------
U0 = zeros(m,1);

for e = 1:nEls
    x1 = nodes(connect(e,1));
    x2 = nodes(connect(e,2));
    nLoc = elDof(e);

    xiSamples = linspace(-1, 1, nLoc);

    Nmat = zeros(nLoc, nLoc);
    fvec = zeros(nLoc, 1);

    for k = 1:nLoc
        xi = xiSamples(k);

        [N, ~] = shape(xi, e, pDeg, pType);
        Nmat(k,:) = N(:).';

        Sx = (x1*(1 - xi) + x2*(1 + xi))/2;
        fvec(k) = max(Sx - K, 0);
    end

    u_loc = Nmat \ fvec;

    for a = 1:nLoc
        gi = dFreedom(e,a);
        U0(gi) = u_loc(a);
    end
end

% -----------------------------
% Crank–Nicolson loop (as in main)
% -----------------------------
U = zeros(m, Nt+1);
U(:,1) = U0;

F = zeros(m,1);                    % as in main

L = M + 0.5*dt*A;
R = M - 0.5*dt*A;

for n = 1:Nt
    tau_n  = n*dt/2;               % <-- matches your main
    bn_raw = R*U(:,n) + dt/2*F;    % <-- matches your main (F=0 anyway)

    [L_bc, bn] = applyBC_time(L, bn_raw, dFreedom, SMax, K, r, tau_n);
    U(:,n+1) = L_bc \ bn;
end

% -----------------------------
% QoI at final time
% -----------------------------
Qh = QoI_point_from_vec(U(:,end), Sstar, nodes, connect, elDof, dFreedom, pDeg, pType);

end
