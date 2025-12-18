function [L_bc, b_bc] = applyBC_MMS(L, b, dFreedom, pDeg, pType)
    L_bc = L;
    b_bc = b;

    nEls = size(dFreedom,1);

    idxL = leftLocalIdx(1,    pDeg, pType);
    idxR = rightLocalIdx(nEls, pDeg, pType);

    cL = dFreedom(1,    idxL);
    cR = dFreedom(nEls, idxR);

    [L_bc, b_bc] = pinDirichlet(L_bc, b_bc, cL, 0);
    [L_bc, b_bc] = pinDirichlet(L_bc, b_bc, cR, 0);
end

function [A, rhs] = pinDirichlet(A, rhs, dof, val)
    rhs = rhs - A(:,dof) * val;
    A(:,dof) = 0;
    A(dof,:) = 0;
    A(dof,dof) = 1;
    rhs(dof) = val;
end
