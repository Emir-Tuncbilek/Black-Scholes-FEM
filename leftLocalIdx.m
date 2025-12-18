function idxL = leftLocalIdx(ne, pDeg, pType)
    [N, ~] = shape(-1, ne, pDeg, pType);   % xi = -1
    [~, idxL] = max(N);
end
