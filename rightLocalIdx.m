function idxR = rightLocalIdx(ne, pDeg, pType)
    [N, ~] = shape(1, ne, pDeg, pType);    % xi = +1
    [~, idxR] = max(N);
end
