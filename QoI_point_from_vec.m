function q = QoI_point_from_vec(uvec, Sstar, nodes, connect, elDof, dFreedom, pDeg, pType)
nEls = size(connect,1);

% locate element
e = -1;
for k = 1:nEls
    x1 = nodes(connect(k,1));
    x2 = nodes(connect(k,2));
    if (Sstar >= x1) && (Sstar <= x2)
        e = k; break;
    end
end
if e < 0, error('Sstar outside domain'); end

x1 = nodes(connect(e,1));
x2 = nodes(connect(e,2));
h  = x2 - x1;
xi = (2*Sstar - (x1 + x2))/h;

[N,~] = shape(xi, e, pDeg, pType);

nLoc = elDof(e);
loc  = dFreedom(e,1:nLoc);

q = 0;
for i = 1:nLoc
    q = q + uvec(loc(i))*N(i);
end
end
