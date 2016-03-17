function z = dimple2(xx)
    global dep l1;
    z = 0*xx;
    z(xx<l1) = dep;
end