function z = dimple1(xx)
    global dep l1;
    dx = xx(2)-xx(1);
    z = 0*xx;
    z(xx <= (l1+dx*0) ) = dep;
end