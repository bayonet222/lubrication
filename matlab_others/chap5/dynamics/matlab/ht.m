function z = ht(t,x)
    global S per;
    z = dimple1(mod(x - S*t,per));
end