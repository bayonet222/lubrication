function z = Fbt(t,beta)
    global S h1 h2;
    x=t*S;
    z = S/2*( 1 + (x*h2^3) ./ ( (h2^3-h1^3)*x+beta*h1^3 ) );
end