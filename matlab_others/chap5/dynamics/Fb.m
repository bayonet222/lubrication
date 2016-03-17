function z=Fb(x,beta)
    global S h1 h2;
    z = S/2*( 1 + (x*h2^3) / ( (h2^3-h1^3)*x+beta*h1^3 ) );
end