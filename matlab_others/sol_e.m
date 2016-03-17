function p = sol_e(XX,beta,x1,betat,h1,h2)
    global S;
    p = XX*0;
    
    pm = x1/h1^3 *2 *(h2-h1) * (S-betat);
    
    II = (XX > 0) & (XX < x1);
    p(II) = XX(II)/x1 * pm;
    
    II = (XX > x1) & (XX < beta);
    p(II) = (beta-XX(II))/(beta-x1) * pm;
end