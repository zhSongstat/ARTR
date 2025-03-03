function x = my_prox_l1(x,t) 
tq = t * 1;
s  = 1 - min( tq./abs(x), 1 );
x  = x .* s;
end