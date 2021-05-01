clear all
clc
format long
syms x n;

%% Init
%f(x) = x-cos(x)
g = cos(x);
TOL = 1.e-5;
x0 = 1;
%% Fixed point Iterative Method
for i=1:150
    xi = vpa(subs(g,x,x0),10);
    sol(1,i) = sym2poly(xi);
    error(1,i) = sym2poly(abs(xi-x0));
    if error(1,i)<TOL
        n = i;
        break;
    end
    x0 = xi;
end
fprintf('fixed point iterative method repeat %d\n',n)
for i=1:n
    if i==1
        fprintf('%dst error is %.10f   ' ,i,error(1,i))
        fprintf('%dst solution is %.10f \n',i,sol(1,i))
    elseif i==2
        fprintf('%dnd error is %.10f   ' ,i,error(1,i))
        fprintf('%dnd solution is %.10f \n',i,sol(1,i))
    else
        fprintf('%dth error is %.10f   ' ,i,error(1,i))
        fprintf('%dth solution is %.10f \n',i,sol(1,i))
    end
end

%% Steffensen Method
x0 = 1;
n = 1;
for i=1:150
    x3i_2 = vpa(subs(g,x,x0),10);
    x3i_1 = vpa(subs(g,x,x3i_2),10);
    xi = x0-((x3i_2-x0)^2)/(x3i_1-2*x3i_2+x0);
    sol(2,i) = sym2poly(xi);
    error(2,i) = sym2poly(abs(xi-x0));
    if error(2,i)<TOL
        n = i;
        break;
    end
    x0 = xi;
end

fprintf('\nSteffensen Method repeat %d\n',n)
for i=1:n
    if i==1
        fprintf('%dst error is %.10f   ' ,i,error(2,i))
        fprintf('%dst solution is %.10f \n',i,sol(2,i))
    elseif i==2
        fprintf('%dnd error is %.10f   ' ,i,error(2,i))
        fprintf('%dnd solution is %.10f \n',i,sol(2,i))
    else
        fprintf('%dth error is %.10f   ' ,i,error(2,i))
        fprintf('%dth solution is %.10f \n',i,sol(2,i))
    end
end