clear all
clc
format long
syms x n;
%% init
f = exp(x)-x-1;
df = diff(f)
d2f = diff(df)
errors = zeros(3,3);
TOL = 1.e-5;
x0 = 1;
n = 5;
%% Newton Method
% f(x) = e^x - x-1 = 0 has multiple root so Rate of Convergence =1
for i=1:n
    fi = vpa(subs(f,x,x0),12);
    dfi = vpa(subs(df,x,x0),12);
    xi = x0 - fi/dfi;
    sol(1,i) = sym2poly(xi);
    error = abs(xi-x0);
    errors(1,i) = sym2poly(error);
    if i==1
        fprintf('%dst newton solution is %.10f   ',i,sol(1,i))
        fprintf('error is %.10f   ' ,errors(1,i))
        fprintf('abs error is %.10f\n',abs(sol(1,i)))
    elseif i==2
        fprintf('%dnd newton solution is %.10f   ',i,sol(1,i))
        fprintf('error is %.10f   ' ,errors(1,i))
        fprintf('abs error is %.10f\n',abs(sol(1,i)))
    else
        fprintf('%dth newton solution is %.10f   ',i,sol(1,i))
        fprintf('error is %.10f   ' ,errors(1,i))
        fprintf('abs error is %.10f\n',abs(sol(1,i)))
    end
    x0 = xi;
end
fprintf('\n')
%% Newton Method Convergence Acceleration
% make Rate of Convergence = 2
x0 = 1;
for i=1:n
    fi = vpa(subs(f,x,x0),12);
    dfi = vpa(subs(df,x,x0),12);
    xi = x0 - 2*fi/dfi;
    sol(2,i) = sym2poly(xi);
    error = abs(xi-x0);
    errors(2,i) = sym2poly(error);
    x0 = xi;
    if i==1
        fprintf('%dst newton accelerating1 solution is %.10f   ',i,sol(2,i))
        fprintf('error is %.10f   ' ,errors(2,i))
        fprintf('abs error is %.10f\n',abs(sol(2,i)))
    elseif i==2
        fprintf('%dnd newton accelerating1 solution is %.10f   ',i,sol(2,i))
        fprintf('error is %.10f   ' ,errors(2,i))
        fprintf('abs error is %.10f\n',abs(sol(2,i)))
    else
        fprintf('%dth newton accelerating1 solution is %.10f   ',i,sol(2,i))
        fprintf('error is %.10f   ' ,errors(2,i))
        fprintf('abs error is %.10f\n',abs(sol(2,i)))
    end
end
fprintf('\n')
%% Convergence Acceleration2
% Rate of Convergence = 2
x0 = 1;
for i=1:n
    fi = vpa(subs(f,x,x0),12);
    dfi = vpa(subs(df,x,x0),12);
    d2fi = vpa(subs(d2f,x,x0),12);
    xi = x0 - (fi*dfi)/(dfi^2-fi*d2fi);
    sol(3,i) = sym2poly(xi);
    error = abs(xi-x0);
    errors(3,i) = sym2poly(error);
    x0 = xi;
    
    if i==1
        fprintf('%dst newton accelerating2 solution is %.10f   ',i,sol(3,i))
        fprintf('error is %.10f   ' ,errors(3,i))
        fprintf('abs error is %.10f\n',abs(sol(3,i)))
    elseif i==2
        fprintf('%dnd newton accelerating2 solution is %.10f   ',i,sol(3,i))
        fprintf('error is %.10f   ',errors(3,i))
        fprintf('abs error is %.10f\n',abs(sol(3,i)))
    else
        fprintf('%dth newton accelerating2 solution is %.10f   ',i,sol(3,i))
        fprintf('error is %.10f   ' ,errors(3,i))
        fprintf('abs error is %.10f\n',abs(sol(3,i)))
    end
end

