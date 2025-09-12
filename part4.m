%Quadratic function with root at the minimum
function [f_val,dfdx] = test_function02(x)
    global input_list;
    input_list(:,end+1) = x;
    f_val = (x-30.879).ˆ2;
    dfdx = 2*(x-30.879);
end

%Example sigmoid function
function [f_val,dfdx] = test_function03(x)
    global input_list;
    input_list(:,end+1) = x;
    a = 27.3; b = 2; c = 8.3; d = -3;
    H = exp((x-a)/b);
    dH = H/b;
    L = 1+H;
    dL = dH;
    f_val = c*H./L+d;
    dfdx = c*(L.*dH-H.*dL)./(L.ˆ2);
end