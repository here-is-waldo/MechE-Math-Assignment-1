function [x_root, exit_flag] = bisection_solver(f_in, L0, R0, max_iter, dx_tol, y_tol)
    
    % check that interval works (diff signs)
    Lf = f_in(L0); Rf = f_in(R0);
    if (Lf<0 && Rf<0) || (Lf>0 && Rf>0)
        disp("Invalid interval. May not contain a root")
    end
    
    % initialize containers
    count = 0;
    x_guess_list = [];

    % assign vals
    L_n = L0;
    R_n = R0;
    M_n = (L_n + R_n)/2;
    M_prev = M_n-2*dx_tol;

    % evaluate at vals
    L_val = f_in(L_n);
    R_val = f_in(R_n);
    M_val = f_in(M_n);

    while abs(M_val) > y_tol ...   # func value isnt close enough to zero
            && abs((M_val-M_prev)) > dx_tol ... % change in guess isn't too small
            && count < max_iter  % within max iteration limit
        
        % update interval
        % temp is arranged so that the negative function
        % value is the first element ie [-, +]
        order = 0;
        if L_val < R_val
            temp = [L_n, R_n];
            order = 1;
        else
            temp = [R_n, L_n];
            order = 2;
        end

        % temp is organized as [n1, n2] such that [f(n1)<0, f(n2)>0]
        % this logic replaces n1 or n2 with M_n according to the sign 
        % of f(M_n) (aka M_val)
        if M_val > 0 
            temp(2) = M_n;
        else
            temp(1) = M_n;
        end
        
        if order == 1
            L_n = temp(1);
            R_n = temp(2);
        elseif order == 2
            L_n = temp(2);
            R_n = temp(1);
        else 
            disp("order variable never updated")
        end

        M_prev = M_n;
        M_n = (L_n + R_n)/2;

        % evaluate at new vals
        L_val = f_in(L_n);
        R_val = f_in(R_n);
        M_val = f_in(M_n);

        count=count+1;
    end

    x_root = M_n;

    % assume success
    exit_flag = 0;

    if abs(M_n-M_prev) > dx_tol && abs(M_val)>y_tol
       exit_flag = 1;
    end

end