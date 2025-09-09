function Assignment1Secant()
    x_range = linspace(-10,38,200);

    [y_vals] = test_func(x_range);

     x0 = -2;
     max_iter = 200;
     dx_tol = 1e-7;
     y_tol = 1e-7;

    [~,x_root1,exit_flag] = secant_solve(@test_func, x0, max_iter, dx_tol,y_tol);
    [f_root1] = test_func(x_root1);
    
    x_root_base = x_root1;

    x0=35;
    [~,x_root2,exit_flag] = secant_solve(@test_func, x0, max_iter, dx_tol,y_tol);
    [f_root2] = test_func(x_root2);

    figure(1);
    hold on;
    plot(x_range,y_vals,'k')
    plot(x_range,x_range*0,'r--')
    plot(x_root1,f_root1,"bo",'MarkerFaceColor','b','MarkerSize',3);
    plot(x_root2,f_root2,"bo",'MarkerFaceColor','b','MarkerSize',3);

    all_iteration_list =[];
    all_error_list = [];

    all_error_current_list = [];
    all_error_next_list = [];

     for n = 1:100
         x0 = x_root_base + (2*rand()-.5);
         [~,~,x_guess_list] = secant_solve(@test_func, x0, max_iter, dx_tol,y_tol);
         iteration_list = 1:length(x_guess_list);
         error_list = abs(x_guess_list-x_root_base);
         all_iteration_list = [all_iteration_list, iteration_list];
         all_error_list = [all_error_list,error_list];

         current_error = error_list(1:end-1);
         next_error = error_list(2:end);

         all_error_current_list = [all_error_current_list,current_error];
         all_error_next_list = [all_error_next_list,next_error];
     end
    
    figure(2);
    plot(abs(x_guess_list-x_root1));

    figure(3);
    loglog(all_error_current_list,all_error_next_list,"or")
end

%%SECANT SOLVER
function [x_root,exit_flag,x_guess_list] = secant_solve(f_in,x0,max_iter,dx_tol,y_tol)
  x_prev = x0-2*dx_tol;
  %x_tprev = x_prev-2*dx_tol;   
  x_val = x0;
  [f_val] = f_in(x_val);
  [f_prev] = f_in(x_prev);
  %[f_tprev] = f_in(x_tprev);
    
  count = 0;
  x_guess_list = x_val;
  x_prev_guess_list = x_prev;

  while abs((x_val-x_prev)) > dx_tol && count<max_iter && abs(f_val) > y_tol && abs(f_prev) > y_tol
        %x_tprev = x_prev;
        x_prev = x_val;
        x_val = x_val - f_val*((x_val-x_prev)/(f_val-f_prev));

        [f_val] = f_in(x_val);
        [f_prev] = f_in(x_prev);
        %[f_tprev] = f_in(x_tprev);
        count=count+1;
        x_guess_list(end+1) = x_val;
        x_prev_guess_list(end+1) = x_prev; 
   end

    x_root = x_val;

   %assume success
   exit_flag = 0;

   if abs(x_val-x_prev) > dx_tol && abs(f_val)>y_tol
       exit_flag = 1;
   end

end

%%IMPORT FUNCTION
function [f_val] = test_func(x_range)
    f_val = (x_range.^3)/100 - (x_range.^2)/8 + 2*x_range + 6*sin(x_range/2+6) -.7 - exp(x_range/6);
end