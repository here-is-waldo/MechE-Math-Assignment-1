%%SECANT SOLVER
function [x_root,exit_flag,x_guess_list] = secant_solve(f_in,x0,x1,max_iter,dx_tol,y_tol)
  x_prev = x0;
  x_tprev = x1;   
  x_val = x_prev+2*dx_tol;

  [f_val] = f_in(x_val);
  % [f_prev] = f_in(x_prev);
  % [f_tprev] = f_in(x_tprev);    
    
  count = 0;

  x_guess_list = x_val;
  x_prev_guess_list = x_prev;

  check = abs((x_val-x_prev));

  while check > dx_tol && count<max_iter && abs(f_val) > y_tol
        [f_val] = f_in(x_val);
        [f_prev] = f_in(x_prev);
        [f_tprev] = f_in(x_tprev);


        x_val = (x_tprev*f_prev-x_prev*f_tprev)/(f_prev-f_tprev);
        check = abs((x_val-x_prev));
        x_tprev=x_prev;        
        x_prev=x_val;

        x_guess_list(end+1) = x_val;
        x_prev_guess_list(end+1) = x_prev; 

        count=count+1;
   end

   x_root = x_val;

   %assume success
   exit_flag = 0;
   
   
   if abs(x_val-x_prev) > dx_tol && abs(f_val)>y_tol
       exit_flag = 1;
   end

end