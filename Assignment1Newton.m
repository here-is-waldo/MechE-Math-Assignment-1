function Assignment1Newton()
    x_range = linspace(-10,38,200);
    input_list = setGlobal(x_range);

    [y_vals,dfdx] = test_func(x_range);

     x0 = -2;
     max_iter = 200;
     dx_tol = 1e-14;
     y_tol = 1e-14;

    [x_root1,exit_flag] = newton_solve(@test_func, x0, max_iter, dx_tol,y_tol);
    [f_root1] = test_func(x_root1);
    
    x_root_base = x_root1;

    x0=35;
    [x_root2,exit_flag] = newton_solve(@test_func, x0, max_iter, dx_tol,y_tol);
    [f_root2] = test_func(x_root2);

    figure(1);
    hold on;
    plot(x_range,y_vals,'k')
    plot(x_range,x_range*0,'r--')
    plot(x_root1,f_root1,"bo",'MarkerFaceColor','b','MarkerSize',3);
    plot(x_root2,f_root2,"bo",'MarkerFaceColor','b','MarkerSize',3);

    all_iteration_list =[];
    all_error_list = [];

    x_current_list = [];
    x_next_list =[];

    index_list = [];

    all_error_current_list = [];
    all_error_next_list = [];

    all_x_current_list = [];
    all_x_next_list = [];

     for n = 1:(max_iter/2)
         x0 = x_root_base + (2*rand()-.5);
         [~,~,x_guess_list] = newton_solve(@test_func, x0, max_iter, dx_tol,y_tol);
         iteration_list = 1:length(x_guess_list);

         %input_list =[];

         error_list = abs(x_guess_list-x_root_base);
         all_iteration_list = [all_iteration_list, iteration_list];
         all_error_list = [all_error_list,error_list];

         x_current_list = [x_current_list,input_list(1:end-1)];
         x_next_list = [x_next_list,input_list(2:end)];

         current_error = error_list(1:end-1);
         next_error = error_list(2:end);

         all_error_current_list = [all_error_current_list,current_error];
         all_error_next_list = [all_error_next_list,next_error];

         all_x_current_list = [all_x_current_list,input_list(1:end-1)];
         all_x_next_list = [all_x_next_list,input_list(2:end)];

         index_list = [index_list,1:length(input_list)-1];
     end


    figure(2);
    
    plot(abs(x_guess_list-x_root1));

    x_regression = []; % e_n
    y_regression = []; % e_{n+1}
    %iterate through the collected data
    error_list0 = all_error_current_list;
    error_list1 = all_error_next_list; 

for n=1:length(index_list)
    if error_list0(n)>1e-14 && error_list0(n)<1e-2 && error_list1(n)>1e-14 && error_list1(n)<1e-2 %&& index_list(n)>2
        x_regression(end+1) = error_list0(n);
        y_regression(end+1) = error_list1(n);
        disp('added')
    end
end
    [p,k] = generate_error_fit(x_regression,y_regression);
    
    fit_line_x = 10.^[-16:.01:1];
    fit_line_y = k*fit_line_x.^p;

    figure(3);
    loglog(all_error_current_list,all_error_next_list,'ro','markerfacecolor',[1,0,0],'markersize',5)
    hold on;
    loglog(x_regression,y_regression,'bo','markerfacecolor',[0,0,1],'markersize',5)
    loglog(fit_line_x,fit_line_y,'k-','linewidth',2)
    hold off;
    xlim([10^-16 10^1])
    ylim([10^-16 10^1])
end

%%NEWTON SOLVER
function [x_root,exit_flag,x_guess_list] = newton_solve(f_in,x0,max_iter,dx_tol,y_tol)
    x_prev = x0-2*dx_tol;    
    x_val = x0;
    [f_val,dfdx_val] = f_in(x_val);

    count = 0;
    x_guess_list = x_val;

  while abs(dfdx_val) > 1e-7 && abs((x_val-x_prev)) > dx_tol && count<max_iter && abs(f_val) > y_tol
        x_prev = x_val;
        x_val = x_val - f_val/dfdx_val;

        [f_val, dfdx_val] = f_in(x_val);
        count=count+1;
        x_guess_list(end+1) = x_val;
   end

    x_root = x_val;

   %assume success
   exit_flag = 0;

   if abs(x_val-x_prev) > dx_tol && abs(f_val)>y_tol
       exit_flag = 1;
   end

end

%%IMPORT FUNCTION
function [f_val,dfdx] = test_func(x_range)
    f_val = (x_range.^3)/100 - (x_range.^2)/8 + 2*x_range + 6*sin(x_range/2+6) -.7 - exp(x_range/6);
    dfdx = 3*(x_range.^2)/100 - 2*x_range/8 + 2 +(6/2)*cos(x_range/2+6) - exp(x_range/6)/6;
end

function input_list = setGlobal(x_range)
    global input_list
    input_list(:,end+1) = x_range;
end

function [p,k] = generate_error_fit(x_regression,y_regression)

    Y = log(y_regression)';
    X1 = log(x_regression)';
    X2 = ones(length(X1),1);
    MEGAX = [X1,X2];
%run the regression
    coeff_vec = regress(Y,MEGAX);
%pull out the coefficients from the fit
    p = coeff_vec(1);
    k = exp(coeff_vec(2));
end