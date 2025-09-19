
function eggxample01()
    clf;

    egg_params = struct();
    egg_params.a = 3; egg_params.b = 2; egg_params.c = .15;

    %specify the position and orientation of the egg
    x0 = 5; y0 = 5; theta = pi/6;
    %set up the axis
    figure(); hold on; axis equal; axis square;
    axis([0,10,0,10])
    
    % s has a period of 1
    s_perimeter = linspace(0, 1, 100);
    [V_vals, ~] = egg_func(s_perimeter, x0, y0, theta, egg_params);

    % plot egg
    plot(V_vals(1, :), V_vals(2, :), 'k', LineWidth=2)
    %plot(V_vals(1, :), V_vals(2, :), 'ro', LineWidth=4)

    % plot bounding box
    [xrange, yrange] = find_bounding_box(x0, y0, theta, egg_params);
    rectangle("Position", [xrange(1), yrange(1), diff(xrange), diff(yrange)])

    % s_tangent = 0.3;
    % [V_tangent, G_tangent] = egg_func(s_tangent, x0, y0, theta, egg_params);
    % 
    % plot(V_tangent(1), V_tangent(2), 'ro', MarkerFaceColor='r', MarkerSize=4)
    % plot(V_tangent(1)+[0, G_tangent(1)], V_tangent(2)+[0, G_tangent(2)], 'k', Linewidth=2)

    % check collision
    figure(); hold on; axis equal; axis square;
    axis([0,15,0,15])
    y_ground = -0.5; x_wall = 14;
    x1 = 0:2:14; y1 = 0.5*ones(8);
    plot(x1, y1)   % plot ground
    plot(14*ones(4), 0.5:2:6.5) % plot the wall

    [t_ground, t_wall] = collision_func(egg_params, y_ground, x_wall);

    % plot egg
    t_plot = min([t_ground, t_wall]);
    [x0, y0, theta] = egg_traj(t_plot);
    [V, ~] = egg_func(s_perimeter, x0, y0, theta, egg_params);

    % plot egg
    plot(V(1, :), V(2, :), 'k', LineWidth=2)
    title("collision test plot")

    % animate the egg
    %animate_egg(egg_params, y_ground, 34)

    % make the video 
    clf;
    video_example(egg_params, y_ground, 34)
end

%Short example demonstrating how to create a MATLAB animation
%In this case, a square moving along an elliptical path
%This version also store the animation in a vide.
function video_example(egg_params, y_ground, x_wall)
    %define location and filename where video will be stored
    %written a bit weird to make it fit when viewed in assignment
    mypath1 = 'C:\Users\lodio\OneDrive - Olin College of Engineering\Desktop\';
    mypath2 = 'Classes\Junior Fall\Orion Math\MechE-Math-Assignment-1';
    fname='egg_animation.avi';
    input_fname = [mypath1,mypath2,fname];
    
    %create a videowriter, which will write frames to the animation file
    writerObj = VideoWriter(input_fname);
    %must call open before writing any frames
    open(writerObj);
    
    % throw in the animate code
    % set up the plotting axis
    fig4 = figure(4); hold on; axis equal; axis square;
    axis_end = 40; %1.2*x_wall;
    axis([-1, axis_end, -1, axis_end])

    % plot boundaries
    n=6; x_vals_g = 0:x_wall/n:x_wall; 
    y_vals_g = y_ground*ones(n+1);

    x_vals_w = x_wall*ones(n+1);
    y_vals_w = y_ground:axis_end/n:axis_end;

    plot(x_vals_g, y_vals_g)    % plot ground
    plot(x_vals_w, y_vals_w)    % plot wall

    % initialize the egg plot
    [x0, y0, theta] = egg_traj(0);
    s_perim = linspace(0, 1, 100);  % s has a period of 1
    [V_vals, ~] = egg_func(s_perim, x0, y0, theta, egg_params);
    egg_plot = plot(V_vals(1, :), V_vals(2, :), 'k', LineWidth=2);

    % calc collision time
    [t_ground, t_wall] = collision_func(egg_params, y_ground, x_wall);
    t_collide = min([t_ground, t_wall]);

    hits_wall_first = 1;
    if t_collide == t_ground
        hits_wall_first = 0;
    end
        

    %iterate through time
    for t = 0:.001:t_collide
        
        % get position of egg at time t
        [x0, y0, theta] = egg_traj(t);
        
        % extract position data
        [V_vals, ~] = egg_func(s_perim, x0, y0, theta, egg_params);

        % assign pos data to vars
        x_plot = V_vals(1, :);
        y_plot = V_vals(2, :);

        % update the coordinates of the egg plot
        set(egg_plot,'xdata', x_plot,'ydata', y_plot);
        %update the actual plotting window
        drawnow;

        %capture a frame (what is currently plotted)
        current_frame = getframe(fig4);
        %write the frame to the video
        writeVideo(writerObj,current_frame);
    end

    %must call close after all frames are written
    close(writerObj);

    % % plot a dot at collision point
    % if hits_wall_first
    %     plot(x_wall, y0, "ro", MarkerSize=3)
    % else
    %     plot(x0, y_ground, "ro", MarkerSize=3)
    % end
end


function animate_egg(egg_params, y_ground, x_wall)
    
    % set up the plotting axis
    figure(); hold on; axis equal; axis square;
    axis_end = 40; %1.2*x_wall;
    axis([-1, axis_end, -1, axis_end])

    % plot boundaries
    n=6; x_vals_g = 0:x_wall/n:x_wall; 
    y_vals_g = y_ground*ones(n+1);

    x_vals_w = x_wall*ones(n+1);
    y_vals_w = y_ground:axis_end/n:axis_end;

    plot(x_vals_g, y_vals_g)    % plot ground
    plot(x_vals_w, y_vals_w)    % plot wall

    % initialize the egg plot
    [x0, y0, theta] = egg_traj(0);
    s_perim = linspace(0, 1, 100);  % s has a period of 1
    [V_vals, ~] = egg_func(s_perim, x0, y0, theta, egg_params);
    egg_plot = plot(V_vals(1, :), V_vals(2, :), 'k', LineWidth=2);

    % calc collision time
    [t_ground, t_wall] = collision_func(egg_params, y_ground, x_wall);
    t_collide = min([t_ground, t_wall]);

    hits_wall_first = 1;
    if t_collide == t_ground
        hits_wall_first = 0;
    end
        

    %iterate through time
    for t = 0:.001:t_collide
        
        % get position of egg at time t
        [x0, y0, theta] = egg_traj(t);
        
        % extract position data
        [V_vals, ~] = egg_func(s_perim, x0, y0, theta, egg_params);

        % assign pos data to vars
        x_plot = V_vals(1, :);
        y_plot = V_vals(2, :);

        % update the coordinates of the egg plot
        set(egg_plot,'xdata', x_plot,'ydata', y_plot);
        %update the actual plotting window
        drawnow;
    end

    
    % plot a dot at collision point
    if hits_wall_first
        plot(x_wall, y0, "ro", MarkerSize=3)
    else
        plot(x0, y_ground, "ro", MarkerSize=3)
    end
end


% Function that computes the collision time for a thrown egg
% INPUTS:
    % traj_fun: a function that describes the [x,y,theta] trajectory
    %           of the egg (takes time t as input)
    % egg_params: a struct describing the hyperparameters of the oval
    % y_ground: height of the ground
    % x_wall: position of the wall
% OUTPUTS:
    % t_ground: time that the egg would hit the ground
    % t_wall: time that the egg would hit the wall
function [t_ground, t_wall] = collision_func(egg_params, y_ground, x_wall)
    % use bisection method

    % find an upper bound (time t when the egg is past collision)
    t_upperx = 1; t_uppery = 1; t_lower = 0;
    xmax = xmax_at_t_wrapper(t_upperx, egg_params); 
    ymin = ymin_at_t_wrapper(t_uppery, egg_params);
    
    while  (ymin > y_ground)
        t_uppery = 2*t_uppery;
        
        ymin = ymin_at_t_wrapper(t_uppery, egg_params);
    end

    while (xmax < x_wall)
        t_upperx = 2*t_upperx;
        xmax = xmax_at_t_wrapper(t_upperx, egg_params); 
    end
    
    % establish the functions i'm finding the roots of
    xmax_over_t = @(t) xmax_at_t_wrapper(t, egg_params) - x_wall; 
    ymin_over_t = @(t) ymin_at_t_wrapper(t, egg_params) - y_ground;
    
    % establish inputs
    dxtol = 1e-14; ytol = 1e-14; max_iter = 200;
    % call solver on functions
    t_wall = bisection_solver(xmax_over_t, t_lower, t_upperx, max_iter, dxtol, ytol);
    t_ground = bisection_solver(ymin_over_t, t_lower, t_uppery, max_iter, dxtol, ytol);

end

function xmax = xmax_at_t_wrapper(t, egg_params)
    [x0, y0, theta] = egg_traj(t);
    [x_range, ~] = find_bounding_box(x0, y0, theta, egg_params);
    xmax = x_range(2);
end

function ymin = ymin_at_t_wrapper(t, egg_params)
    [x0, y0, theta] = egg_traj(t);
    [~, y_range] = find_bounding_box(x0, y0, theta, egg_params);
    ymin = y_range(1);
end


function [x_range, y_range] = find_bounding_box(x0, y0, theta, egg_params)
    % given the egg position and orientation, find roots of Gx and Gy
    % (where G = [Gx, Gy] and defines the tangent line to the egg)
    % which correspond to the points which define the bounding box

    % wrapper funcs return the Gx and Gy component, respectively, in 
    % accordance with the point defined by the inputs into egg_func
    egg_wrapper_func2_x = @(s_in) egg_wrapper_func1_x(s_in, x0, y0, theta, egg_params);
    egg_wrapper_func2_y = @(s_in) egg_wrapper_func1_y(s_in, x0, y0, theta, egg_params);

    % make containers
    x_list = []; y_list = []; 

    % generate evenly spaced guesses (s ranges 0-1)
    s_guess_list = 0:0.2:1;
    % establish params (not used with fzero, used in our solvers)
    dxtol = 1e-14; ytol = 1e-14; max_iter = 200; dfdxmin = 1e-8;

    for s_guess = s_guess_list
        % s_rootx = secant_solver(egg_wrapper_func2_x, s_guess, dxtol, ytol, max_iter, defdxmin);
        % use fzero for now
        % s_rootx = fzero(egg_wrapper_func2_x, s_guess);
        % s_rooty = fzero(egg_wrapper_func2_y, s_guess);
        s_rootx = secant_solve(egg_wrapper_func2_x, s_guess, s_guess+0.05, max_iter, dxtol, ytol);
        s_rooty = secant_solve(egg_wrapper_func2_y, s_guess, s_guess+0.05, max_iter, dxtol, ytol);

        % eval egg_func at the roots, store relevant pos param
        [V, ~] = egg_func(s_rootx, x0, y0, theta, egg_params);
        x_list(end+1) = V(1);
        [V, ~] = egg_func(s_rooty, x0, y0, theta, egg_params);
        y_list(end+1) = V(2);
    end

    % extract max and min (the bounding vals)
    x_range = [min(x_list), max(x_list)];
    y_range = [min(y_list), max(y_list)];
end

function Gx = egg_wrapper_func1_x(s, x0, y0, theta, egg_params)
    [~, G] = egg_func(s, x0, y0, theta, egg_params);
    Gx = G(1);
end

function Gy = egg_wrapper_func1_y(s, x0, y0, theta, egg_params)
    [~, G] = egg_func(s, x0, y0, theta, egg_params);
    Gy = G(2);
end


% Example parabolic trajectory
function [x0, y0, theta] = egg_traj(t)
    x0 = 7*t + 8;
    y0 = -6*t.^2 + 20*t + 6;
    theta = 5*t;
end

% This function generates the parametric curve describing an oval
% INPUTS:
    % s: the curve parametr. s is a number from 0 to 1. The curve fcn has
    %       a period of 1, so s=0.3 and s=1.3 will generate the same output
    %       s can also be a list (row vector) of numbers
    % theta: rotation of the oval. theta is a number from 0 to 2*pi.
    %        Increasing theta rotates the oval counterclockwise
    % x0: horizontal offset of the oval
    % y0: vertical offset of the oval
    % egg_params: a struct describing the hyperparameters of the oval
    %       egg_params has three variables, a,b, and c
    %       without any rotation/translation, the oval satisfies the eqn:
    %       xˆ2/aˆ2 + (yˆ2/bˆ2)*eˆ(c*x) = 1
    %       tweaking a,b,c changes the shape of the oval
% OUTPUTS:
    % V: the position of the point on the oval given the inputs
    %       If s is a single number, then V will have the form of a
    %       column vector [x_out;y_out] where (x_out, y_out) are the coords 
    %       of the point on the oval. If the input t is a list of numbers 
    %       (a row vector) i.e.: s = [s_1,...,s_N]
    %       then V will be an 2xN matrix: [x_1,...,x_N; y_1,...,y_N]
    %       where (x_i,y_i) correspond to input s_i
    % G: the gradient of V taken with respect to s
    %       If s is a single number, then G will be the col vector [dx/ds; dy/ds]
    %       If s is the list [s_1,...,s_N], then G will be the 2xN matrix:
    %       [dx_1/ds_1,...,dx_N/ds_N; dy_1/ds_1,...,dy_N/ds_N]

function [V, G] = egg_func(s,x0,y0,theta,egg_params)
    %unpack the struct
    a=egg_params.a;
    b=egg_params.b;
    c=egg_params.c;
    %compute x (without rotation or translation)
    x = a*cos(2*pi*s);
    %useful intermediate variable
    f = exp(-c*x/2);
    %compute y (without rotation or translation)
    y = b*sin(2*pi*s).*f;
    %compute the derivatives of x and y (without rotation or translation)
    dx = -2*pi*a*sin(2*pi*s);
    df = (-c/2)*f.*dx;
    dy = 2*pi*b*cos(2*pi*s).*f + b*sin(2*pi*s).*df;
    %rotation matrix corresponding to theta
    R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
    %compute position and gradient for rotated + translated oval
    V = R*[x;y]+[x0*ones(1,length(theta));y0*ones(1,length(theta))];
    G = R*[dx;dy];
end
