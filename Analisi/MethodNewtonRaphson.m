% Newton - Rhapson's Method

function [theta_est, i] = MethodNewtonRaphson(theta_est, epsilon, max_iterations, f, j)

    if (nargin < 4)
        error('Too few input arguments');
    end
    
    if (nargout < 2)
        error('Too few output arguments');
    end

    if (epsilon < 0 || max_iterations < 0)
        error('Lengths must be positive numbers')
    end

    i = 0;

    while 1
        f_i = f(theta_est);
        j_i = j(theta_est);
        if (rank(j_i) ~= size(j_i,1))
            error("The jacobian can not be inverted");
        end
    
        delta_theta = -j_i\f_i;
        if (norm(delta_theta, 2) <= epsilon || i >= max_iterations)
            break;
        end
        theta_est(2) = theta_est(2) + delta_theta(1);
        theta_est(3) = theta_est(3) + delta_theta(2);

        i = i + 1;
    end
end

