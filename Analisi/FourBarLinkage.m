% Articulated Quadrilateral's Closing Equations

function [obj, jacob] = FourBarLinkage(L1,L2,L3,L4)
    
    if (nargin < 4)
        error('Too few input arguments');
    end

    if (nargout < 2)
        error('Too few output arguments');
    end

    if (L1 < 0 || L2 < 0 || L3 < 0 || L4 < 0)
        error('Lengths must be positive numbers')
    end
    
    function f = funct(x)
        if (size(x) ~= 4)
            error('The input argument has wrong dimension');
        end
        f = zeros(2,1);
        f(1) = L1*cos(x(1)) + L2*cos(x(2)) + L3*cos(x(3)) + L4*cos(x(4));
        f(2) = L1*sin(x(1)) + L2*sin(x(2)) + L3*sin(x(3)) + L4*sin(x(4));
    end
    obj = @funct;

    function j = jacobian(x)
        if (size(x,1) ~= 4)
            error('The input argument has wrong dimension');
        end
        
        j = zeros(2,2);
        j(1,1) = -L2*sin(x(2));
        j(1,2) = -L3*sin(x(3));
        j(2,1) = L2*cos(x(2));
        j(2,2) = L3*cos(x(3));
    end
    jacob = @jacobian;
end

