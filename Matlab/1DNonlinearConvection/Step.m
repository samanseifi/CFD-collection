function [ u_init ] = Step( x, nx)
    
    for i = 1:nx
        if (x(i) >= 0) && (x(i) <= 1)
            u_init(i) = 1.0;
        else
            u_init(i) = 0;
        end
    end

end

