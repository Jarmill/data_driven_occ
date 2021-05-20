function [p] = line_range(coeff, xrange, yrange)
    %find extreme points of line inside bounding box xrange, yrange
    %line is coeff'*[x y 1] = 0 = ax + by + c
    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    
    p={};
    
    %really bad code
    
    %when x is low
    y_low = (-c-a*xrange(1))/b;
    
    if (yrange(1) <= y_low) && (y_low <= yrange(2))
        p{end+1} = [xrange(1), y_low];
    end
    
    %x is high
    y_high = (-c-a*xrange(2))/b;
    
    if (yrange(1) <= y_high) && (y_high <= yrange(2))
        p{end+1} = [xrange(2), y_high];
    end
    
    %when x is low
    if length(p) < 2
        x_low = (-c-b*yrange(1))/a;
        if (xrange(1) <= x_low) && (x_low <= xrange(2))
            p{end+1} = [x_low, yrange(1)];
        end
    end
    %x is high
    if length(p) < 2
        x_high = (-c-b*yrange(2))/a;
        if (xrange(1) <= x_high) && (x_high <= xrange(2))
            p{end+1} = [x_high, yrange(2)];
        end    
    end
    
end