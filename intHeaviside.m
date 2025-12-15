function [int_h] = intHeaviside(t)
    if t>=1
        int_h=t;
        else
        int_h = 0;
    end
end