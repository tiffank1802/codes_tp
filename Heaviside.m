function [h] = Heaviside(t)
    if t>=1
        h=1;
    else
        h = 0;
    end
end