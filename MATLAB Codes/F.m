function [out] = F(x,y,Ps)
    int = int32(x == y);
    if int == 1
        out = (1-Ps);
    else
        out = Ps;
    end
end