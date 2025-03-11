function res = max_star(a,b,c,d,e,f,g,h)

switch nargin

case 2
    res = max(a,b) + log(1+exp(-abs(a-b)));
case 3
    res = max_star(a,b);
    res = max(res,c) + log(1+exp(-abs(res-c)));
case 4
    res = max_star(a,b);
    res = max(res,c) + log(1+exp(-abs(res-c)));
    res = max(res,d) + log(1+exp(-abs(res-d)));
case 5
    res = max_star(a,b);
    res = max(res,c) + log(1+exp(-abs(res-c)));
    res = max(res,d) + log(1+exp(-abs(res-d)));
    res = max(res,e) + log(1+exp(-abs(res-e)));
case 6
      res = max_star(a,b);
    res = max(res,c) + log(1+exp(-abs(res-c)));
    res = max(res,d) + log(1+exp(-abs(res-d)));
    res = max(res,e) + log(1+exp(-abs(res-e)));
    res = max(res,f) + log(1+exp(-abs(res-f)));
case 7
     res = max_star(a,b);
    res = max(res,c) + log(1+exp(-abs(res-c)));
    res = max(res,d) + log(1+exp(-abs(res-d)));
    res = max(res,e) + log(1+exp(-abs(res-e)));
    res = max(res,f) + log(1+exp(-abs(res-f)));
    res = max(res,g) + log(1+exp(-abs(res-g)));
case 8
    res = max_star(a,b);
    res = max(res,c) + log(1+exp(-abs(res-c)));
    res = max(res,d) + log(1+exp(-abs(res-d)));
    res = max(res,e) + log(1+exp(-abs(res-e)));
    res = max(res,f) + log(1+exp(-abs(res-f)));
    res = max(res,g) + log(1+exp(-abs(res-g)));
    res = max(res,h) + log(1+exp(-abs(res-h)));
otherwise
    res = 0;
end

end


