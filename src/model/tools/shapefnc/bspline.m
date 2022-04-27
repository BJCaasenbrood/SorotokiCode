function y = bspline(x, j, d, t, cur_d )

    function [ f ] = bspline_impl(x, j, d, t, cur_d)
        f = zeros(size(x));
        if cur_d == 0
            if (j + 1) >= numel(t) - d
                f(x >= t(j)) = 1;
            else
                f(x >= t(j) & x < t(j + 1)) = 1;
            end
        else
            if t(j + cur_d) ~= t(j)
                leftc = (x - t(j)) / (t(j + cur_d) - t(j));
                f = f + leftc.*bspline_impl(x, j, d, t, cur_d - 1);
            end
            
            if t(j + cur_d + 1) ~= t(j + 1)
                rightc = (t(j + cur_d + 1) - x) / (t(j + cur_d + 1) - t(j + 1));
                f = f + rightc.*bspline_impl(x, j + 1, d, t, cur_d - 1);
            end
        end
    end

    y = bspline_impl(x, j, d, t, cur_d);
end