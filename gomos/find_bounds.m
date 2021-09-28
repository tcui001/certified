function j = find_bounds(V, ref_H, tol, max)

VH = ref_H*V;
j = max;
for i = 1:max
    bound = trace(ref_H - VH(:,1:i)*V(:,1:i)' );
    if bound < tol
        j = i;
        break;
    end
end

end