function ev = get_edges(t, x)
% x is logical or 0/1 numeric
x = logical(x(:));
t = double(t(:));

dx = diff(x);
ev.rise_idx = find(dx == 1) + 1;
ev.fall_idx = find(dx == -1) + 1;

ev.rise_t = t(ev.rise_idx);
ev.fall_t = t(ev.fall_idx);
end
