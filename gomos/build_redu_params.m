
%[~, new_D2, ~] = svd(new_H2);


[U, new_D, ~] = svd(new_H);
new_P = U(:, 1:nps(end));

[U, lis_D, ~] = svd(lis_H);
lis_P = U(:, 1:nps(end));

[U, as_D, ~] = svd(as_H);
as_P = U(:, 1:nps(end));

[U, prlis_D, ~] = svd(prlis_H);
prlis_P = U(:, 1:nps(end));

[U, lap_D, ~] = svd(lap_H);
lap_P = U(:, 1:nps(end));

[U, laplis_D, ~] = svd(laplis_H);
laplis_P = U(:, 1:nps(end));