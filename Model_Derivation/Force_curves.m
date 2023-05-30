% Define the parameters
slip_ratio = 0.2;
ca = 150000;
cs = 150000;
mu = 1;
Fz = 5000;

% Create an array of slip angles
slip_angles = linspace(-0.5, 0.5, 100);

% Preallocate an array to hold the results
Fy_bars = zeros(size(slip_angles));

% Loop over the slip angles and compute the corresponding Fy_bar values
for i = 1:length(slip_angles)
    slip_angle = slip_angles(i);
    [Fy_bar, ca, ~] = DugoffLinearizer(slip_angle, slip_ratio, ca, cs, mu, Fz);
    Fy_bars(i) = ca;
end

% Plot the results
plot(slip_angles, Fy_bars);
xlabel('Slip angle');
ylabel('Fy_bar');
title('Dugoff tire model');
