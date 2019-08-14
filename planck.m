%% Define planck function

h = 6.626e-34; % planck constant (J*s)
c = 299792458; % speed of light in a vacuum (m/s)
kB = 1.38065e-23; % Boltzmann constant (J/K)

T = linspace(1,6000,6000)'; % temperature array
dT = 1/(T(end)-T(1));
yl = linspace(1e-7, 2e-4, 10000)'; % lambda wavelength array

B = nan(length(T), length(yl)); % preallocate 2D Planck function array 

yl1 = 14e-6;
yl2 = 16e-6;
loc = find(yl >= yl1 & yl <= yl2); % define wavelength region of interest

for i = 1 : length(T)  
    for j = 1 : length(yl)        
        B(i,j) = 2*h*c^2/yl(j)^5 / (exp( h*c / (kB * T(i) * yl(j))) -1) ;       
    end
end

%% Integrate Planck function

P = nan(length(T),1); % total power emitted
p = nan(length(T),1); % power emitted through small window

for i = 1 : length(T)
    P(i) = sum(B(i,:))*dT;
    p(i) = sum(B(i,loc))*dT;
end

frac = p./P;

%% Plot Planck function

temp = 200;
numTemps = 10;
leg = nan(16,1);

figure(1)
clf
hold on
for i = 1 : numTemps
    temp = temp + 10;
    leg(i) = temp;
    plot(yl*1e6, B(temp,:), 'linewidth', 2)
end
temp = 200;
for i = 1 : numTemps
    leg(i+numTemps) = "";
    temp = temp + 10;
    plot(yl(loc)*1e6, B(temp,loc), 'w-','linewidth', 3);
end
plot( [yl1*1e6 yl1*1e6], [ 0 max(B(temp,loc))] , 'k-')
plot( [yl2*1e6 yl2*1e6], [ 0 max(B(temp,loc))] , 'k-')
xlim([0 80])
ylabel(('Radiant Intensity (W \cdot sr^{-1} \cdot m^{-3})'))
xlabel('Wavelength ({\mu}m)')
set(gca, 'fontsize', 18)
legend(strcat(num2str(leg(1:numTemps)), " K"))
hold off

%% Plot fractional emission of total power through window 

figure(2)
clf
hold on
plot(T, frac*100, 'linewidth', 2)
ylabel(['% total power emitted between [', num2str(yl1*1e6), '{\mu}m , ', num2str(yl2*1e6), '{\mu}m ]'])
xlabel('Temperature (K)')
set(gca, 'fontsize', 18)
xlim([200, 300])
hold off









