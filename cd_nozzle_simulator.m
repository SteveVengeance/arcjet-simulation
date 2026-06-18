%% CD Nozzle Flow Simulator with Electric Arc Heating
% Author : Arif Hossain - 10/10/2015
% Modified: Multi-gas support, 1D/2D flow properties, back pressure analysis,
%           electric arc heating, forward simulation and inverse design modes

clear all; close all

%% =========================================================
%  MODE SELECTION
%  =========================================================
fprintf('=========================================\n')
fprintf('  CD Nozzle Simulator with Arc Heating\n')
fprintf('=========================================\n')
fprintf('Select operating mode:\n')
fprintf('  1. Forward Simulation  — specify nozzle design and simulate\n')
fprintf('  2. Inverse Design      — specify target exit conditions, find optimal parameters\n')
sim_mode = input('Enter mode (1 or 2): ');

%% =========================================================
%  GAS SELECTION (both modes)
%  =========================================================
fprintf('\nSelect gas:\n')
fprintf('  1. Air\n')
fprintf('  2. Nitrogen (N2)\n')
fprintf('  3. Oxygen (O2)\n')
fprintf('  4. Helium (He)\n')
fprintf('  5. Argon (Ar)\n')
fprintf('  6. Hydrogen (H2)\n')
fprintf('  7. Carbon Dioxide (CO2)\n')
gas_choice = input('Enter gas number: ');

gas_names = {'Air', 'Nitrogen (N2)', 'Oxygen (O2)', 'Helium (He)', ...
             'Argon (Ar)', 'Hydrogen (H2)', 'Carbon Dioxide (CO2)'};
gama_list = [1.400, 1.400, 1.400, 1.667, 1.667, 1.405, 1.289];
R_list    = [287.0, 296.8, 259.8, 2077.0, 208.1, 4124.0, 188.9]; % J/(kg*K)

gama     = gama_list(gas_choice);
R        = R_list(gas_choice);
gas_name = gas_names{gas_choice};
Cp       = gama * R / (gama - 1);
B        = (gama-1)/2;
fprintf('Selected: %s  (gamma = %.3f,  R = %.1f J/kg/K)\n\n', gas_name, gama, R)

%% =========================================================
%  MODE 1: FORWARD SIMULATION — user specifies all design parameters
%  =========================================================
if sim_mode == 1

    fprintf('--- Nozzle Parameters ---\n')
    design    = input('Design Mach number (0<M<5):  ');
    L         = input('Nozzle length (mm):  ');
    d_inlet   = input('Throat diameter (mm):  ');
    T0        = input('Stagnation temperature T0 (degC):  ') + 273.15;
    P0        = input('Stagnation pressure P0 (kPa):  ') * 1000;
    P_back    = input('Back pressure P_back (kPa):  ')   * 1000;

    fprintf('\n--- Electric Arc Configuration ---\n')
    V_arc     = input('Arc voltage (V):  ');
    I_arc     = input('Arc current (A):  ');
    x_cathode = input('Cathode center position (mm):  ');
    w_cathode = input('Cathode width (mm):  ');
    x_anode   = input('Anode center position (mm):  ');
    w_anode   = input('Anode width (mm):  ');

%% =========================================================
%  MODE 2: INVERSE DESIGN — user specifies targets, script finds parameters
%  =========================================================
elseif sim_mode == 2

    fprintf('--- Target Exit Conditions ---\n')
    V_exit_target   = input('Target exit velocity (m/s):  ');
    T_exit_target_C = input('Target exit temperature (degC):  ');
    P_back          = input('Ambient / back pressure P_back (kPa):  ') * 1000;

    fprintf('\n--- Gas Supply Conditions ---\n')
    T_supply_C = input('Gas supply temperature (degC):  ');
    T_supply_K = T_supply_C + 273.15;

    fprintf('\n--- Nozzle Geometry ---\n')
    d_inlet   = input('Throat diameter (mm):  ');
    L         = input('Nozzle length (mm):  ');

    fprintf('\n--- Arc Electrode Geometry ---\n')
    x_cathode = input('Cathode center position (mm):  ');
    w_cathode = input('Cathode width (mm):  ');
    x_anode   = input('Anode center position (mm):  ');
    w_anode   = input('Anode width (mm):  ');

    % --- Back-calculate design parameters from targets ---
    T_exit_K = T_exit_target_C + 273.15;
    a_exit   = sqrt(gama * R * T_exit_K);       % speed of sound at exit
    design   = V_exit_target / a_exit;           % required design Mach number

    % Stagnation temperature required after arc heating
    T0_required = T_exit_K * (1 + B * design^2);

    % Supply pressure required for perfectly expanded nozzle (P_exit = P_back)
    P0       = P_back * (1 + B * design^2)^(gama/(gama-1));
    T0       = T_supply_K;                       % supply stagnation (before arc)

    % Arc heating required
    dT0_needed = T0_required - T_supply_K;

    fprintf('\n========= Inverse Design Results =========\n')
    fprintf('Required design Mach number  : %.4f\n', design)
    fprintf('Required stagnation temp     : %.1f degC\n', T0_required - 273.15)
    fprintf('Required supply pressure P0  : %.2f kPa\n', P0/1000)

    if design < 0 || ~isreal(design)
        error('Target velocity and temperature are not physically consistent.')
    end
    if design > 5
        warning('Required Mach number (%.2f) exceeds simulation range (M=5). Results may be extrapolated.', design)
    end
    if design < 1
        fprintf('Note: design Mach < 1 — subsonic nozzle design\n')
    end

    if dT0_needed <= 0
        fprintf('Arc heating required         : NONE\n')
        fprintf('Supply temperature is already sufficient (excess: %.1f degC)\n', -dT0_needed)
        V_arc = 0;
        I_arc = 0;
    else
        fprintf('Arc temperature rise needed  : %.1f degC\n', dT0_needed)

        % Compute mass flow rate at throat using supply conditions
        T_star_inv   = T0 * 2/(gama+1);
        P_star_inv   = P0 * (2/(gama+1))^(gama/(gama-1));
        rho_star_inv = P_star_inv / (R * T_star_inv);
        V_star_inv   = sqrt(gama * R * T_star_inv);
        A_star_m2    = pi * (d_inlet/2 * 1e-3)^2;
        m_dot_inv    = rho_star_inv * V_star_inv * A_star_m2;

        P_arc_required = m_dot_inv * Cp * dT0_needed;
        fprintf('Required arc power           : %.1f W  (%.2f kW)\n', P_arc_required, P_arc_required/1000)
        fprintf('Mass flow rate               : %.5f kg/s\n', m_dot_inv)

        fprintf('\n--- Arc Power Setting ---\n')
        fprintf('  1. Specify voltage, compute current\n')
        fprintf('  2. Specify current, compute voltage\n')
        arc_input_mode = input('Select (1 or 2): ');

        if arc_input_mode == 1
            V_arc = input('Arc voltage (V):  ');
            I_arc = P_arc_required / V_arc;
            fprintf('Required arc current : %.2f A\n', I_arc)
        else
            I_arc = input('Arc current (A):  ');
            V_arc = P_arc_required / I_arc;
            fprintf('Required arc voltage : %.2f V\n', V_arc)
        end
    end

    fprintf('==========================================\n\n')
    fprintf('Running simulation with computed parameters...\n\n')

else
    error('Invalid mode selection. Please enter 1 or 2.')
end

%% =========================================================
%  AREA-MACH RELATION  (common to both modes)
%  =========================================================
A_coeff = 2/(gama+1);
C       = (gama+1)/(gama-1);
M       = (0.01:0.001:5);

AR = sqrt((1./M.^2).*(A_coeff*(1+B*M.^2)).^C);
i  = find(M>=1,      1);
k  = find(M>=design, 1);

figure(1)
if design < 1
    plot(AR(1:i), M(1:i), 'linewidth', 2); hold on
    scatter(AR(k), M(k), 'r', 'linewidth', 2)
    xlabel('A/A^*','fontsize',14); ylabel('M','fontsize',14)
    set(gca,'fontsize',13)
    title(['Area Mach Relation (Subsonic design) - ', gas_name])
    grid on; legend('Subsonic branch','Design point')
    xlim([0,10]); ylim([0,1])
else
    plot(AR(i:end), M(i:end), 'linewidth', 2); hold on
    scatter(AR(k), M(k), 'r', 'linewidth', 2)
    xlabel('A/A^*','fontsize',14); ylabel('M','fontsize',14)
    set(gca,'fontsize',13)
    title(['Area Mach Relation (Supersonic design) - ', gas_name])
    grid on; legend('Supersonic branch','Design point','location','southeast')
end

%% =========================================================
%  NOZZLE WALL PROFILE (5th-order polynomial)
%  Boundary conditions:
%    y(0)  = d_inlet,  y(L)  = d_outlet
%    y'(0) = y''(0)  = y'(L) = y''(L) = 0
%  =========================================================
j        = find(M>=design, 1);
d_outlet = sqrt(AR(j)) * d_inlet;
x        = linspace(0, L, 1000);

a_mat = [L^5 L^4 L^3; 5*L^4 4*L^3 3*L^2; 20*L^3 12*L^2 6*L];
c_mat = [d_outlet; 0; 0];
b_mat = a_mat\c_mat;
y = b_mat(1).*x.^5 + b_mat(2).*x.^4 + b_mat(3).*x.^3 + d_inlet;

%% =========================================================
%  ARC ZONE GEOMETRY
%  =========================================================
P_arc = V_arc * I_arc;

x_cat_start = max(0, x_cathode - w_cathode/2);
x_cat_end   = min(L, x_cathode + w_cathode/2);
x_an_start  = max(0, x_anode   - w_anode/2);
x_an_end    = min(L, x_anode   + w_anode/2);
x_arc_start = min(x_cat_start, x_an_start);
x_arc_end   = max(x_cat_end,   x_an_end);

%% =========================================================
%  MASS FLOW RATE AND ARC ENERGY BALANCE
%  =========================================================
T_star    = T0 * 2/(gama+1);
P_star    = P0 * (2/(gama+1))^(gama/(gama-1));
rho_star  = P_star / (R * T_star);
V_star    = sqrt(gama * R * T_star);
A_star_m2 = pi * (d_inlet/2 * 1e-3)^2;
m_dot     = rho_star * V_star * A_star_m2;
dT0_arc   = P_arc / (m_dot * Cp);

fprintf('========= Arc & Flow Summary =========\n')
fprintf('Arc power              : %.1f W  (%.2f kW)\n', P_arc, P_arc/1000)
fprintf('Mass flow rate         : %.5f kg/s\n', m_dot)
fprintf('Stagnation temp rise   : %.1f degC\n', dT0_arc)
fprintf('T0 before / after arc  : %.1f degC  /  %.1f degC\n', T0-273.15, T0+dT0_arc-273.15)
if sim_mode == 2
    fprintf('--- Verification of target ---\n')
    T0_after  = T0 + dT0_arc;
    T_exit_check = T0_after / (1 + B*design^2);
    V_exit_check = design * sqrt(gama * R * T_exit_check);
    fprintf('Achieved exit temperature : %.1f degC  (target: %.1f degC)\n', ...
        T_exit_check-273.15, T_exit_target_C)
    fprintf('Achieved exit velocity    : %.1f m/s  (target: %.1f m/s)\n', ...
        V_exit_check, V_exit_target)
end
fprintf('======================================\n\n')

%% =========================================================
%  STAGNATION TEMPERATURE PROFILE T0(x) WITH ARC HEATING
%  =========================================================
T0_x = T0 * ones(size(x));
arc_mask = (x >= x_arc_start) & (x <= x_arc_end);
if any(arc_mask)
    arc_ramp            = (x(arc_mask) - x_arc_start) / (x_arc_end - x_arc_start);
    T0_x(arc_mask)      = T0 + dT0_arc .* arc_ramp;
end
T0_x(x > x_arc_end) = T0 + dT0_arc;

%% =========================================================
%  1D ISENTROPIC FLOW PROPERTIES ALONG NOZZLE
%  =========================================================
AR_x = (y / d_inlet).^2;

if design >= 1
    M_x = interp1(AR(i:end), M(i:end), AR_x, 'linear', 'extrap');
else
    M_x = interp1(fliplr(AR(1:i)), fliplr(M(1:i)), AR_x, 'linear', 'extrap');
end

% Without arc
T_x     = T0   ./ (1 + B.*M_x.^2);
P_x     = P0   .* (1 + B.*M_x.^2).^(-gama/(gama-1));
rho_x   = P_x  ./ (R .* T_x);
V_x     = M_x  .* sqrt(gama .* R .* T_x);
T_x_C   = T_x  - 273.15;

% With arc
T_x_arc   = T0_x ./ (1 + B.*M_x.^2);
P_x_arc   = P0   .* (1 + B.*M_x.^2).^(-gama/(gama-1));
rho_x_arc = P_x_arc ./ (R .* T_x_arc);
V_x_arc   = M_x .* sqrt(gama .* R .* T_x_arc);
T_x_arc_C = T_x_arc - 273.15;

%% =========================================================
%  BACK PRESSURE ANALYSIS
%  =========================================================
shock_in_nozzle = false;
x_shock         = NaN;

if design >= 1
    P_exit_design = P_x(end);
    M_design      = M_x(end);
    P_NS_exit     = P_exit_design * (2*gama*M_design^2 - (gama-1)) / (gama+1);

    fprintf('======= Nozzle Operating Condition =======\n')
    fprintf('Isentropic exit pressure  : %8.2f kPa\n', P_exit_design/1000)
    fprintf('Normal shock at exit (P2) : %8.2f kPa\n', P_NS_exit/1000)
    fprintf('Back pressure (P_back)    : %8.2f kPa\n', P_back/1000)
    fprintf('------------------------------------------\n')

    if P_back < P_exit_design
        fprintf('Regime: UNDEREXPANDED\n')
        fprintf('        Expansion fans form outside nozzle\n')
    elseif abs(P_back - P_exit_design)/P0 < 0.005
        fprintf('Regime: PERFECTLY EXPANDED (design condition)\n')
    elseif P_back <= P_NS_exit
        fprintf('Regime: OVEREXPANDED\n')
        if abs(P_back - P_NS_exit)/P0 < 0.01
            fprintf('        Normal shock at exit plane\n')
        else
            fprintf('        Oblique shocks outside nozzle\n')
        end
    else
        fprintf('Regime: NORMAL SHOCK INSIDE NOZZLE\n')
        shock_in_nozzle = true;

        AR_exit          = AR_x(end);
        M_exit_sub       = interp1(fliplr(AR(1:i)), fliplr(M(1:i)), AR_exit, 'linear', 'extrap');
        P_exit_sub_ratio = (1 + B*M_exit_sub^2)^(-gama/(gama-1));
        P02_P0_req       = (P_back/P0) / P_exit_sub_ratio;

        M1     = M_x;
        M2_ns  = sqrt(((gama-1).*M1.^2 + 2) ./ (2*gama.*M1.^2 - (gama-1)));
        P2_P1  = (2*gama.*M1.^2 - (gama-1)) ./ (gama+1);
        P1_P0  = (1 + B.*M1.^2).^(-gama/(gama-1));
        P02_P2 = (1 + B.*M2_ns.^2).^(gama/(gama-1));
        P02_P0 = P02_P2 .* P2_P1 .* P1_P0;

        idx_shock = find(P02_P0 <= P02_P0_req, 1);
        if ~isempty(idx_shock)
            x_shock = x(idx_shock);
            M_shock = M_x(idx_shock);
            fprintf('        Shock location : x = %.1f mm  (M = %.3f)\n', x_shock, M_shock)
        end
    end
    fprintf('==========================================\n\n')
end

%% =========================================================
%  FIGURE 2: NOZZLE PROFILE WITH ARC ELECTRODES
%  =========================================================
figure(2)

fill([x, fliplr(x)], [y, fliplr(-y)], [0.75 0.88 1.0], ...
    'EdgeColor','b','linewidth',1.5,'DisplayName','Nozzle wall'); hold on

% Arc zone highlight
x_arc_vec   = x(x >= x_arc_start & x <= x_arc_end);
y_arc_upper = interp1(x, y, x_arc_vec);
if ~isempty(x_arc_vec)
    fill([x_arc_vec, fliplr(x_arc_vec)], [y_arc_upper, fliplr(-y_arc_upper)], ...
        [1.0 0.85 0.2], 'FaceAlpha', 0.35, 'EdgeColor', 'none', 'DisplayName', 'Arc zone')
end

t_elec = 0.10 * max(y);

% Cathode (red)
x_cat_vec = linspace(x_cat_start, x_cat_end, 60);
y_cat_vec = interp1(x, y, x_cat_vec);
fill([x_cat_vec, fliplr(x_cat_vec)], ...
     [y_cat_vec, fliplr(y_cat_vec + t_elec)], ...
     'r', 'EdgeColor','k','linewidth',0.5,'DisplayName','Cathode (-)')
fill([x_cat_vec, fliplr(x_cat_vec)], ...
     [-y_cat_vec, fliplr(-y_cat_vec - t_elec)], ...
     'r', 'EdgeColor','k','linewidth',0.5,'HandleVisibility','off')

% Anode (gold)
x_an_vec = linspace(x_an_start, x_an_end, 60);
y_an_vec = interp1(x, y, x_an_vec);
fill([x_an_vec, fliplr(x_an_vec)], ...
     [y_an_vec, fliplr(y_an_vec + t_elec)], ...
     [1 0.72 0], 'EdgeColor','k','linewidth',0.5,'DisplayName','Anode (+)')
fill([x_an_vec, fliplr(x_an_vec)], ...
     [-y_an_vec, fliplr(-y_an_vec - t_elec)], ...
     [1 0.72 0], 'EdgeColor','k','linewidth',0.5,'HandleVisibility','off')

y_cat_label = interp1(x, y, x_cathode) + t_elec * 2.5;
y_an_label  = interp1(x, y, x_anode)  + t_elec * 2.5;
text(x_cathode, y_cat_label, sprintf('(-)\n%.0f V', V_arc), ...
    'Color','r','FontSize',10,'HorizontalAlignment','center','FontWeight','bold')
text(x_anode, y_an_label, sprintf('(+)\n%.0f A', I_arc), ...
    'Color',[0.6 0.4 0],'FontSize',10,'HorizontalAlignment','center','FontWeight','bold')
text((x_arc_start+x_arc_end)/2, 0, sprintf('%.1f kW', P_arc/1000), ...
    'Color',[0.5 0.35 0],'FontSize',10,'HorizontalAlignment','center','FontWeight','bold')

if sim_mode == 2
    text(L*0.98, y(end)*1.3, ...
        sprintf('V_{exit}=%.0f m/s\nT_{exit}=%.0f^{\\circ}C', V_exit_target, T_exit_target_C), ...
        'FontSize', 9, 'HorizontalAlignment', 'right', 'Color', [0 0.45 0.7])
end

xlim([0, L]); ylim([-y(end)*1.6, y(end)*1.6])
xlabel('x (mm)','fontsize',14); ylabel('Radius (mm)','fontsize',14)
set(gca,'fontsize',13)
if sim_mode == 1
    title(['Nozzle Profile with Arc Electrodes - ', gas_name])
else
    title(sprintf('Nozzle Profile — Inverse Design (M=%.2f) — %s', design, gas_name))
end
legend('show','location','northwest'); grid on

%% =========================================================
%  FIGURE 3: 1D FLOW PROPERTIES — NO ARC vs WITH ARC
%  =========================================================
figure(3)

% Safe ylim helper (handles negative Celsius values)
safe_ylim = @(v) [min(v) - 0.05*range(v), max(v) + 0.05*range(v)];

ax_ylims = {safe_ylim([T_x_C,   T_x_arc_C  ]), ...
            safe_ylim([V_x,     V_x_arc    ]), ...
            safe_ylim([rho_x,   rho_x_arc  ]), ...
            [0, max(P_x/1000)*1.1]};

prop_no_arc = {T_x_C,     V_x,     rho_x,     P_x/1000};
prop_arc    = {T_x_arc_C, V_x_arc, rho_x_arc, P_x_arc/1000};
ylabels     = {'T (degC)', 'V (m/s)', '\rho (kg/m^3)', 'P (kPa)'};
titls_f3    = {'Static Temperature', 'Flow Velocity', 'Gas Density', 'Static Pressure'};
colors      = {'r','b','g','k'};

for s = 1:4
    subplot(4,1,s)

    yl = ax_ylims{s};
    patch([x_arc_start x_arc_end x_arc_end x_arc_start], ...
          [yl(1) yl(1) yl(2) yl(2)], ...
          [1.0 0.92 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
          'DisplayName', 'Arc zone'); hold on

    plot(x, prop_no_arc{s}, [colors{s},'--'], 'linewidth', 1.5, 'DisplayName', 'No arc')
    plot(x, prop_arc{s},     colors{s},        'linewidth', 2,   'DisplayName', 'With arc')

    % Target markers for inverse mode
    if sim_mode == 2
        if s == 1
            yline(T_exit_target_C, 'k:', 'linewidth', 1.5, 'DisplayName', 'Target T_{exit}')
        elseif s == 2
            yline(V_exit_target, 'k:', 'linewidth', 1.5, 'DisplayName', 'Target V_{exit}')
        end
    end

    if s == 4
        yline(P_back/1000, 'r--', 'linewidth', 1.5, 'DisplayName', 'P_{back}')
        if shock_in_nozzle && ~isnan(x_shock)
            xline(x_shock, 'm--', 'linewidth', 1.5, 'DisplayName', 'Shock location')
        end
    end

    ylim(yl)
    xlabel('x (mm)','fontsize',12); ylabel(ylabels{s},'fontsize',12)
    title(titls_f3{s},'fontsize',13); grid on
    legend('location','best','fontsize',9)
end

if sim_mode == 1
    sgtitle(sprintf('1D Flow Properties  |  Arc: %.0f V x %.0f A = %.1f kW  |  %s', ...
        V_arc, I_arc, P_arc/1000, gas_name), 'fontsize', 13)
else
    sgtitle(sprintf('1D Flow Properties  |  Inverse Design: %.0f m/s, %.0f degC  |  %s', ...
        V_exit_target, T_exit_target_C, gas_name), 'fontsize', 13)
end

%% =========================================================
%  FIGURE 4: 2D FLOW PROPERTY DISTRIBUTION INSIDE NOZZLE
%  =========================================================
figure(4)

nx_g = 600;
ny_g = 250;
x_g  = linspace(0, L, nx_g);
y_g  = linspace(-max(y)*1.15, max(y)*1.15, ny_g);
[X_g, Y_g] = meshgrid(x_g, y_g);

T_1d   = interp1(x, T_x_arc_C, x_g, 'linear');
V_1d   = interp1(x, V_x_arc,   x_g, 'linear');
rho_1d = interp1(x, rho_x_arc, x_g, 'linear');
y_wall = interp1(x, y,         x_g, 'linear');

T_2d   = repmat(T_1d,   ny_g, 1);
V_2d   = repmat(V_1d,   ny_g, 1);
rho_2d = repmat(rho_1d, ny_g, 1);

Y_wall_2d = repmat(y_wall, ny_g, 1);
outside   = abs(Y_g) > Y_wall_2d;
T_2d(outside)   = NaN;
V_2d(outside)   = NaN;
rho_2d(outside) = NaN;

x_cat_v = linspace(x_cat_start, x_cat_end, 60);
y_cat_v = interp1(x, y, x_cat_v);
x_an_v  = linspace(x_an_start,  x_an_end,  60);
y_an_v  = interp1(x, y, x_an_v);

props_2d = {T_2d,       V_2d,    rho_2d};
clbls    = {'T (degC)', 'V (m/s)', '\rho (kg/m^3)'};
titls_f4 = {'Static Temperature', 'Flow Velocity', 'Gas Density'};
cmaps    = {'hot',      'turbo',  'winter'};

for s = 1:3
    subplot(3,1,s)

    pcolor(X_g, Y_g, props_2d{s}); shading interp
    colormap(gca, cmaps{s}); hold on

    plot(x_g,  y_wall, 'k-', 'linewidth', 2, 'HandleVisibility', 'off')
    plot(x_g, -y_wall, 'k-', 'linewidth', 2, 'HandleVisibility', 'off')

    xline(x_arc_start, 'w--', 'linewidth', 1.2, 'HandleVisibility', 'off')
    xline(x_arc_end,   'w--', 'linewidth', 1.2, 'HandleVisibility', 'off')

    % Cathode
    fill([x_cat_v, fliplr(x_cat_v)], [ y_cat_v, fliplr( y_cat_v + t_elec)], ...
        'r', 'EdgeColor','k','linewidth',0.5,'HandleVisibility','off')
    fill([x_cat_v, fliplr(x_cat_v)], [-y_cat_v, fliplr(-y_cat_v - t_elec)], ...
        'r', 'EdgeColor','k','linewidth',0.5,'HandleVisibility','off')

    % Anode
    fill([x_an_v, fliplr(x_an_v)], [ y_an_v, fliplr( y_an_v + t_elec)], ...
        [1 0.72 0], 'EdgeColor','k','linewidth',0.5,'HandleVisibility','off')
    fill([x_an_v, fliplr(x_an_v)], [-y_an_v, fliplr(-y_an_v - t_elec)], ...
        [1 0.72 0], 'EdgeColor','k','linewidth',0.5,'HandleVisibility','off')

    cb = colorbar;
    cb.Label.String   = clbls{s};
    cb.Label.FontSize = 11;
    xlabel('x (mm)', 'fontsize', 12)
    ylabel('y (mm)', 'fontsize', 12)
    title(titls_f4{s}, 'fontsize', 13)
    xlim([0, L]); axis tight
    set(gca, 'fontsize', 11)
end

if sim_mode == 1
    sgtitle(sprintf('2D Flow Distribution  |  Arc: %.0f V x %.0f A = %.1f kW  |  %s', ...
        V_arc, I_arc, P_arc/1000, gas_name), 'fontsize', 13)
else
    sgtitle(sprintf('2D Flow Distribution  |  Inverse Design: %.0f m/s, %.0f degC  |  %s', ...
        V_exit_target, T_exit_target_C, gas_name), 'fontsize', 13)
end
