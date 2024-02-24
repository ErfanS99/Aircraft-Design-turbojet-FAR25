clc
clear
close all;
tic
%% Finding A & B , ploting Regression Line
W_E_e = xlsread('AD', 1, 'R2:R42');
W_TO_to = xlsread('AD', 1, 'Q2:Q42');
name = xlsread('AD', 1, 'P2:P42');
show = [name, W_E_e, W_TO_to];
%
W_E = log10(W_E_e)';
W_TO = log10(W_TO_to)';
%
[R, B, A] = regression(W_TO, W_E);  
format long
fprintf("A: %d\n",A)
fprintf("B: %d\n",B)
fprintf("R: %d\n",R)
figure, p = plotregression(W_E, W_TO);
xlim([4.7 5.5]);
ylim([4.85 5.8]);
legend off
title("\fontsize{14}\fontname{Times}R^2: ",R);
xlabel("\fontsize{14}\fontname{Times}Take off Weight");
ylabel("\fontsize{14}\fontname{Times}Empty Weight");
%% Error Estimating
Error = 1;
W_TO = 0;
format long;
for W_TO_guess= 50000:0.01:200000
    W_E1 = regAns(W_TO_guess);
    W_E2 = berguetAns(W_TO_guess);
    delta = abs(W_E1 - W_E2);
    if W_E1 > W_E2
        x = delta / W_E1;
    else
        x = delta / W_E2;
    end
    if x < 0.005
        if x < Error
            Error = x;
            W_TO = W_TO_guess;
        end
    end
end
fprintf('\nError: %d\n',Error) ;
fprintf('Max Take-off Weight: %d\n',W_TO) ;
fprintf('Empty Weight: %d\n',berguetAns(W_TO));
%% Sensitivities
% Defining Variables
syms W_TO y W_PL dW_TO dy Sen_W_PL Sen_Range Sen_Endurance Sen_Speed Sen_SFC Sen_L_Over_D
M_ff = 0.791743433225269;
M_res = 0;
M_tfo = 0.005;
W_PL = (150 * 175) + (150 * 40) + (5 * 30);
W_crew = 5 * 175;
R = 2205;
E = 0.5;
C_j_cr = 0.5;
C_j_ltr = 0.6;
L_over_D_cr = 16; 
L_over_D_ltr = 15;
V_cr = 458.88;
A = 0.0933;
B = 1.0329;
C_1 = (1 - (1 + M_res) * (1 - M_ff) - M_tfo);
C_2 = (M_ff * (1 + M_res) - M_tfo -M_res);
D = (W_PL + W_crew);
W_TO = vpasolve(A + B * log10(C_1 * W_TO - D) - log10(W_TO)==0,W_TO);
F = - B * (W_TO ^ 2) * (1 / (C_2 * W_TO * (1 - B) - D)) * (1 + M_res) * M_ff;
%%
% Calculating Sensitivities
Sen_W_PL = (B * W_TO)/((D - C_1 *(1 - B) * W_TO));
Sen_W_E = (B * W_TO) / regAns(W_TO);
Sen_Range = (F * C_j_cr) / (V_cr * L_over_D_cr);
Sen_Endurance = F * C_j_ltr * (1 / L_over_D_ltr);
Sen_C_j_cr = F * R * (1 / (V_cr * L_over_D_cr));
Sen_C_j_ltr = (F * E) / L_over_D_ltr;
Sen_Speed = (F * (-R) * C_j_ltr) / (power(V_cr, 2) * L_over_D_cr);
Sen_L_Over_D_cr = F * (-R) * C_j_cr * (1 / (V_cr * (L_over_D_cr ^ 2)));
Sen_L_Over_D_ltr = F * (-E) * C_j_ltr * (1 /(L_over_D_ltr ^ 2));
%% 
% Tabeling The Results
Name = ["Sen_W_PL", "Sen_W_E","Sen_Range","Sen_Endurance","Sen_C_j_cr","Sen_C_j_ltr","Sen_Speed","Sen_L_Over_D_cr","Sen_L_Over_D_ltr"]';
Value = [Sen_W_PL, Sen_W_E, Sen_Range, Sen_Endurance, Sen_C_j_cr, Sen_C_j_ltr, Sen_Speed, Sen_L_Over_D_cr, Sen_L_Over_D_ltr]';
Value = double(Value);
format long
Sensitivities = table(Name, Value);
disp(Sensitivities);
%% Stall Sizing
V_S_clean = 421.916; % Unit: ft/s 250 knot 0 deg
V_S_L = 273.4016; % Unit: ft/s  162 knot 40 deg
V_S_TO = 354.4094; % Unit: ft/s  210 knot  10 deg
rho_cruise = 0.000632184; % Unit: slug/ft^3
rho_TOandL = 0.00226042; % Unit: slug/ft^3
C_L_max_clean = [1.2, 1.4, 1.6, 1.8];
C_L_max_takeoff = [1.6, 1.8, 2, 2.2];
C_L_max_landing = [1.8, 2, 2.2, 2.4, 2.6, 2.8];
figure;
disp("                           ""For Cruise Phase""         ");

subplot(3,1,1), hold on
for C_L = C_L_max_clean
W_over_S = 0.5 * rho_cruise * (V_S_clean ^ 2) * C_L;
fprintf('For C_L_max_clean = %d    ---->    W_over_S = %d\n',C_L,W_over_S);
xline(W_over_S, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{14}\fontname{Times}Stall Speed Sizing for V_s Clean (Flaps up)");
xlim([65, 105]);
xlabel("\fontsize{14}\fontname{Times}W/S [psf]");
ylabel("\fontsize{14}\fontname{Times}T/W [{lbf}/{lbs}]");
disp("                           ""For Take off Phase""         ");

subplot(3,1,2), hold on
for C_L = C_L_max_takeoff
W_over_S = 0.5 * rho_TOandL * (V_S_TO ^ 2) * C_L;
fprintf('For C_L_max_Takeoff = %d    ---->    W_over_S = %d\n',C_L,W_over_S);
xline(W_over_S, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{14}\fontname{Times}Stall Speed Sizing for V_s Take off (Flaps Down)");
xlim([220, 320]);
xlabel("\fontsize{14}\fontname{Times}W/S [psf]");
ylabel("\fontsize{14}\fontname{Times}T/W [{lbf}/{lbs}]");
disp("                           ""For Landing Phase""         ");

subplot(3,1,3), hold on
for C_L = C_L_max_landing
W_over_S = 0.5 * rho_TOandL * (V_S_L ^ 2) * C_L;
fprintf('For C_L_max_Landing = %d    ---->    W_over_S = %d\n',C_L,W_over_S);
xline(W_over_S, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{14}\fontname{Times}Stall Speed Sizing for V_s Landing (Flaps Down)");
xlim([150, 240]);
xlabel("\fontsize{14}\fontname{Times}W/S [psf]");
ylabel("\fontsize{14}\fontname{Times}T/W [{lbf}/{lbs}]");
%% Sizing of Take-off Field length
C_L_max_takeoff = [1.6, 1.8, 2, 2.2];
s_tofl = 6850;
sigma = 0.950885;
W_over_S = 1:0.001:500;
figure, hold on
for C_L_max = C_L_max_takeoff
    T_over_W = W_over_S./((s_tofl/37.5) .* sigma .* C_L_max);
    plot(W_over_S, T_over_W,"LineWidth",3);
    grid on
end
legend('(Cl_m_a_x)_T_O = 1.4', '(Cl_m_a_x)_T_O = 1.6', '(Cl_m_a_x)_T_O = 1.8', '(Cl_m_a_x)_T_O = 2.0');
title("\fontsize{14}\fontname{Times}Thrust Loading to Wing Loading (Take-off)");
xlabel("W/S");
ylabel("T/W");
hold off
%% Sizing of Landing Field length (1)
rho_TOandL = 0.00226042; % Unit: slug/ft^3
W_TO = 1.444106e+05;
W_L = [0.65 0.84 1.00] * W_TO;
fprintf('W_L = %d\n', W_L);
C_L_max_landing = [1.8, 2, 2.2, 2.4, 2.6, 2.8];
S_FL = 6850;
S_L = S_FL * 0.6;
V_A = sqrt(S_FL/0.3);
V_S_L = V_A/1.3;
sva = 0:0.01:40; % Square of V approach
sfl = (6/20) * sva;   % Square of V approach to FAR25 Landing Field Length
figure, plot(sva, sfl, ':','LineWidth', 3);
title("\fontsize{16}\fontname{Times}Square of V approach to FAR25 Landing Field Length");
xlabel("\fontsize{14}\fontname{Times}(V_A)^2  [(KN^2)10^-^3]");
ylabel("\fontsize{14}\fontname{Times}FAR25 Landing Field Length - S_F_L  [ft(10^-^3)]");
grid on
hold on
v_a_p = (V_A ^ 2) * (10 ^ (-3));
sfl_p = S_FL * (10 ^ (-3));
plot(v_a_p, sfl_p,'O','MarkerSize',7, 'MarkerFaceColor', 'r');
legend('\fontsize{14}\fontname{Times}Average Line', '\fontsize{14}\fontname{Times}Our Airplane');
%% Sizing of Landing Field length (2)
disp("                           ""For Landing Wing Loading""         ");
V_S_L = V_S_L * 1.6878098571;   % Changing the unit from KN to ft/s 
figure, hold on
for C_L = C_L_max_landing
    W_over_S_Landing = 0.5 * rho_TOandL * (V_S_L ^ 2) * C_L;
    fprintf('For C_L_max_Landing = %d    ---->    W_over_S = %d\n',C_L,W_over_S_Landing);
    xline(W_over_S_Landing, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{16}\fontname{Times}Wing Loading in Landing");
xlabel("\fontsize{16}\fontname{Times}(W/S)_L_a_n_d_i_n_g");
ylabel("\fontsize{16}\fontname{Times}(T/W)_L_a_n_d_i_n_g");
% xlim([60 120]);
%% Sizing of Landing Field length (3)
rho_TOandL = 0.00226042; % Unit: slug/ft^3
C_L_max_landing = [1.8, 2, 2.2, 2.4, 2.6, 2.8];
figure,
disp("                           ""Medium Aircraft Take-off Wing Loading""         ");
subplot(3,1,1), hold on
for C_L = C_L_max_landing
W_over_S_TakeOff = ((0.5 * rho_TOandL * (V_S_L ^ 2))/0.65) * C_L;
fprintf('For C_L_max_clean = %d    ---->    W_over_S_Take-off = %d\n',C_L,W_over_S_TakeOff);
xline(W_over_S_TakeOff, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{16}\fontname{Times}Landing Sizing for Medium Aircraft");
% xlabel("\fontsize{14}\fontname{Times}W/S [psf]");
ylabel("\fontsize{14}\fontname{Times}T/W");
xlim([115, 190]);

disp("                           ""Average Aircraft Take-off Wing Loading""         ");
subplot(3,1,2), hold on
for C_L = C_L_max_landing
W_over_S_TakeOff = ((0.5 * rho_TOandL * (V_S_L ^ 2))/0.84) * C_L;
fprintf('For C_L_max_Takeoff = %d    ---->    W_over_S_Take-off = %d\n',C_L,W_over_S_TakeOff);
xline(W_over_S_TakeOff, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{16}\fontname{Times}Landing Sizing for Average Aircraft");
% xlabel("\fontsize{14}\fontname{Times}W/S [psf]");
ylabel("\fontsize{14}\fontname{Times}T/W");
% xlim([90, 190]);

disp("                           ""Maximum Aircraft Take-off Wing Loading""         ");
subplot(3,1,3), hold on
for C_L = C_L_max_landing
W_over_S_TakeOff = ((0.5 * rho_TOandL * (V_S_L ^ 2))/1) * C_L;
fprintf('For C_L_max_Landing = %d    ---->    W_over_S_Take-off = %d\n',C_L,W_over_S_TakeOff);
xline(W_over_S_TakeOff, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{16}\fontname{Times}Landing Sizing for Maximum Aircraft");
xlabel("\fontsize{14}\fontname{Times}W/S [psf]");
ylabel("\fontsize{14}\fontname{Times}T/W");
% xlim([90, 125]);
%% ُClimb Sizing
W_TO = 1.444106e+05;
S = 1340.968;
a = -2.5229;
b = 1;
c = 0.0199;
d = 0.7531;
S_wet = power(10,c + (d * log10(W_TO)));
f = power(10, a + (b * log10(S_wet)));
C_D_0 = f / S;
m = 0.0100;
x = 100:0.1:10000;
x_1 = 0;
x_2 = -6.4;
x_3 = -12.7;
x_4 = -19;
y_1 = m * x + x_1;
y_2 = m * x + x_2;
y_3 = m * x + x_3;
y_4 = m * x + x_4;
p_1 = m * S_wet + x_1;
p_2 = m * S_wet + x_2;
p_3 = m * S_wet + x_3;
p_4 = m * S_wet + x_4;
figure, hold on
plot(x, y_1, 'LineWidth',3);
plot(x, y_2, 'LineWidth',3);
plot(x, y_3, 'LineWidth',3);
plot(x, y_4, 'LineWidth',3);
ylim([40 100]);
xlim([6000 10000]);
xline(S_wet, 'LineWidth',3,'Label','\fontsize{16}\fontname{Times}S_W_e_t','Color','r');
plot(S_wet, [p_1, p_2, p_3, p_4],'O','MarkerSize',9, 'MarkerFaceColor', 'r','MarkerEdgeColor','black');
text(S_wet,p_4,cellstr(append('    ',num2str(p_4))));
text(S_wet,p_3,cellstr(append('    ',num2str(p_3))));
text(S_wet,p_2,cellstr(append('    ',num2str(p_2))));
text(S_wet,p_1,cellstr(append('    ',num2str(p_1))));
legend('\fontsize{16}\fontname{Times}C_f = 0.0020', '\fontsize{16}\fontname{Times}C_f = 0.0030','\fontsize{16}\fontname{Times}C_f = 0.0040', '\fontsize{16}\fontname{Times}C_f = 0.0050');
xlabel('\fontsize{16}\fontname{Times}Wetted Area S_W_e_t ft^2');
ylabel('\fontsize{16}\fontname{Times}(Equivalent Parasite Area ft^2)*10^2');
grid on
%% Drag Polar
syms C_L
AR = 9.45;
e = [0.85 0.8 0.75 0.8 0.75];
deltaC_D_0 = [0 0.020 0.075 0.025 0.025];
C_D_0 = (C_D_0 + deltaC_D_0);
CC = 1./(pi.* AR .* e);
fprintf('For Clean Phase: %d+(%d)C_L^2\n', C_D_0(1), CC(1));
fprintf('For Take-ff (gear up) Phase: %d+(%d)C_L^2\n', C_D_0(2), CC(2));
fprintf('For Landing (gear up) Phase: %d+(%d)C_L^2\n', C_D_0(3), CC(3));
fprintf('For Take-off (gear down) Phase: %d+(%d)C_L^2\n', C_D_0(4), CC(4));
fprintf('For Landing (gear down) Phase: %d+(%d)C_L^2\n', C_D_0(5), CC(5));
disp('');
%% Climb requirements
C_L_max_clean = [1.2, 1.4, 1.6, 1.8];
C_L_max_takeoff = [1.6, 1.8, 2, 2.2];
C_L_max_landing = [1.8, 2, 2.2, 2.4, 2.6, 2.8];
% V_S_TO = 354.4094; % Unit: ft/s  210 knot  10 deg
% rho_airport = 0.00226042; % Unit: slug/ft^3
sigma = 0.950885;
C_D = (C_D_0 + (CC .* (C_L.^2)))';
N = 2;
z = 0.9;
figure;
% FAR 25.111 (OEI)
subplot(3,2,1), hold on;
for i = C_L_max_takeoff
    a = 1.2;
    L_over_D = double(subs(C_D(2) ./ C_L, C_L, i/a));
    CGR = 0.012;
    T_over_W = ((N / (N - 1)) * (L_over_D + CGR)) / sigma;
    yline(T_over_W/(sigma^z),"LineWidth", 1.5,"Label", i);
end
hold off;
title('\fontsize{16}\fontname{Times}FAR 25.111 (OEI)');
xlabel('\fontsize{14}\fontname{Times}Wing Loading');
ylabel('\fontsize{14}\fontname{Times}Thrust Loading');

% FAR 25.121 (OEI) (1)
subplot(3,2,2), hold on;
for i = C_L_max_takeoff
    a = 1.1;
    L_over_D = double(subs(C_D(3) ./ C_L, C_L, i/(a^2)));
    CGR = 0;
    T_over_W = ((N / (N - 1)) * (L_over_D + CGR)) / sigma;
    yline(T_over_W/(sigma^z), "Color",rand(1,3),"LineWidth",3,"Label", i);
end
hold off;
title('\fontsize{16}\fontname{Times}FAR 25.121 (1) (OEI)');
xlabel('\fontsize{14}\fontname{Times}Wing Loading');
ylabel('\fontsize{14}\fontname{Times}Thrust Loading');

% FAR 25.121 (OEI) (2)
subplot(3,2,3), hold on;
for i = C_L_max_takeoff
    a = 1.2;
    L_over_D = double(subs(C_D(2) ./ C_L, C_L, i/a));
    CGR = 0.024;
    T_over_W = ((N / (N - 1)) * (L_over_D + CGR)) / sigma;
    yline(T_over_W/(sigma^z), "Color",rand(1,3),"LineWidth",3,"Label", i);
end
hold off;
title('\fontsize{16}\fontname{Times}FAR 25.121 (2) (OEI)');
xlabel('\fontsize{14}\fontname{Times}Wing Loading');
ylabel('\fontsize{14}\fontname{Times}Thrust Loading');

% FAR 25.121 (OEI) (3)
subplot(3,2,4), hold on;
for i = C_L_max_clean
    a = 1.25;
    L_over_D = double(subs(C_D(1) ./ C_L, C_L, i/(a ^ 2)));
    CGR = 0.012;
    T_over_W = ((N / (N - 1)) * (L_over_D + CGR)) / sigma;
    yline(T_over_W/0.94/(sigma^z), "Color",rand(1,3),"LineWidth",3,"Label", i)
end
hold off;
title('\fontsize{16}\fontname{Times}FAR 25.121 (3) (OEI)');
xlabel('\fontsize{14}\fontname{Times}Wing Loading');
ylabel('\fontsize{14}\fontname{Times}Thrust Loading');

% FAR 25.119 (AEI)
subplot(3,2,5), hold on;
for i = C_L_max_landing
    a = 1.3;
    L_over_D = double(subs(C_D(5) ./ C_L, C_L, i/(a ^ 2)));
    CGR = 0.032;
    T_over_W = ((N / (N - 1)) * (L_over_D + CGR)) / sigma;
    yline((T_over_W/(W_L(2)/W_TO)/(sigma^z)), "Color",rand(1,3),"LineWidth",3,"Label", i)
end
hold off;
title('\fontsize{16}\fontname{Times}FAR 25.119 (AEI)');
xlabel('\fontsize{14}\fontname{Times}Wing Loading');
ylabel('\fontsize{14}\fontname{Times}Thrust Loading');

% FAR 25.121 (OEI)
subplot(3,2,6), hold on;
for i = C_L_max_landing
    a = 1.5;
    L_over_D = double(subs(C_D(5) ./ C_L, C_L, i / (a ^ 2)));
    CGR = 0.021;
    T_over_W = ((N / (N - 1)) * (L_over_D + CGR)) / sigma;
    yline(T_over_W/(W_L(2)/W_TO)/(sigma^z), "Color",rand(1,3),"LineWidth",3,"Label", i);
end
hold off;
title('\fontsize{16}\fontname{Times}FAR 25.119 (AEI)');
xlabel('\fontsize{14}\fontname{Times}Wing Loading');
ylabel('\fontsize{14}\fontname{Times}Thrust Loading');

title('\fontsize{18}\fontname{Times}Climb Sizing');
%% Linear Relation Between RC and h
% h = 37000; 
% t_cl = 15; % Unit min
% h_abs = 0.045; % Unit ft
% RC_0 = (h_abs ./ t_cl) .* log(1 ./ (1 - (h ./ h_abs)));
% h = 0:1:45000;
% RC_h = RC_0 .* (1 - (h ./ h_abs));
% sigma = 0.842744;
% figure, plot(RC_h, h, 'LineWidth', 3, 'Color', 'red');
% title('\fontsize{18}\fontname{Times}Linearized Rate-of-Climb With Altitude');
% xlabel('\fontsize{16}\fontname{Times}Rate-of-Climb [fpm]');
% ylabel('\fontsize{16}\fontname{Times}Altitude [ft]');
%% Time-to-Climb Sizing
% C_L_max_clean = [1.2, 1.4, 1.6, 1.8];
% h = 37000;
% RC_0 = (h_abs ./ t_cl) .* log(1 ./ (1 - (h ./ h_abs)));
% RC = RC_0 .* (1 - (h ./ h_abs));
% RCP = RC / 33000;
% W_over_S = 0:0.1:100;
% etha_p = 0.72;
% figure, hold on;
% for C_L = C_L_max_clean
%     C_D_clean = 1.645775e-02 + (4.160914e-02) .* C_L .^ 2;
%     temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./ C_D_clean) .* sqrt(sigma));
%     W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
%     plot(W_over_S, W_over_P, 'LineWidth', 2);
% end
% xlim([10 100])
% legend('\fontsize{14}\fontname{Times}C_L_{max} = 1.2' , '\fontsize{14}\fontname{Times}C_L_{max} = 1.4' , '\fontsize{14}\fontname{Times}C_L_{max} = 1.6' , '\fontsize{14}\fontname{Times}C_L_{max} = 1.8')
% title('\fontsize{18}\fontname{Times} Time to Climb Sizing');
% xlabel('\fontsize{16}\fontname{Times}W/S [psf]');
% ylabel('\fontsize{16}\fontname{Times}W/P [lb/hsp]');
% hold off
%% Ceiling
RC = 500/60;
W_over_S = W_TO/S;
CD0 = 0.018;
e = 0.8;
figure;
L_over_D_max = 0.5 * sqrt((pi * AR * e)/CD0);
v = sqrt((2 * (W_over_S)) / (rho_cruise * sqrt(CD0 * pi * AR * e)));
T_over_W = (RC/v) + (1 / L_over_D_max);
yline(T_over_W, 'LineWidth', 2);
xlim([10 100])
title('\fontsize{18}\fontname{Times} Ceiling Sizing');
xlabel('\fontsize{16}\fontname{Times}W/S [psf]');
ylabel('\fontsize{16}\fontname{Times}T/W [lbf/lbs]');
%% Manuvering
W_over_S = 0:0.001:200;
rho_cruise = 0.000632184;
AR = 9.45;
e = 0.85;
M = 0.8;
V_cr = 774.3948;
q = 0.5 * rho_cruise * (V_cr ^ 2);
n_max = 2.5;
T_over_W = ((C_D_0(1) .* q) ./ W_over_S) + ((W_over_S * (n_max^2))/(pi*AR*e*q));
figure, plot(W_over_S, T_over_W, 'LineWidth', 1.5, "Color", 'yellow');
ylim([0, 1]);
title('\fontsize{18}\fontname{Times}Sizing to Maneuvering');
xlabel('\fontsize{16}\fontname{Times} W/S [psf]');
ylabel('\fontsize{16}\fontname{Times} T/W [lbf/lbs]');
grid on;
%% Cruise
rho_cruise = 0.000632184;
W_over_S = 0:0.001:200;
AR = 9.45;
e = 0.85;
V_cr = 774.3948;
q = 0.5 * rho_cruise * (V_cr ^ 2);
T_over_W = C_D_0(1) .* q .* (1 ./ W_over_S) ...
    + (W_over_S./(pi .* AR .* e .* q));
figure, plot(W_over_S, T_over_W, 'LineWidth', 3, "Color", 'red');
title('\fontsize{18}\fontname{Times}Cruise Sizing');
xlabel('\fontsize{16}\fontname{Times} W/S [psf]');
ylabel('\fontsize{16}\fontname{Times} T/W [lb/hsp]');
xlim([0, 120]);
ylim([0, 0.4]);
%% Matching Diagram
V_S_clean = 421.916; % Unit: ft/s 250 knot 0 deg
% V_S_L = 273.4016; % Unit: ft/s  162 knot 40 deg
% V_S_TO = 354.4094; % Unit: ft/s  210 knot  10 deg
rho_cruise = 0.000632184; % Unit: slug/ft^3
rho_TOandL = 0.00226042; % Unit: slug/ft^3
C_L_max_clean = [1.2, 1.4, 1.6, 1.8];
C_L_max_takeoff = [1.6, 1.8, 2, 2.2];
C_L_max_landing = [1.8, 2, 2.2, 2.4, 2.6, 2.8];
W_TO = 1.444106e+05;
W_L = W_TO * 0.84;
rho_airport = 0.00226042; % Unit: slug/ft^3

figure, hold on;

% Stall
for C_L = C_L_max_clean
    W_over_S = 0.5 * rho_cruise * (V_S_clean ^ 2) * C_L;
    xline(W_over_S, "Color",'black',"LineWidth", 1.5,"Label",C_L);
end

% for C_L = C_L_max_takeoff
%     W_over_S = 0.5 * rho_TOandL * (V_S_TO ^ 2) * C_L;
%     xline(W_over_S, "Color", 'black', "LineWidth", 1.5, "Label",C_L);
% end

% for C_L = C_L_max_landing
%     W_over_S = 0.5 * rho_TOandL * (V_S_L ^ 2) * C_L;
%     xline(W_over_S, "Color",'black',"LineWidth",1.5,"Label",C_L);
% end

% Take off
s_tofl = 6850;
sigma = 1.0001;
W_over_S = 1:0.001:200;
for C_L_max = C_L_max_takeoff
    T_over_W = W_over_S./((s_tofl/37.5) .* sigma .* C_L_max);
    plot(W_over_S, T_over_W,"LineWidth", 1.5);
    grid on
end

% Landing
S_FL = 6850;
V_A = sqrt(S_FL/0.3);
V_S_L = V_A/1.3;
V_S_L = V_S_L * 1.6878098571;
for C_L = C_L_max_landing
W_over_S_TakeOff = ((0.5 * rho_TOandL * (V_S_L ^ 2))/0.84) * C_L;
fprintf('For C_L_max_Takeoff = %d    ---->    W_over_S_Take-off = %d\n',C_L,W_over_S_TakeOff);
xline(W_over_S_TakeOff, "Color",'blue',"LineWidth",1.5,"Label",C_L);
end

% Climb
syms C_L
C_D(1) = (C_D_0(1) + (CC(1) * (C_L^2)));
C_D(2) = (C_D_0(2) + (CC(2) * (C_L^2)));
C_D(3) = (C_D_0(3) + (CC(3) * (C_L^2)));
C_D(4) = (C_D_0(4) + (CC(4) * (C_L^2)));
C_D(5) = (C_D_0(5) + (CC(5) * (C_L^2)));
N = 2;

% FAR 25.111 (OEI)
z = 0.9;
a = 1.2;
L_over_D = double(subs(C_D(2) ./ C_L, C_L, max(C_L_max_takeoff) / a));
CGR = 0.012;
T_over_W = ((N / (N - 1)) * (L_over_D + CGR)) / sigma;
yline(T_over_W/(sigma^z),"LineWidth", 1.5,"Label", 'FAR 25.111 (OEI)')

% FAR 25.121 (OEI) (1)
a = 1.1;
L_over_D = double(subs(C_D(3) ./ C_L, C_L,  max(C_L_max_takeoff) / (a^2)));
CGR = 0;
T_over_W = ((N / (N - 1)) * (L_over_D + CGR)) / sigma;
yline(T_over_W/(sigma^z),"LineWidth", 1.5,"Label", 'FAR 25.121 (OEI) Transition')

% FAR 25.121 (OEI) (2)
a = 1.2;
L_over_D = double(subs(C_D(2) ./ C_L, C_L,  max(C_L_max_takeoff) / a));
CGR = 0.024;
T_over_W = ((N / (N - 1)) * (L_over_D + CGR)) / sigma;
yline(T_over_W/(sigma^z),"LineWidth", 1.5,"Label", 'FAR 25.121 (OEI) Second Segment')

% FAR 25.121 (OEI) (3)
a = 1.25;
L_over_D = double(subs(C_D(1) ./ C_L, C_L,  max(C_L_max_clean) / (a ^ 2)));
CGR = 0.012;
T_over_W = ((N / (N - 1)) * (L_over_D + CGR)) / sigma;
yline(T_over_W/0.94/(sigma^z),"LineWidth", 1.5,"Label", 'FAR 25.121 (OEI) En-Route')

% FAR 25.119 (AEI)
a = 1.3;
L_over_D = double(subs(C_D(5) ./ C_L, C_L, max(C_L_max_landing) / (a ^ 2)));
CGR = 0.032;
T_over_W = ((N / (N - 1)) * (L_over_D + CGR)) / sigma;
yline(T_over_W/(W_L/W_TO)/(sigma^z),"LineWidth", 1.5,"Label", 'FAR 25.119 (AEI) Approach')

% FAR 25.121 (OEI)
a = 1.5;
L_over_D = double(subs(C_D(5) ./ C_L, C_L, max(C_L_max_landing) / (a ^ 2)));
CGR = 0.021;
T_over_W = ((N / (N - 1)) * (L_over_D + CGR)) / sigma;
yline(T_over_W/(W_L/W_TO)/(sigma^z),"LineWidth", 1.5,"Label", 'FAR 25.121 (OEI) Landing')

% Cruise
e = 0.85;
V_cr = 774.3948;
q = 0.5 * rho_cruise * (V_cr ^ 2);
T_over_W = C_D_0(1) .* q .* (1 ./ W_over_S) ...
    + (W_over_S./(pi .* AR .* e .* q));
plot(W_over_S, T_over_W, 'LineWidth', 1.5, "Color", 'red');

% Ceiling
RC = 500/60;
W_over_S = W_TO/S;
CD0 = 0.018;
e = 0.8;
L_over_D_max = 0.5 * sqrt((pi * AR * e)/CD0);
v = sqrt((2 * (W_over_S)) / (rho_cruise * sqrt(CD0 * pi * AR * e)));
T_over_W = (RC/v) + (1 / L_over_D_max);
yline(T_over_W, 'LineWidth', 2);

% Manuvering
delta = 37000;
V_cr = 774.3948;
q = 0.5 * rho_cruise * (V_cr ^ 2);
W_over_S = 0:0.001:200;
n_max = 2.5;
T_over_W = ((C_D_0(1) .* q) ./ W_over_S) + ((W_over_S * (n_max^2))/(pi*AR*e*q));
plot(W_over_S, T_over_W, 'LineWidth', 1.5, "Color", 'yellow');

% Design Point
thrust_loading_1 = 0.31495;
wing_loading_1 = 101.2832;

thrust_loading_2 = 0.31495;
wing_loading_2 = 90.029632;

thrust_loading_3 = 0.31495;
wing_loading_3 = 78.7759284;

plot(wing_loading_1, thrust_loading_1, 'Marker','O', 'MarkerFaceColor', 'red', 'MarkerSize', 10);

plot(wing_loading_2, thrust_loading_2, 'Marker','O', 'MarkerFaceColor', 'blue', 'MarkerSize', 10);

plot(wing_loading_3, thrust_loading_3, 'Marker','O', 'MarkerFaceColor', 'blue', 'MarkerSize', 10);

text(wing_loading_1, thrust_loading_1, '\fontname{Times}   Design Point 1');
text(wing_loading_2, thrust_loading_2, '\fontname{Times}   Design Point 2');
text(wing_loading_3, thrust_loading_3, '\fontname{Times}   Design Point 3');

plot(106.9940271620881, 0.2542981542765046, 'Marker','diamond', 'MarkerFaceColor', 'blue', 'MarkerSize', 10);
text(106.9940271620881, 0.2542981542765046, '\fontname{Times}   Target Airplane');

fprintf('Design Point:\n');
fprintf('       Wing Loading: %d\n', wing_loading_1);
fprintf('       Thrust Loading: %d\n', thrust_loading_1);
xlabel('\fontsize{16}\fontname{Times}Wing Loading [psf]');
ylabel('\fontsize{16}\fontname{Times}Thrust Loading [lbf/lbs]');
title('\fontsize{18}\fontname{Times} Matching Diagram');

name_1 = "Stall CLmax = 1.2";
name_2 = "Stall CLmax = 1.4";
name_3 = "Stall CLmax = 1.6";
name_4 = "Stall CLmax = 1.8";
name_5 = "Take-off CLmax = 1.6";
name_6 = "Take-off CLmax = 1.8";
name_7 = "Take-off CLmax = 2";
name_8 = "Take-off CLmax = 2.2";
name_9 = "Landing Clmax = 1.8";
name_10 = "Landing Clmax = 2";
name_11 = "Landing Clmax = 2.2";
name_12 = "Landing Clmax = 2.4";
name_13 = "Landing Clmax = 2.6";
name_14 = "Landing Clmax = 2.8";
name_15 = "FAR 25.111 (OEI)";
name_16 = "FAR 25.121 (OEI) (1)";
name_17 = "FAR 25.121 (OEI) (2)";
name_18 = "FAR 25.121 (OEI) (3)";
name_19 = "FAR 25.119 (AEI)";
name_20 = "FAR 25.121 (OEI)";
name_21 = "Cruise";
name_22 = "Ceiling";
name_23 = "Manuvering";

legend(name_1, name_2, name_3, name_4, name_5, name_6, name_7, name_8, name_9, name_10, name_11, name_12, ...
    name_13, name_14, name_15, name_16, name_17, name_18, name_19, name_20, name_21, name_22, name_23);

ylim([0 1]);
xlim([0 200]);
grid on;

toc