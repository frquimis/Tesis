clear all;
clc; 

%valores iniciales
M=[8 0 0;0 2 0;0 0 5];
C=[25 -15 -10;-15 35 -20;-10 -20 30];
K=[3600 -1200 -2400;-1200 2700 -1500;-2400 -1500 3900];
q0=[0;0.03;0];

% vector que pertenece al espacio nuo de la matriz K
Uo = null(K);
%frecuencia Wo 
Wo= 0;

Uo_normalized = Uo / sqrt(transpose(Uo) * M * Uo);

 
%matriz de restricciones S
S = [1 0; 0 1; -1.6 -0.4];

% el sistema positivo definido reducido es M'q'+ C'q'+K'q'
Mp =S.' * M * S;
Cp=S.' * C * S;
Kp=S.' * K * S;

B0=Uo_normalized.'* M * [0 ; 0.03; 0];


qrbpun=[0;0;0];

qrbt  =Uo_normalized*B0;
qrb0=qrbt; 

qnrb0=q0-qrb0;
qnrb0pun=[0;0;0];

parte1=-inv(Mp)*Kp;
parte2=-inv(Mp)*Cp;

A=[zeros(2, 2) eye(2);parte1 parte2];


% Valores propios proporcionados

% Calcular los vectores propios correspondientes
[V, D] = eig(A);

x=eig(A);


% Extraer los vectores propios derechos correspondientes a los valores propios proporcionados
X = V;
% Calcular los vectores propios izquierdos correspondientes
Y = transpose(inv(X));
%variale=transpose(Y)*A *X;
% Verificar la relación Y' X = I
%identity_matrix = transpose(Y)* X;



x0 = [qnrb0(1:2); zeros(2,1)];

syms t
matrizdiagonal=diag([exp(x(1)*t),exp(x(2)*t),exp(x(3)*t),exp(x(4)*t)]);


resultado2 = X * matrizdiagonal * transpose(Y) *x0;

qnrbt= S * resultado2(1:2);

% Evaluar las expresiones en un valor específico de t y convertir en decima

cifras_significativas = 4;

% Convierte la matriz simbólica a decimales con vpa
decimal_matrix = vpa(qnrbt, cifras_significativas);

% Muestra la matriz en formato decimal
disp(decimal_matrix);
% Extraer las partes real e imaginaria de la matriz original
decimal_matrix_real = real(decimal_matrix);
decimal_matrix_imag = imag(decimal_matrix);

% Aplicar la condición usando isAlways para evaluar con precisión
umbral = 1e-4; % Ajusta el umbral según sea necesario
mask_real = isAlways(abs(decimal_matrix_real) < umbral);
mask_imag = isAlways(abs(decimal_matrix_imag) < umbral);

% Redondear los valores cercanos a cero
decimal_matrix_real(mask_real) = 0;
decimal_matrix_imag(mask_imag) = 0;

% Unir las partes real e imaginaria nuevamente
decimal_matrix_combined = decimal_matrix_real + 1i * decimal_matrix_imag;

% Mostrar la matriz resultante
disp(decimal_matrix_combined);



% Muestra más resultados si es necesario
qt=qrbt+qnrbt;

% Define los valores de t
t_values = 0:0.1:0.6;

% Inicializa matrices para almacenar los resultados
qt_values = zeros(length(t_values), 3);

% Itera sobre los valores de t
for i = 1:length(t_values)
    % Calcula el valor de qt para el valor de t actual
    t = t_values(i);
    matrizdiagonal=diag([exp(x(1)*t),exp(x(2)*t),exp(x(3)*t),exp(x(4)*t)]);
    resultado2 = X * matrizdiagonal * transpose(Y) *x0;
    qnrbt= S * resultado2(1:2);
    qt_values(i, :) = (qrbt + qnrbt).';
end

% Interpolación spline para suavizar la curva
t_interp = linspace(min(t_values), max(t_values), 100); % 100 puntos para la interpolación
qt_interp = interp1(t_values, qt_values, t_interp, 'spline');

% Grafica los componentes de qt interpolados
figure;
plot(t_interp, qt_interp(:, 1), 'r-', 'LineWidth', 2);
hold on;
plot(t_interp, qt_interp(:, 2), 'g-', 'LineWidth', 2);
plot(t_interp, qt_interp(:, 3), 'b-', 'LineWidth', 2);
xlabel('t');
ylabel('Valor de qt');
title('Gráfica de qt en función de t');
legend('q1', 'q2', 'q3');
grid on;
