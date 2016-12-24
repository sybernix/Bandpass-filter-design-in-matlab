%Design of FIR Digital Filter
clc;
close all;
%parameters
Ap = 0.14; %dB maximum passband ripple
Aadesirable = 52; % dB minimum stopband attenuation
wp1 = 500; %rad/s lower passband edge
wp2 = 800; %rad/s upper passband edge
wa1 = 400; % rad/s lower stopband edge
wa2 = 950; %rad/s upper stopband edge
ws = 2400; %rad/s sampling frequency
T = 2*pi/ws;

Bt = min(wp1-wa1, wa2-wp2);
wc1 = wp1 - Bt/2;
wc2 = wp2 + Bt/2;

dp = ((10^(0.05*Ap))-1)/((10^(0.05*Ap))+1);
da = 10^(-0.05*Aadesirable);
d = min(dp, da);

Aa = -20*log10(d); %actual Aa of the filter

%kaiser window
if Aa<=21
    alpha = 0;
elseif Aa<=50
        alpha = (0.5842*((Aa-21)^0.4)) + (0.07886*(Aa-21));
else
    alpha = 0.1102*(Aa - 8.7);
end

if Aa<=21
    D = 0.9222;
else
    D = (Aa - 7.95)/14.36;
end

N = ceil((ws*D/Bt)+1);

if mod(N,2) == 0
    N=N+1;
end

wk = zeros(N,1);

for n = -(N-1)/2:(N-1)/2
   beta = alpha * (1 - (2*n/(N-1))^2)^0.5;
   numerator = myBessel(beta);
   denominator = myBessel(alpha);
   wk(n+(N-1)/2+1) = numerator/denominator;
end

stem(wk);

h = zeros(N,1);
h(38) = (2/ws)*(wc2-wc1);
for n = -(N-1)/2:(N-1)/2
    if n==0
        h(n+(N-1)/2+1) = (2/ws)*(wc2-wc1);
    else
        h(n+(N-1)/2+1) = (1/(n*pi)) * (sin(wc2*n*T) - sin(wc1*n*T));
    end
end

figure;
stem(h);

fil = h.*wk;
filim =figure;
stem(fil);
title('Impulse Response');
xlabel('n + (N-1)/2');
ylabel('h[n]');
grid on;
saveas(filim, 'q4.png');

[amp, digiFreq] = freqz(fil);

analogFreq = digiFreq*ws/(2*pi);
ampdb = 20*log10(abs(amp));

fr = figure;
plot(analogFreq, ampdb);
axis([wp1 wp2 -0.05 0.05]);
title('Magnitude Response in Passband');
xlabel('w in rad/s');
ylabel('|H(w)|');
grid on;
saveas(fr, 'q3.png');

x = zeros(300,1);
w1 = wa1/2;
w2 = (wp2+wp1)/2;
w3 = (ws/2+wa2)/2;

for n = 1:300
    x(n) = sin (w1*n*T) + sin (w2*n*T) + sin (w3*n*T);
end

[amp, digiFreq] = freqz(x);

analogFreq = digiFreq*ws/(2*pi);
ampdb = 20*log10(abs(amp));

fr = figure;
plot(analogFreq, ampdb);
title('Input Signal in Frequency Domain');
xlabel('w in rad/s');
ylabel('|X(w)|');
grid on;
saveas(fr, 'q6.png');

lenin = length(x);
lenh = length(fil);
lenfft = lenin+lenh-1;

IN = fft(x, lenfft);
H = fft(fil, lenfft);
OUT = H.*IN;
out = ifft(OUT, lenfft);

threeplots = figure;
subplot(3,1,1);
stem(x);
title('Input Signal');
xlabel('n');
ylabel('x[n]');
grid on;


e = zeros(300,1);

for n = 1:300
    e(n) = sin (w2*n*T);
end
subplot(3,1,2);
stem(e);
axis([0 300 -2 2]);
title('Expected Output');
xlabel('n');
ylabel('e[n]');
grid on;

subplot(3,1,3);
stem(out);
axis([0 300 -2 2]);
title('Output');
xlabel('n');
ylabel('Y[n]');
grid on;

figure;

[amp, digiFreq] = freqz(out);

analogFreq = digiFreq*ws/(2*pi);
ampdb = 20*log10(abs(amp));



w = ws*(1-lenfft/2:lenfft/2)/lenfft;
IN1 = abs(fftshift(IN));
OUT1 = abs(fftshift(OUT));
E = fft(e, lenfft);
E1 = abs(fftshift(E));
figure;
subplot(3,1,1);
plot(w, IN1);
title('Frequency Spectrum of Input Signal');
xlabel('w in rad/s');
ylabel('|X(w)|');
grid on;

subplot(3,1,2);
plot(w, E1);
title('Frequency Spectrum of Expected Output Signal');
xlabel('w in rad/s');
ylabel('|E(w)|');
grid on;

subplot(3,1,3);
plot(w, OUT1);
title('Frequency Spectrum of Output Signal');
xlabel('w in rad/s');
ylabel('|Y(w)|');
grid on;




