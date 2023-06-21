%% main parameters
sample = 6e1 %the size of the sample for the random variable (h1,h2)
noIterations = 3000 % I run the programm 100 so that I can find the variance of different types of correlations. 
% channelModel = 'Ricean'
channelModel = 'nakgmi'
parametersRice = [0,0,1] % for Ricean channels. (m1, m2, sigma^2):
                    % m1=E(x), m2=E(y), and sigma^2 = E[(x-m1)^2]=E[(y-m2)^2] 
parameternakgmi = 2
%% main
R = zeros(noIterations, 1);
Rsqr = zeros(noIterations, 1);
Renv = zeros(noIterations, 1);
Rpower = zeros(noIterations, 1);

for i =1:noIterations
    % Collect samples
    h1 = zeros(sample,1); % initialisation
    h2 = zeros(sample,1);

    for j = [1:sample]

      if channelModel == 'Ricean'
            assert(length(parametersRice) == 3, 'error in parameters')
            m1 = parametersRice(1); m2 = parametersRice(2); sigma2 = parametersRice(3);
            x1 = random('Normal', m1, (sigma2/2)^.5); 
            y1 = random('Normal', m2, (sigma2/2)^.5); % the channel coefficient is h = x + yi
            h1(j) = complex(x1, y1);
            x2 = random('Normal',m1,(sigma2/2)^.5); 
            y2 = random('Normal',m2,(sigma2/2)^.5); % the channel coefficient is h = x + yi
            h2(j) = complex(x2, y2);
            a = (['Sample Size = ', num2str(sample), '       Ricean: K = ', num2str((m1^2+m2^2)/(sigma2))]);
      elseif channelModel == 'nakgmi'
          m = parameternakgmi
          phase1 = unifrnd(-pi,pi); phase2 = unifrnd(-pi,pi);
          gain1 = (1/(2*m)) * chi2rnd(2*m);   gain2= (1/(2*m)) * chi2rnd(2*m);
          h1(j) =sqrt(gain1) * exp(1i*phase1);
          h2(j) =sqrt(gain2) * exp(1i*phase2);
          a = (['Sample Size = ', num2str(sample), '       Nakagami: mu = ', num2str(m)]);
        end
    end
    R(i) = corr(h1,h2);
    Rsqr(i) = abs(R(i))^2;
    Renv(i) = corr(abs(h1),abs(h2));
    Rpower(i) = corr(abs(h1).^2,abs(h2).^2);
end
%% plots and results
figure

subplot(2,2,1)
histogram(real(R),'Normalization','probability')
legend(['\sigma^2 = ', num2str(var(real(R)))])
title(a)
subtitle('real(R)')

subplot(2,2,2)
histogram(Rsqr,'Normalization','probability')
legend(['\sigma^2 = ', num2str(var(Rsqr))])
subtitle('Rsqr')

% subplot(2,2,1)
% histogram(imag(R))
% title('image(R)')
subplot(2,2,3)
histogram(Renv,'Normalization','probability')
legend(['\sigma^2 = ', num2str(var(Renv))])
subtitle('Renv')
subplot(2,2,4)
histogram(Rpower,'Normalization','probability')
legend(['\sigma^2 = ', num2str(var(Rpower))])
subtitle('Rpower')


