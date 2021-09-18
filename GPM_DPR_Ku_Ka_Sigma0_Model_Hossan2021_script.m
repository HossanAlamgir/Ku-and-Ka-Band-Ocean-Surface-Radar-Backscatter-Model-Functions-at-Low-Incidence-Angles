%%% 
% Hossan, A.; Jones, W.L. Ku- and Ka-Band Ocean Surface Radar Backscatter Model Functions
% at Low-Incidence Angles Using Full-Swath GPM DPR Data. Remote Sens. 2021, 13, 1569. https://doi.org/10.3390/rs13081569
% The model functions presented in this paper along with their corresponding bin average measurements,
% polynomial coefficients and SST correction factors are publicly available in
% https://github.com/HossanAlamgir/Ku-and-Ka-Band-Ocean-Surface-Radar-Backscatter-Model-Functions-at-Low-Incidence-Angles.git,
%%%
%% Load the Model Coefficients (either from Table A1-A4 of Appendix A in the paper of from the GitHub repository)
clear
load('GPM_DPR_Ku_Ka_Sigma0_Model_Hossan2021.mat')
clear modku modka
%% Define Wind Speeds and Compute A0, A1,and A2 Coefficients as Function of WS Using the Equations 4,6, and 7 respectively 
% Let x represent ws
x=[1:1:20]'; % for viz, I used 1 deg steps, but the smaller, better
xl=log10(x);
% for a0 poly3 vs log(ws); f(x) = p1*x^3 + p2*x^2 + p3*x + p4  dB
fa0 = a01.*xl.^3 + a02.*xl.^2 + a03.*xl + a04;
fc0 = c01.*xl.^3 + c02.*xl.^2 + c03.*xl + c04;

% for a1 poly3 vs ws; f(x) = p1*x^3 + p2*x^2 + p3*x + p4
fa1 = a11.*x.^3 + a12.*x.^2 + a13.*x + a14;
fc1 = c11.*x.^3 + c12.*x.^2 + c13.*x + c14;

% revA: f(x) = p1*x^7 + p2*x^6 + p3*x^5 + p4*x^4 + p5*x^3 + p6*x^2+p7*x+p8
fa2 = a21.*x.^7 + a22.*x.^6 + a23.*x.^5 + a24.*x.^4 + a25.*x.^3 + a26.*x.^2 + a27.*x + a28;
fc2 = c21.*x.^7 + c22.*x.^6 + c23.*x.^5 + c24.*x.^4 + c25.*x.^3 + c26.*x.^2 + c27.*x + c28;

%% Generate the Model (either the residual using eq. 5, or the full model including eq. 2 in the paper)
% The Fourier coefficients were developed in dB
for c=1:361
    chi=c-1;
modku(:,:,c) =  fa0+fa1.*cosd(chi)+fa2.*cosd(2*chi);
modka(:,:,c) =  fc0+fc1.*cosd(chi)+fc2.*cosd(2*chi);
end
modku=permute(modku,[1 3 2]); % Ku model for 1-20 m/s WS (1 m/s step) X 0-360 deg. rel. WD Chi (1 deg step) X 1-25 PR beams (~0.76 deg EIA step, 18-0 deg.)
modka=permute(modka,[1 3 2]); % Ka model for 1-20 m/s WS (1 m/s step) X 0-360 deg. rel. WD Chi (1 deg step) X 1-25 PR beams (~0.76 deg EIA step, 18-0 deg.)

mod_anom_ku=modku-nanmean(modku,2); % Ku model anomaly for 1-20 m/s WS (1 m/s step) X 0-360 deg. rel. WD Chi (1 deg step) X 1-25 PR beams (~0.76 deg EIA step, 18-0 deg.)
mod_anom_ka=modka-nanmean(modka,2); % Ka model anomaly for 1-20 m/s WS (1 m/s step) X 0-360 deg. rel. WD Chi (1 deg step) X 1-25 PR beams (~0.76 deg EIA step, 18-0 deg.)


%% Visualizatioina and Verification
% Case 1
lw=1.5;
ms=12;
ls='-';
chi=10:10:350;
chim=0:1:360;
figure
plot(chim,mod_anom_ku(10,:,20),'c--','LineWidth',lw,'DisplayName',[num2str(mean_eia_ku(20)),char(176)])
grid
hold
plot(chim,mod_anom_ku(10,:,17),'b-','LineWidth',lw,'DisplayName',[num2str(mean_eia_ku(17)),char(176)])
plot(chim,mod_anom_ku(10,:,14),'r-','LineWidth',lw,'DisplayName',[num2str(mean_eia_ku(14)),char(176)])
plot(chim,mod_anom_ku(10,:,11),'k-','LineWidth',lw,'DisplayName',[num2str(mean_eia_ku(11)),char(176)])
legend;
xlim([0 360])
set(gca,'XTick',[0:90:360] );
set(gca,'XTickLabel',[0:90:360] );   
title('Ku Sig0 Model Anomaly @ WS = 10 m/s')
ylabel('\sigma\circ Residual [dB]','FontSize',12,'Fontweight','bold');
xlabel('Relative Wind Direction [deg.]','FontSize',10,'Fontweight','bold');
% Thank you!
% Alamgir Hossan, CFRSL, UCF, 9/18/2021 
% ah@knights.ucf.edu

  