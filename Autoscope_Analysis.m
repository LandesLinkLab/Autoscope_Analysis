%THIS MATLAB CODE CAN BE USED WHEN THE AUTOMATED BACKGROUND MEASUREMENT PER PARTICLE WAS PERFORMED.


%Clear all previous files/variables in the workspace and close previous figures
clear all
close all

%Define path, where the input data is taken from and addpath, where the
%Lorentzian function is placed
path = '\\smb.researchdata.illinois.edu\illinois-flandes\2Data\Autoscope\Ovi\Trial\New_Data_Old_WL_DC';
savepath = '\\smb.researchdata.illinois.edu\illinois-flandes\2Data\Autoscope\Ovi\Trial\New_Data_Old_WL_DC';


cd(path) %changes the working directory to the specified location

b = 3; %First number of the spectra you are looking at
f = 3; %Last number of the spectra you are looking at
d = '092024'; %Date of the spectrum
w = '001'; %Number of White Light spectrum
dc = '001'; %Number of Darkcurrent spectrum
region = 'all'; %Dataset region being analysed


%Load Input data, convert to cell format and transpose.
%Average the dark current, and white counts (if need be)

WhiteLightName = [d,'White_',w,'.txt'];
DarkCurrentName = [d,'Darkcurrent',dc,'.txt'];
White = mean(load(WhiteLightName)',2);
White = White(1:1022,:); %Extracts the firs 1022 values
Darkcurrent = mean(load(DarkCurrentName)',2);
Darkcurrent = Darkcurrent(1:1022,:);
NumberOfValues = length(White);
L = zeros([f-b+1, 1]); %Initialize the array for storing resonance positions

%Set the range of the wavelength that was looked at and define the
%increments for the x-axis
MinWave = 435.778;
MaxWave = 860.09;
Range = MaxWave - MinWave;
Steps = Range / NumberOfValues;
Steps = round(Steps*10000)/10000;
MinWaveNew1 = MinWave;
MaxWaveNew1 = MaxWave;
MinWaveNew2 = round(MinWaveNew1*10000)/10000;
MaxWaveNew2 = round(MaxWaveNew1*10000)/10000; %Rounds values to avoid floating point precision errors


%Start to analyze the data: Each spectrum file is analyzed using the
%background file with the same number.

    for z = b:f   %Iterates over each spectrum from b to f 
    h = z-b+1;
        
            
    if z < 10 
        filename = [d,'Spec_00',num2str(z),'.txt'];  %Constructs filenames for spectra data files
    elseif z < 100 
        filename = [d,'Spec_0',num2str(z),'.txt'];    
    else
        filename = [d,'Spec_',num2str(z),'.txt'];
    end
     
    if(~exist(filename, 'file')), continue; end   %Check if the file exists, if not, skips it
    
    Data = load(filename)';  %Loads and transposes the data
    NAcqu = size(Data,2)/3;
    Bck1Max = 3*(NAcqu-1)+1;
    PtcMax = 3*(NAcqu-1)+2;
    Bck2Max = 3*(NAcqu-1)+3;
    
    Background1 = mean(Data(:,1:3:Bck1Max),2);  %First column of data
    Particle = mean(Data(:,2:3:PtcMax),2);       %Second column of data
    Background2 = mean(Data(:,3:3:Bck2Max),2);    %Third column of data

    
    Check = mean(Background1) < mean(Background2);  %Checks which bakcground signal is lower
        
    Sample = ['Particle_', num2str(filename)];
    
%     %Calculate the signal
if mean(Background1) < mean(Background2)
    
    Signal = (Particle - Background1)./(White - Darkcurrent);  %Computes signal by normalizing the particle intensity using white light and dark current
    else
        Signal = (Particle - Background2)./(White - Darkcurrent);
end       
    %Definition of the axes
    x = MinWaveNew2:Steps:MaxWaveNew2;
    y = Signal;
    x = 1239.8./x;   %converting nm to eV
    %save('ephoton_interp','x');
%   

    %Quick check if irregularities in the background causes weird spectral shapes
        close all
        figure(1)
    plot(x,Background1)
    hold on
    plot(x,Background2)
    hold on
    plot(x,Particle)
    hold off
    legend('back1','back2','sig')
    figure(2)
    plot(x,White)
    hold on
    plot(x,Darkcurrent)
    hold on
    plot(x,White-Darkcurrent)
    legend('WL','DC','WL-DC')
%       
% Use this code when the Lorentzian fit is off:
low = 550;
high = 660;
lowCut = round((low-MinWaveNew2)/Steps);
highCut = round((MaxWaveNew2-high)/Steps);
x1 = x(lowCut:length(x)-highCut);
y1 = y(lowCut:length(y)-highCut);

%code to choose fitting window based on spectra file
vals = [50,51,52,53,73,74,87,97,115,121,122,126,130,153,158]; %spectra files numbers that need refiting 
    if ismember(z,vals)
        xlow = 1.7;
        xhigh = 2.15;
    else
        xlow = 1.5;
        xhigh = 2.4;
    end

%     Calculating the Lorentzian fit
% code to close in fitting window
 [xvallow ,xindlow] = (min(abs(x-xlow)));
 [xvalhigh ,xindhigh] = (min(abs(x-xhigh)));
%     


    [param_1,param_2]= fn_lorentz_fit(x',y,1,1); %Uses fn_lorentz_fit function to fit the spectrum with a Lorentzian curve.
    a1 = param_1.a1;
    b1 = param_1.b1;
    c1 = param_1.c1;
    resonance = b1;  %Extracts resonance energy (b1) and full-width at half-maximum (FWHM, c1)
    FWHM = c1;
    r_list = param_2.rsquare;
    lorentz_fit =(2*a1/pi).*(c1./(4*(x'-b1).^2+c1.^2)); %Computes Lorentzian fit function
      

    L(h) = (resonance);
    F(h) = FWHM;


    %Calculate the SNR (according to Hyperspectral code):
    AllSigyNoi =[];
    diff = y-lorentz_fit;
    Noi = std(diff);
    [Notneeded, IndiMax]=min(abs(x-resonance));
    Sigy = y(IndiMax);
    SnN= Sigy/Noi;
    AllSigyNoi = [AllSigyNoi,[Noi;Sigy;SnN]];


    % %Calculate the noise and the SNB (not SNR:
    % Backgrounddifference = (Background1-Background2)./(White - Darkcurrent);
    % Noise = std2(Backgrounddifference);
    % SNR_hm = max(Signal)./min(Noise);
    % SNR = max(Signal)./std(Signal);

    %Plotting the Graph, consistent of the actual data and the fit.
    %Creating a legend with resonance wavelenght, FWHM and Signal to noise
    %ratio.
    figure1 = figure;
    hold all
    plot(x,y,'b','linewidth',3)
    plot(x,lorentz_fit,'k--','linewidth',3)
    xlabel('Energy (eV)','fontsize',32)
    ylabel('Scattering','fontsize',32)
    set(gca,'FontSize',28,'box','on')
    xlim([1.4 2.5])
   
    position = (MaxWave-round(resonance));
    if position > 200

    text(0.65,0.9,['\lambda_m_a_x = ',num2str(round(resonance,2)), ' eV'],'fontsize',15,'Units','normalized')
    text(0.65,0.78,['\Gamma = ',num2str(round(FWHM,2)), ' eV'],'fontsize',15,'Units','normalized')
    % text(0.65,0.66,['S/N = ',num2str(round(SNR))],'fontsize',15,'Units','normalized')
    text(0.65,0.65,['S/N (HS)= ',num2str(round(SnN))],'fontsize',15,'Units','normalized')
    else
    text(0.05,0.9,['\lambda_m_a_x = ',num2str(round(resonance,2)), ' eV'],'fontsize',15,'Units','normalized')
    text(0.05,0.78,['\Gamma = ',num2str(round(FWHM,2)), ' eV'],'fontsize',15,'Units','normalized')
    % text(0.05,0.66,['S/N = ',num2str(round(SNR))],'fontsize',15,'Units','normalized')
    text(0.65,0.65,['S/N (HS)= ',num2str(round(SnN))],'fontsize',15,'Units','normalized')
    end
cd(savepath)
  
    hgsave(Sample,'-v6')
    saveas(figure1,[Sample,'.jpg'])

    % Convert resonance and linewidth back to nm
    resonance_nm = 1239.8 / resonance;
    FWHM_nm = (1239.8 * FWHM) / (resonance^2);
x_wavelength = 1239.8 ./ x; % Convert back to wavelength

figure4= figure;
plot(x_wavelength, y, 'r', 'linewidth', 3)
hold on
plot(x_wavelength, lorentz_fit, 'k--', 'linewidth', 3)
xlabel('Wavelength (nm)', 'fontsize', 32)
ylabel('Scattering', 'fontsize', 32)
set(gca, 'FontSize', 28, 'box', 'on')
xlim([450 850])
hgsave(Sample,'-v6')
    saveas(figure1,[Sample,'.jpg'])

% Display resonance and linewidth in nm
text(0.05, 0.9, ['\lambda_{max} = ', num2str(round(resonance_nm,2)), ' nm'], 'fontsize', 15, 'Units', 'normalized')
text(0.05, 0.78, ['\Gamma = ', num2str(round(FWHM_nm,2)), ' nm'], 'fontsize', 15, 'Units', 'normalized')
% Save the figure
saveas(figure4, 'Scattering_Wavelength_Energy.jpg')

    
%     
exnumber = z-b+3;
exname=['ScatteringData' d '_',region,  '.xlsx'];
xlswrite(exname,(resonance),'Sheet1',['R' num2str(exnumber)])  %Saves extracted resonance and FWHM values to an Excel file
xlswrite(exname,(FWHM),'Sheet1',['S' num2str(exnumber)])
xlswrite(exname,round(SNR),'Sheet1',['T' num2str(exnumber)])
xlswrite(exname,(param_2.rsquare),'Sheet1',['U' num2str(exnumber)])
%xlswrite(exname,round(SNR_hm),'Sheet1',['V' num2str(exnumber)])
% processed_spectra.x(:,h) = x;
% processed_spectra.y(:,h) = y;
% processed_spectra.fit(:,h) = lorentz_fit;


 cd(path)    
    end

L = nonzeros(L)';
F = nonzeros(F)';
% save('processed_spectra','-struct','processed_spectra');   

% Make image for the resonance position
figure2 = figure;
av = mean(L);
st = std(L);
h = histogram(L,10);
h.BinWidth = 0.0269;
c = h.BinWidth;
xlabel('Resonance Position (eV)','fontsize',30)
ylabel('Occurence','fontsize',30)
text(0.65,0.66,['Average: ',num2str(av)],'fontsize',15,'Units','normalized')
text(0.65,0.54,['STD: ',num2str(st)],'fontsize',15,'Units','normalized')
set(gca,'fontsize',25)
hgsave(strcat('resonanceposition_',region),'-v6')
saveas(figure2,['resonanceposition_',region,'.jpg'])
%save(L)


%Make image for the fwhm
figure3 = figure;
av1 = mean(F);
st1 = std(F);
h2 = histogram(F,10);
h2.BinWidth = 0.0390;
c2 = h2.BinWidth;
xlabel('FWHM (eV)','fontsize',30)
ylabel('Occurence','fontsize',30)
text(0.65,0.66,['Average: ',num2str(av1)],'fontsize',15,'Units','normalized')
text(0.65,0.54,['STD: ',num2str(st1)],'fontsize',15,'Units','normalized')
set(gca,'fontsize',25)
hgsave(strcat('fwhm_',region),'-v6')
saveas(figure3,['fwhm_',region,'.jpg'])