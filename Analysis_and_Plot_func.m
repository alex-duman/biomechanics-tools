function Table = Analysis_and_Plot_func(trial,Plot)
% This function brings in raw EMG and kinematic data, processes it and will
% produce 2 figures if Plot term is set to 'Y' or 'Yes' (otherwise graphs will
% not be output), it will always output a table containing variables of
% interst.

% INPUTS
% trial - string for the name of the trial (e.g. Rmar15_01)
%  Plot - determines whether figures will be outputed, will only output if
%           set to 'Y' or 'Yes'

% OUTPUTS
% Table - table containing variables of interest with variable names as
%           column headings
% *Figures - (optional; require Plot = 'Y' or 'Yes')
%            first figure plots EMG, Motor data, and Ankle Angle through
%            time
%            second figure provides Fast Fourier Transform of raw &
%            filtered EMG data

%% Read in and Find Data of Interest
Igor = readtable([trial,'_Igor.csv']);
I = table2array(Igor);
% column 1 is EMG, 2 is Force, 3 is Length, and 4 is Trigger
Kin = readtable([trial,'_xypts.csv']);
K = table2array(Kin);
% [x-,y-Pt1, x-,y-Pt2, x-,y-Pt3, x-,y-Pt4]
% looking for angle between AB and CD where A = Pt1, B = Pt2, C = Pt3, & D
% = Pt4

%Trim Igor data to 1.5 second length
I = I(1:15000,:); % recorded at 10,000Hz
%Trim Igor data to only pre-trigger
while I(end,end)>2
    I(end,:) = [];
end
while I(end,end)<1
    I(end,:) = [];
end
while length(I)>10000
    I(1,:) = [];
end
%Determine number of frames to trim off front of video (if necessary)
n = round((1-length(I)/10000)*(length(K))); %ASSUMES filmed for 1 second!
if n>0
    for i = 1:n
        K(2,:) = [];
    end
end
    
v1 = [(K(1,1)-K(1,3)),(K(1,2)-K(1,4))];
v2 = [(K(:,7)-K(:,5)),(K(:,8)-K(:,6))];
u1 = v1/sqrt(dot(v1,v1));
for i = 1:length(K)
    u2(i,:) = v2(i,:)/sqrt(dot(v2(i,:),v2(i,:)));
    theta(i,1) = acosd(dot(u1,u2(i,:)));
end
fc = 200; %cutoff frequency
fs = 1964.64; %sample frequency in Hz
Wn = fc/(fs/2); % normalized cutt-off frequency (fc/nyguist = fc/(fs/2))
Order = 6;
[b,a] = butter(Order,Wn,'low'); % high-pass butterworth filter
theta_f = filtfilt(b,a,theta);
thresh = theta_f(1,1) - (theta_f(1,1)-min(theta_f))/100; % threshold for when movement is initiated (move 1% of max displacement)

% Correcting Force to be in Newtons & Motor Length to mm
I(:,2) = 2*I(:,2);
I(:,3) = 0.5*I(:,3);
% Creating times for each data point
t = linspace(0.0001,(length(I)/10000),length(I));
t_k = linspace(0,(length(I)/10000),length(K));

% Determining WHEN flexion movement starts (determined by threshold)
idx_flex = nnz(theta_f > thresh);
StartMove = t_k(1,idx_flex);

% Determine when start of extension first occurs and ends
theta_f_Dot = diff(theta_f);
theta_f_Dot(1:(idx_flex-1),1) = zeros((idx_flex-1),1); % sets initial velocity to zero
sign_theta_Dot = sign(theta_f_Dot); % sign of velocity
idx_ext = 1;
while sign_theta_Dot(1,1) < 1
    idx_ext = idx_ext+1; % index when ankle extension begins to occur
    sign_theta_Dot(1,:) = [];
end
idx_end = idx_ext; % creating index when ankle extension stops (ends)
while sign_theta_Dot(1,1) > 0
    idx_end = idx_end+1; % correcting index when ankle extension stops (ends)
    sign_theta_Dot(1,:) = [];
end

%% FILTERING Igor DATA
% Filtering EMG with High-pass Buttersworth
fc = 500; %cutoff frequency
fs = 10000; %sample frequency in Hz
Wn = fc/(fs/2); % normalized cutt-off frequency (fc/nyguist = fc/(fs/2))
Order = 2;
[b,a] = butter(Order,Wn,'high'); % high-pass butterworth filter
EMG_filt = filtfilt(b,a,I(:,1));

EMG_r = abs(EMG_filt(:,1)); % Absolute Value of Filtered EMG
% EMG_f(1:((StartMove-0.002)*10000),1) = abs(I(1:((StartMove-0.002)*10000),1)); %overwrites filter prior to movement with raw data (filter adds signal where there was none)

% Filtering Other Data
Force_f = lowpass(I(:,2),300,10000);
Length_f = lowpass(I(:,3),300,10000);

%% Analysis of EMG Delays and Intensity

% Quantifying background noise in EMG signal
Idx_Flex = round(StartMove*10000); % row index of EMG signal where angular displacement is triggered
Idx_Ext = round(t_k(1,idx_ext)*10000); % row index of EMG signal where angular displacement is triggered
Idx_End = round(t_k(1,idx_end)*10000); % row index of EMG signal where angular displacement is triggered
BkgndNoise = EMG_r((Idx_Flex-1050):(Idx_Flex-50),1); % Background Noise for rectified EMG signal from 0.1second interval starting 0.105seconds prior to movement intitiation
AvgNoise = mean(BkgndNoise); % average background noise during 0.1s interval
SDnoise = std(BkgndNoise);   % standard deviation of background noise during 0.1s interval
MaxNoise = max(BkgndNoise);  % maximum background noise during 0.1s interval

% Normalize rectified EMG signal to magnitude of the average rectified noise
AvgRectNoise = mean(EMG_r((Idx_Flex-1050):(Idx_Flex-50),1)); % average rectified background noise during 0.1s interval
EMG_N = EMG_r/AvgRectNoise;

% Measure delay between (start?) of sensory signals and motor response
%   Assumes sensory signal occurs while ankle is flexing & motor response (if
%       any) occurs during extension of ankle
%   Activation onset is defined as an EMG activity greater than 2 SD above the mean
%       background noise [EMG on > mean(Noise) + 2*std(Noise)]
EMG_thresh = (AvgNoise + 5*SDnoise); % 5SD captures 99.9999% of data, so expected value of false positive in 0.1s interval is below 1
[On_motor,EMG_s] = EMGOnsetTime(t,EMG_r,[Idx_Flex,Idx_End],EMG_thresh,5);
MotorDelay = On_motor - t(Idx_Flex);

% Defining EMG segments of Flexion & Extension
EMG_Flex = EMG_r(Idx_Flex:(Idx_Ext-1));
EMG_Ext = EMG_r(Idx_Ext:Idx_End);

% Calculating proportion of time activity above threshold during flexion &
% extension phases
t_on_flex = sum(EMG_Flex >= EMG_thresh)/length(EMG_Flex);
t_on_ext = sum(EMG_Ext >= EMG_thresh)/length(EMG_Ext);
Ratio_t_on_ext_flex = t_on_ext/t_on_flex; % ratio of the proportion of time muscle active during extension phase over proportion of time muscle active during flexion phase 

% Calculating intensity (area under rectified EMG curve) for both sensory &
% motor activity separately
I_flexion = sum(EMG_Flex(EMG_Flex >= EMG_thresh))*0.0001; % area under rectified EMG signal during flexion
I_extension = sum(EMG_Ext(EMG_Ext >= EMG_thresh))*0.0001; % area under rectified EMG signal during extension
Ratio_I_ext_flex = I_extension/I_flexion; % ratio of EMG intensity during extension over intensity during flexion

%% Collecting variables of interest into table
% Individual Number & Condition
if strcmp('Empty',trial(1:5))==1   % Empty Setup trials
    Individual = NaN;
    Condition = {'Empty Setup'};
    
else
    Individual = str2double(trial(5:6));
    % Condition
    if strcmp('d',trial(end-1))==1 % post-mortem trials
        Condition = {'Post-mortem'};
        
    elseif length(trial)==9
        Condition = {'Pre-op'};
        
    elseif strcmp('sham',trial(8:11))==1 % sham surgery
        Condition = {'Sham'};
        
    elseif strcmp('1wk',trial(8:10))==1 % 1wk post-op
        Condition = {'1wk Post-op'};
        
    elseif strcmp('3mo',trial(8:10))==1 % 3mo post-op
        Condition = {'3mo Post-op'};
        
    elseif strcmp('6mo',trial(8:10))==1 % 6mo post-op
        Condition = {'6mo Post-op'};
        
    else % Likely a problem
        Condition = {''};
        error([trial ' name is not structured correctly!']);
        
    end
end


% Trial Number
Trial = str2double(trial(end));

Table = table(Individual,Condition,Trial,EMG_thresh,I_flexion,I_extension,Ratio_I_ext_flex,t_on_flex,t_on_ext,Ratio_t_on_ext_flex);

%% Plotting (only happens when Plot = 'Y' or 'Yes'
if strcmp(Plot,'Y')==1 || strcmp(Plot,'Yes')==1
    % Plot Dataset to Visualize & Interpret Filters
    M_idx = t == On_motor;        
    
    figure('Name',['Stretch Reflex (' trial ')'])
    subplot(3,1,1)
    plot(t,abs(I(:,1)),'Color',[0.8,0.8,0.8]);
    hold on
    plot(t,EMG_r,'Color',[0,0,0.8]);
    plot(t,EMG_s,'r:');
    xline(t(Idx_Flex),'k:');
    xline(t(Idx_Ext),'k:');
    xline(t(Idx_End),'k:');
    yline((EMG_thresh),'k--','threshold','LabelHorizontalAlignment','right');
    if isnan(On_motor)==0
        plot(On_motor,EMG_s(M_idx),'ro');
    end
    hold off
    legend('Raw','Filt','','','','','Mot. Act','Location','northwest')
    % ylim([0, 0.8])
    % yline((AvgNoise+4*SDnoise),'--','Color',[0.5,0.5,0.5]);
    % xline(StartMove,'--','Color',[.5,.5,.5]);
    ylabel('Plant Act (V)') %y-axis label
    xlim([(StartMove-0.1), (StartMove+0.15)]) % only plots data from 0.1sec prior to movement to 1.5sec after start of movement
    
    subplot(3,1,2)
    yyaxis left
    plot(t,I(:,2),'Color',[0.8,0,0])
    xlim([(StartMove-0.1), (StartMove+0.15)])
    ylabel('Force (N)') %left y-axis label
    yyaxis right
    plot(t,I(:,3),'Color',[0.7,0.4,0])
    ylabel('Length motor (mm)')%right y-axis label
    
    
    subplot(3,1,3)
    plot(t_k,theta,'Color',[0,0.8,0])
    hold on
    xline(t(Idx_Flex),'k:');
    xline(t(Idx_Ext),'k:');
    xline(t(Idx_End),'k:');
    % ylim([40, 100]);
    xlim([(StartMove-0.1), (StartMove+0.15)])
    ylabel('Joint Ang (\circ)') %y-axis label
    xlabel('time (s)') %x-axis label
    hold off
    
    %% Plotting Fast Fourier Transform of EMG
    FFT_raw = abs(fft(I(:,1)));
    FFT_filt = abs(fft(EMG_filt(:,1)));
    t_fft = (0:(length(I(:,1))-1))/fs;
    f = t_fft*fs;
    
    figure('Name','Fourier Transform')
    subplot(2,1,1)
    plot(t_fft,I(:,1),'Color',[0.8,0,0])
    hold on
    plot(t_fft,EMG_filt(:,1),'Color',[0,0,0.8])
    hold off
    legend('Raw','Filtered','Location','northwest')
    ylabel('EMG (V)')
    xlabel('time (s)')
    
    subplot(2,1,2)
    plot(f,FFT_raw,'Color',[0.8,0,0])
    hold on
    plot(f,FFT_filt,'Color',[0,0,0.8])
    hold off
    legend('Raw','Filtered','Location','northwest')
    xlim([0,1000]) % limits x-axis to 0-1000Hz
    xlabel('frequency (Hz)')
    
    %% FFT for both flexion and extension phases
    Flex_fft = abs(fft(I(Idx_Flex:Idx_Ext,1)));
    Ext_fft = abs(fft(I(Idx_Ext:Idx_End,1)));
    
    t_Ffft = (0:(length(I(Idx_Flex:Idx_Ext))-1))/fs;
    fF = t_Ffft*fs;
    t_Efft = (0:(length(I(Idx_Ext:Idx_End))-1))/fs;
    fE = t_Efft*fs;
    
    
    figure('Name','Flexion & Extension FFT')
    subplot(2,1,1)
    plot(t,I(:,1),'k-')
    hold on
    plot(t(Idx_Flex:Idx_Ext),I(Idx_Flex:Idx_Ext,1),'b-');
    plot(t(Idx_Ext:Idx_End),I(Idx_Ext:Idx_End,1),'r-');
    ylabel('EMG (V)')
    yyaxis right
    plot(t_k,theta_f,'g')
    hold off
    legend('raw EMG','Flexion (sensory)','Extension (motor)','ankle anlge','Location','best');
    xlabel('time (s)')
    
    subplot(2,1,2)
    plot(fF,Flex_fft,'b-')
    ylabel('Flexion')
    hold on
    yyaxis right
    plot(fE,Ext_fft,'r-')
    ylabel('Extension')
    hold off
    xlabel('frequency (Hz)')
end
end