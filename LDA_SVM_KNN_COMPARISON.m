%% first filter
[s, h] = sload('A09T.gdf', 0, 'OVERFLOWDETECTION:OFF');

s=s';
sampleRate = 250; % Hz
lowEnd = 3; % Hz
highEnd = 30; % Hz
filterOrder = 2; % Filter order 
[s_1, s_2] = butter(filterOrder, [lowEnd highEnd]/(sampleRate/2)); % Generate filter coefficients
s_f = filtfilt(s_1, s_2, s); % Apply filter to data using zero-phase filtering
s_f=s_f';

position=transpose(h.EVENT.POS);
labels=transpose(h.EVENT.TYP);

%% EOG CORRECTION
noise=s_f(1:75000,23:25); % EOG channels during the 5 minutes given for EOG artifact detection before the trials
swoutEOG=s_f(1:75000,1:22);% Noisy signal without the EOG channels during the 5 minutes given for EOG artifact detection before the trial
autocov=inv((transpose(noise))*noise);% Inverse auto-covariance matrix of the EOG channels
crosscov=(transpose(noise))*swoutEOG;% cross-covariance between the EOG channels and EEG channels
b=autocov*crosscov;% weighting coefficients
swoutEOG2=s_f(:,1:22);% Complete noisy signal
noise2=s_f(:,23:25);% Complete EOG signal
S=swoutEOG2-noise2*b;
signal=transpose(S);%Clean Signal

%% band pass filter

lowEnd = 8; % Hz
highEnd = 13; % Hz
filterOrder = 2; % Filter order 
[Lower_limit, Upper_limit] = butter(filterOrder, [lowEnd highEnd]/(sampleRate/2)); % Generate filter coefficients
band_pass_filtered_data = filtfilt(Lower_limit, Upper_limit, signal); % Apply filter to data using zero-phase filtering

%% REJECTED TRIALS

index_reject=find(labels==1023);% Variable storing the indexes of the event types marked as rejected trials
index_reject=index_reject+1;% Moving the indexes one cell forward to be positioned in the indexes corresponding to the rejected tasks
labels(index_reject)=1023;% Replace the number stored in that index position for another so we no longer have one corresponding to a task
%% EPOCHING
epochs = zeros(size(band_pass_filtered_data,1),1250,length(labels));% Creates the empty array with the size necessary to store the epochs
for i=1:length(labels)% This for loop will create the epochs based on the positions of the events stored in the indexes of the "position" variable
epochs(:,:,i) = band_pass_filtered_data(:,position(i):position(i)+1249);
end

%% SEPARATE CLASSES
index1=find(labels==769);% First we find the indexes of only the left motor imagery taks
LEFTepochs=epochs(:,:,index1);% Now we extract and store the epochs of those specific indexes from the "epochs" array

index2=find(labels==770);% First we find the indexes of only the right motor imagery taks
RIGHTepochs=epochs(:,:,index2);% Now we extract and store the epochs of those specific indexes from the "epochs" array

index3=find(labels==771);% First we find the indexes of only the right motor imagery taks
FOOTepochs=epochs(:,:,index3);% Now we extract and store the epochs of those specific indexes from the "epochs" array

index4=find(labels==772);% First we find the indexes of only the right motor imagery taks
TONGUEepochs=epochs(:,:,index4);% Now we extract and store the epochs of those specific indexes from the "epochs" array
%% SPATIAL FILTER

load EEG_LEFT
LEFTepochs_spatially_filtered=zeros(size(LEFTepochs));
for i=1:size(LEFTepochs,3)
EEG.data=LEFTepochs(:,:,i);
LEFTepochs_spatially_filtered(:,:,i) = eeg_laplac(EEG, 1);
end

left_mean_spatial=mean(LEFTepochs_spatially_filtered,3);
left_mean_spatial=left_mean_spatial';


RIGHTepochs_spatially_filtered=zeros(size(RIGHTepochs));
for i=1:size(RIGHTepochs,3)
EEG.data=RIGHTepochs(:,:,i);
RIGHTepochs_spatially_filtered(:,:,i) = eeg_laplac(EEG, 1);
end

right_mean_spatial=mean(RIGHTepochs_spatially_filtered,3);
right_mean_spatial=right_mean_spatial';




FOOTepochs_spatially_filtered=zeros(size(FOOTepochs));
for i=1:size(FOOTepochs,3)
EEG.data=FOOTepochs(:,:,i);
FOOTepochs_spatially_filtered(:,:,i) = eeg_laplac(EEG, 1);
end

foot_mean_spatial=mean(FOOTepochs_spatially_filtered,3);
foot_mean_spatial=foot_mean_spatial';

TONGUEepochs_spatially_filtered=zeros(size(TONGUEepochs));
for i=1:size(TONGUEepochs,3)
EEG.data=TONGUEepochs(:,:,i);
TONGUEepochs_spatially_filtered(:,:,i) = eeg_laplac(EEG, 1);
end

tongue_mean_spatial=mean(TONGUEepochs_spatially_filtered,3);
tongue_mean_spatial=tongue_mean_spatial';
%% PSD

F_LEFT_EPOCHS=zeros(751,size(LEFTepochs_spatially_filtered,1),size(LEFTepochs_spatially_filtered,3));
for i=1:size(F_LEFT_EPOCHS,3)
    F_LEFT_EPOCHS(:,:,i)= LEFTepochs_spatially_filtered(:,250:1000,i)';
end
nfft = 1024; 
for i=1:size(F_LEFT_EPOCHS,3)
    for n=1:size(F_LEFT_EPOCHS,2)
[Lpxx(:,n,i),f_l] = pwelch(detrend(F_LEFT_EPOCHS(:,n,i)),200,50,nfft,250);
    end
end



F_RIGHT_EPOCHS=zeros(751,size(RIGHTepochs_spatially_filtered,1),size(RIGHTepochs_spatially_filtered,3));
for i=1:size(F_RIGHT_EPOCHS,3)
    F_RIGHT_EPOCHS(:,:,i)= RIGHTepochs_spatially_filtered(:,250:1000,i)';
end
nfft = 1024; 
for i=1:size(F_RIGHT_EPOCHS,3)
    for n=1:size(F_RIGHT_EPOCHS,2)
[Rpxx(:,n,i),f_r] = pwelch(detrend(F_RIGHT_EPOCHS(:,n,i)),200,50,nfft,250);
    end
end




F_FOOT_EPOCHS=zeros(751,size(FOOTepochs_spatially_filtered,1),size(FOOTepochs_spatially_filtered,3));
for i=1:size(F_FOOT_EPOCHS,3)
    F_FOOT_EPOCHS(:,:,i)= FOOTepochs_spatially_filtered(:,250:1000,i)';
end
nfft = 1024; 
for i=1:size(F_FOOT_EPOCHS,3)
    for n=1:size(F_FOOT_EPOCHS,2)
[Fpxx(:,n,i),f_f] = pwelch(detrend(F_FOOT_EPOCHS(:,n,i)),200,50,nfft,250);
    end
end

F_TONGUE_EPOCHS=zeros(751,size(TONGUEepochs_spatially_filtered,1),size(TONGUEepochs_spatially_filtered,3));
for i=1:size(F_TONGUE_EPOCHS,3)
    F_TONGUE_EPOCHS(:,:,i)= TONGUEepochs_spatially_filtered(:,250:1000,i)';
end
nfft = 1024; 
for i=1:size(F_TONGUE_EPOCHS,3)
    for n=1:size(F_TONGUE_EPOCHS,2)
[Tpxx(:,n,i),f_t] = pwelch(detrend(F_TONGUE_EPOCHS(:,n,i)),200,50,nfft,250);
    end
end
%% feature extraction


alpha_idx_l = find((f_l<=13)&(f_l>=8));% frequency index of alpha band power
alpha_idx_r = find((f_r<=13)&(f_r>=8));
alpha_idx_f = find((f_f<=13)&(f_f>=8));
alpha_idx_t = find((f_t<=13)&(f_t>=8));



left_alpha(:,1,:)=Lpxx(alpha_idx_l,19,:);
right_alpha(:,1,:)=Rpxx(alpha_idx_r,19,:);
foot_alpha(:,1,:)=Fpxx(alpha_idx_f,19,:);
tongue_alpha(:,1,:)=Tpxx(alpha_idx_t,19,:);

left_alpha(:,2,:)=Lpxx(alpha_idx_l,11,:);
right_alpha(:,2,:)=Rpxx(alpha_idx_r,11,:);
foot_alpha(:,2,:)=Fpxx(alpha_idx_f,11,:);
tongue_alpha(:,2,:)=Tpxx(alpha_idx_t,11,:);



left_class=zeros((size(left_alpha,1)*size(left_alpha,2)),size(left_alpha,3));
x=1;
y=size(left_alpha,1);
for i=1:size(left_class,2)
   for j=1:size(left_alpha,2)
       
          left_class(x:y,i)=left_alpha(:,j,i);
          x=x+size(left_alpha,1);
          y=y+size(left_alpha,1);
   end
    x=1;
    y=size(left_alpha,1);
end
left_class=left_class';


right_class=zeros((size(right_alpha,1)*size(right_alpha,2)),size(right_alpha,3));
x=1;
y=size(right_alpha,1);
for i=1:size(right_class,2)
   for j=1:size(right_alpha,2)
       
          right_class(x:y,i)=right_alpha(:,j,i);
          x=x+size(right_alpha,1);
          y=y+size(right_alpha,1);
   end
    x=1;
    y=size(right_alpha,1);
end

right_class=right_class';


foot_class=zeros((size(foot_alpha,1)*size(foot_alpha,2)),size(foot_alpha,3));
x=1;
y=size(foot_alpha,1);
for i=1:size(foot_class,2)
   for j=1:size(foot_alpha,2)
       
          foot_class(x:y,i)=foot_alpha(:,j,i);
          x=x+size(foot_alpha,1);
          y=y+size(foot_alpha,1);
   end
    x=1;
    y=size(foot_alpha,1);
end

foot_class=foot_class';

tongue_class=zeros((size(tongue_alpha,1)*size(tongue_alpha,2)),size(tongue_alpha,3));
x=1;
y=size(tongue_alpha,1);
for i=1:size(tongue_class,2)
   for j=1:size(tongue_alpha,2)
       
          tongue_class(x:y,i)=tongue_alpha(:,j,i);
          x=x+size(tongue_alpha,1);
          y=y+size(tongue_alpha,1);
   end
    x=1;
    y=size(tongue_alpha,1);
end

tongue_class=tongue_class';
%% TOPOPLOTS WITH CHANNELS SIGNALS
lmt_p2=left_mean_spatial(:,19)';
lmt_fc3=left_mean_spatial(:,11)';



left_alpha_mean=mean(Lpxx,3);
right_alpha_mean=mean(Rpxx,3);
foot_alpha_mean=mean(Fpxx,3);
tongue_alpha_mean=mean(Tpxx,3);

lim_neg_caxis_left=-10;
lim_positive_caxis_left=10;
ylim_neg_left=-15;
ylim_pos_left=15;
ylim_neg_psd_left=-20;
ylim_pos_psd_left=20;


lim_neg_caxis_right=-5;
lim_positive_caxis_right=5;
ylim_neg_right=-15;
ylim_pos_right=15;
ylim_neg_psd_right=-20;
ylim_pos_psd_right=20;


lim_neg_caxis_foot=-5;
lim_positive_caxis_foot=5;
ylim_neg_foot=-15;
ylim_pos_foot=15;
ylim_neg_psd_foot=-20;
ylim_pos_psd_foot=20;

lim_neg_caxis_tongue=-5;
lim_positive_caxis_tongue=5;
ylim_neg_tongue=-15;
ylim_pos_tongue=15;
ylim_neg_psd_tongue=-20;
ylim_pos_psd_tongue=20;

figure()


subplot(3,7,1)
topoplot(left_mean_spatial(250,:), 'miscanales.ced');colorbar;
caxis([lim_neg_caxis_left lim_positive_caxis_left])
title('Elec. reading at 1 S')

subplot(3,7,2)
topoplot(left_mean_spatial(375,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_left lim_positive_caxis_left])
title('Elec. reading at 1.5 S')

subplot(3,7,3)
topoplot(left_mean_spatial(500,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_left lim_positive_caxis_left])
title('Elec. reading at 2 S')

subplot(3,7,4)
topoplot(left_mean_spatial(625,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_left lim_positive_caxis_left])
title('Elec. reading at 2.5 S')

subplot(3,7,5)
topoplot(left_mean_spatial(750,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_left lim_positive_caxis_left])
title('Elec. reading at 3 S')

subplot(3,7,6)
topoplot(left_mean_spatial(875,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_left lim_positive_caxis_left])
title('Elec. reading at 3.5 S')

subplot(3,7,7)
topoplot(left_mean_spatial(1000,:), 'miscanales.ced');colorbar;
caxis([lim_neg_caxis_left lim_positive_caxis_left])
title('Elec. reading at 4 S')

subplot(3,7,8:14)
p2=plot(lmt_p2(1,250:1000));M1 = "Channel P1";
xlabel('Time')
ylabel('Amplitud')
set(gca,'XTick',0:125:750)
set(gca,'XTickLabel',1:0.5:4)
ylim([ylim_neg_left ylim_pos_left])
xlim([0 750])
grid on

hold on
fc3=plot(lmt_fc3(1,250:1000));M2 = "Channel C2";
xlabel('Time')
ylabel('Amplitud')
set(gca,'XTick',0:125:750)
set(gca,'XTickLabel',1:0.5:4)
ylim([ylim_neg_left ylim_pos_left])
xlim([0 750])
grid on

title('Electrical activity in channels P1, C2')
legend([p2,fc3], [M1, M2]);



subplot(3,7,15:21)
psdc2=plot(f_l,10*log10(left_alpha_mean(:,19)));
hold on
psdfc3=plot(f_l,10*log10(left_alpha_mean(:,11)));
hold on



xlabel('Frequency (Hz)')
ylabel('dB/Hz')
ylim([ylim_neg_psd_left ylim_pos_psd_left])
xlim([8 13])
title('Left task PSD')
legend([psdc2,psdfc3], [M1, M2]);
grid on
sgtitle('LEFT TASK')




rmt_p2=right_mean_spatial(:,19)';
rmt_fc3=right_mean_spatial(:,11)';





figure()


subplot(3,7,1)
topoplot(right_mean_spatial(250,:), 'miscanales.ced');colorbar;
caxis([lim_neg_caxis_right lim_positive_caxis_right])
title('Elec. reading at 1 S')

subplot(3,7,2)
topoplot(right_mean_spatial(375,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_right lim_positive_caxis_right])
title('Elec. reading at 1.5 S')

subplot(3,7,3)
topoplot(right_mean_spatial(500,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_right lim_positive_caxis_right])
title('Elec. reading at 2 S')

subplot(3,7,4)
topoplot(right_mean_spatial(625,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_right lim_positive_caxis_right])
title('Elec. reading at 2.5 S')

subplot(3,7,5)
topoplot(right_mean_spatial(750,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_right lim_positive_caxis_right])
title('Elec. reading at 3 S')

subplot(3,7,6)
topoplot(right_mean_spatial(875,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_right lim_positive_caxis_right])
title('Elec. reading at 3.5 S')

subplot(3,7,7)
topoplot(right_mean_spatial(1000,:), 'miscanales.ced');colorbar;
caxis([lim_neg_caxis_right lim_positive_caxis_right])
title('Elec. reading at 4 S')

subplot(3,7,8:14)
rp2=plot(rmt_p2(1,250:1000));rM1 = "Channel P1";
xlabel('Time')
ylabel('Amplitud')
set(gca,'XTick',0:125:750)
set(gca,'XTickLabel',1:0.5:4)
ylim([ylim_neg_right ylim_pos_right])
xlim([0 750])
grid on

hold on
rfc3=plot(rmt_fc3(1,250:1000));rM2 = "Channel C2";
xlabel('Time')
ylabel('Amplitud')
set(gca,'XTick',0:125:750)
set(gca,'XTickLabel',1:0.5:4)
ylim([ylim_neg_right ylim_pos_right])
xlim([0 750])
grid on





title('Electrical activity in channels P1, C2')
legend([rp2,rfc3], [rM1, rM2]);



subplot(3,7,15:21)
rpsdc2=plot(f_r,10*log10(right_alpha_mean(:,19)));
hold on
rpsdfc3=plot(f_r,10*log10(right_alpha_mean(:,11)));
hold on







xlabel('Frequency (Hz)')
ylabel('dB/Hz')
ylim([ylim_neg_psd_right ylim_pos_psd_right])
xlim([8 13])
title('Right task PSD')
legend([rpsdc2,rpsdfc3], [rM1, rM2]);
grid on
sgtitle('RIGHT TASK')

fmt_p2=foot_mean_spatial(:,19)';
fmt_fc3=foot_mean_spatial(:,11)';




figure()


subplot(3,7,1)
topoplot(foot_mean_spatial(250,:), 'miscanales.ced');colorbar;
caxis([lim_neg_caxis_foot lim_positive_caxis_foot])
title('Elec. reading at 1 S')

subplot(3,7,2)
topoplot(foot_mean_spatial(375,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_foot lim_positive_caxis_foot])
title('Elec. reading at 1.5 S')

subplot(3,7,3)
topoplot(foot_mean_spatial(500,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_foot lim_positive_caxis_foot])
title('Elec. reading at 2 S')

subplot(3,7,4)
topoplot(foot_mean_spatial(625,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_foot lim_positive_caxis_foot])
title('Elec. reading at 2.5 S')

subplot(3,7,5)
topoplot(foot_mean_spatial(750,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_foot lim_positive_caxis_foot])
title('Elec. reading at 3 S')

subplot(3,7,6)
topoplot(foot_mean_spatial(875,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_foot lim_positive_caxis_foot])
title('Elec. reading at 3.5 S')

subplot(3,7,7)
topoplot(foot_mean_spatial(1000,:), 'miscanales.ced');colorbar;
caxis([lim_neg_caxis_foot lim_positive_caxis_foot])
title('Elec. reading at 4 S')

subplot(3,7,8:14)
fp2=plot(fmt_p2(1,250:1000));fM1 = "Channel P1";
xlabel('Time')
ylabel('Amplitud')
set(gca,'XTick',0:125:750)
set(gca,'XTickLabel',1:0.5:4)
ylim([ylim_neg_foot ylim_pos_foot])
xlim([0 750])
grid on

hold on
ffc3=plot(fmt_fc3(1,250:1000));fM2 = "Channel C2";
xlabel('Time')
ylabel('Amplitud')
set(gca,'XTick',0:125:750)
set(gca,'XTickLabel',1:0.5:4)
ylim([ylim_neg_foot ylim_pos_foot])
xlim([0 750])
grid on




title('Electrical activity in channels P1, C2')
legend([fp2,ffc3], [fM1, fM2]);



subplot(3,7,15:21)
fpsdc2=plot(f_f,10*log10(foot_alpha_mean(:,19)));
hold on

fpsdfc3=plot(f_f,10*log10(foot_alpha_mean(:,11)));
hold on






xlabel('Frequency (Hz)')
ylabel('dB/Hz')
ylim([ylim_neg_psd_foot ylim_pos_psd_foot])
xlim([8 13])
title('Foot task PSD')
legend([fpsdc2,fpsdfc3], [fM1, fM2]);
grid on
sgtitle('FOOT TASK')


tmt_p2=tongue_mean_spatial(:,19)';
tmt_fc3=tongue_mean_spatial(:,11)';


figure()


subplot(3,7,1)
topoplot(tongue_mean_spatial(250,:), 'miscanales.ced');colorbar;
caxis([lim_neg_caxis_tongue lim_positive_caxis_tongue])
title('Elec. reading at 1 S')

subplot(3,7,2)
topoplot(tongue_mean_spatial(375,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_tongue lim_positive_caxis_tongue])
title('Elec. reading at 1.5 S')

subplot(3,7,3)
topoplot(tongue_mean_spatial(500,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_tongue lim_positive_caxis_tongue])
title('Elec. reading at 2 S')

subplot(3,7,4)
topoplot(tongue_mean_spatial(625,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_tongue lim_positive_caxis_tongue])
title('Elec. reading at 2.5 S')

subplot(3,7,5)
topoplot(tongue_mean_spatial(750,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_tongue lim_positive_caxis_tongue])
title('Elec. reading at 3 S')

subplot(3,7,6)
topoplot(tongue_mean_spatial(875,:), 'miscanales.ced');colorbar; 
caxis([lim_neg_caxis_tongue lim_positive_caxis_tongue])
title('Elec. reading at 3.5 S')

subplot(3,7,7)
topoplot(tongue_mean_spatial(1000,:), 'miscanales.ced');colorbar;
caxis([lim_neg_caxis_tongue lim_positive_caxis_tongue])
title('Elec. reading at 4 S')

subplot(3,7,8:14)
tp2=plot(tmt_p2(1,250:1000));tM1 = "Channel P1";
xlabel('Time')
ylabel('Amplitud')
set(gca,'XTick',0:125:750)
set(gca,'XTickLabel',1:0.5:4)
ylim([ylim_neg_tongue ylim_pos_tongue])
xlim([0 750])
grid on

hold on
tfc3=plot(tmt_fc3(1,250:1000));tM2 = "Channel C2";
xlabel('Time')
ylabel('Amplitud')
set(gca,'XTick',0:125:750)
set(gca,'XTickLabel',1:0.5:4)
ylim([ylim_neg_tongue ylim_pos_tongue])
xlim([0 750])
grid on




title('Electrical activity in channels P1, C2')
legend([tp2,tfc3], [tM1, tM2]);



subplot(3,7,15:21)
tpsdc2=plot(f_t,10*log10(tongue_alpha_mean(:,19)));
hold on
tpsdfc3=plot(f_t,10*log10(tongue_alpha_mean(:,11)));
hold on







xlabel('Frequency (Hz)')
ylabel('dB/Hz')
ylim([ylim_neg_psd_tongue ylim_pos_psd_tongue])
xlim([8 13])
title('Tongue task PSD')
legend([tpsdc2,tpsdfc3], [tM1, tM2]);
grid on
sgtitle('TONGUE TASK')

%% Classifiers


all_samples=[left_class;right_class;foot_class;tongue_class];
all_labels=[1+zeros(size(left_class,1),1);2+zeros(size(right_class,1),1);3+zeros(size(foot_class,1),1);4+zeros(size(tongue_class,1),1)];


K= 10;
tp_left=zeros(size(K));
tn_left=zeros(size(K));
fp_left=zeros(size(K));
fn_left=zeros(size(K));
tp_right=zeros(size(K));
tn_right=zeros(size(K));
fp_right=zeros(size(K));
fn_right=zeros(size(K));
tp_foot=zeros(size(K));
tn_foot=zeros(size(K));
fp_foot=zeros(size(K));
fn_foot=zeros(size(K));
tp_tongue=zeros(size(K));
tn_tongue=zeros(size(K));
fp_tongue=zeros(size(K));
fn_tongue=zeros(size(K));

lfnr=zeros(size(K));
lfnf=zeros(size(K));
lfnt=zeros(size(K));

rfnl=zeros(size(K));
rfnf=zeros(size(K));
rfnt=zeros(size(K));

ffnl=zeros(size(K));
ffnr=zeros(size(K));
ffnt=zeros(size(K));

tfnl=zeros(size(K));
tfnr=zeros(size(K));
tfnf=zeros(size(K));






indices=crossvalind('Kfold',all_labels,K);
for k=1:K
    cv_test_idx = find(indices == k); % indices for test samples in one trial of validation
    cv_train_idx = find(indices ~= k); % indices for training samples in one trial of validation
    Mdl_lda = fitcdiscr(all_samples(cv_train_idx,:),all_labels(cv_train_idx),'DiscrimType','linear','Gamma',1);%
    cv_classout = predict(Mdl_lda,all_samples(cv_test_idx,:));
    cv_acc(k) = mean(cv_classout==all_labels(cv_test_idx));
    TP_left_class = sum((cv_classout==all_labels(cv_test_idx))&(cv_classout==1));
    TN_left_class = sum((cv_classout==all_labels(cv_test_idx))&(cv_classout~=1));
    FN_left_class = sum((cv_classout~=1)&(all_labels(cv_test_idx)==1));
    FP_left_class = sum((cv_classout==1)&(all_labels(cv_test_idx)~=1));
    
   
    LFNR = sum((cv_classout~=all_labels(cv_test_idx))&(cv_classout==1)&(all_labels(cv_test_idx)==2));
    LFNF = sum((cv_classout~=all_labels(cv_test_idx))&(cv_classout==1)&(all_labels(cv_test_idx)==3));
    LFNT = sum((cv_classout~=all_labels(cv_test_idx))&(cv_classout==1)&(all_labels(cv_test_idx)==4));
    

    TP_right_class = sum((cv_classout==all_labels(cv_test_idx))&(cv_classout==2));
    TN_right_class = sum((cv_classout==all_labels(cv_test_idx))&(cv_classout~=2));
    FN_right_class = sum((cv_classout~=2)&(all_labels(cv_test_idx)==2));
    FP_right_class = sum((cv_classout==2)&(all_labels(cv_test_idx)~=2));
    
    RFNL = sum((cv_classout~=all_labels(cv_test_idx))&(cv_classout==2)&(all_labels(cv_test_idx)==1));
    RFNF = sum((cv_classout~=all_labels(cv_test_idx))&(cv_classout==2)&(all_labels(cv_test_idx)==3));
    RFNT = sum((cv_classout~=all_labels(cv_test_idx))&(cv_classout==2)&(all_labels(cv_test_idx)==4));
    

    TP_foot_class = sum((cv_classout==all_labels(cv_test_idx))&(cv_classout==3));
    TN_foot_class = sum((cv_classout==all_labels(cv_test_idx))&(cv_classout~=3));
    FN_foot_class = sum((cv_classout~=3)&(all_labels(cv_test_idx)==3));
    FP_foot_class = sum((cv_classout==3)&(all_labels(cv_test_idx)~=3));
    
    FFNL = sum((cv_classout~=all_labels(cv_test_idx))&(cv_classout==3)&(all_labels(cv_test_idx)==1));
    FFNR = sum((cv_classout~=all_labels(cv_test_idx))&(cv_classout==3)&(all_labels(cv_test_idx)==2));
    FFNT = sum((cv_classout~=all_labels(cv_test_idx))&(cv_classout==3)&(all_labels(cv_test_idx)==4));


    TP_tongue_class = sum((cv_classout==all_labels(cv_test_idx))&(cv_classout==4));
    TN_tongue_class = sum((cv_classout==all_labels(cv_test_idx))&(cv_classout~=4));
    FN_tongue_class = sum((cv_classout~=4)&(all_labels(cv_test_idx)==4));
    FP_tongue_class = sum((cv_classout==4)&(all_labels(cv_test_idx)~=4));
    
    TFNL = sum((cv_classout~=all_labels(cv_test_idx))&(cv_classout==4)&(all_labels(cv_test_idx)==1));
    TFNR = sum((cv_classout~=all_labels(cv_test_idx))&(cv_classout==4)&(all_labels(cv_test_idx)==2));
    TFNF = sum((cv_classout~=all_labels(cv_test_idx))&(cv_classout==4)&(all_labels(cv_test_idx)==3));

    cv_sensitivity_left(k) = TP_left_class/(TP_left_class+FN_left_class);%Tells what fraction of predictions as left class were actually left
    cv_specificity_left(k) = TN_left_class/(TN_left_class+FP_left_class);%Tells what fraction of all not left samples are correctly predicted as not left by the classifier.

    cv_sensitivity_right(k) = TP_right_class/(TP_right_class+FN_right_class);%Tells what fraction of predictions as right class were actually right
    cv_specificity_right(k) = TN_right_class/(TN_right_class+FP_right_class);%Tells what fraction of all not right samples are correctly predicted as not right by the classifier.

    cv_sensitivity_foot(k) = TP_foot_class/(TP_foot_class+FN_foot_class);%Tells what fraction of predictions as foot class were actually foot
    cv_specificity_foot(k) = TN_foot_class/(TN_foot_class+FP_foot_class);%Tells what fraction of all not foot samples are correctly predicted as not foot by the classifier.

    cv_sensitivity_tongue(k) = TP_tongue_class/(TP_tongue_class+FN_tongue_class);%Tells what fraction of predictions as tongue class were actually tongue 
    cv_specificity_tongue(k) = TN_tongue_class/(TN_tongue_class+FP_tongue_class);%Tells what fraction of all not tongue samples are correctly predicted as not tongue by the classifier.
    tp_left(k)=TP_left_class;
    tn_left(k)=TN_left_class;
    fp_left(k)=FP_left_class;
    fn_left(k)=FN_left_class;
    
    lfnr(k)=LFNR;
    lfnf(k)=LFNF;
    lfnt(k)=LFNT;
    
    
    
    tp_right(k)=TP_right_class;
    tn_right(k)=TN_right_class;
    fp_right(k)=FP_right_class;
    fn_right(k)=FN_right_class;
    
    rfnl(k)=RFNL;
    rfnf(k)=RFNF;
    rfnt(k)=RFNT;
    
    tp_foot(k)=TP_foot_class;
    tn_foot(k)=TN_foot_class;
    fp_foot(k)=FP_foot_class;
    fn_foot(k)=FN_foot_class;
    
    ffnl(k)=FFNL;
    ffnr(k)=FFNR;
    ffnt(k)=FFNT;
    
    tp_tongue(k)=TP_tongue_class;
    tn_tongue(k)=TN_tongue_class;
    fp_tongue(k)=FP_tongue_class;
    fn_tongue(k)=FN_tongue_class;
    
    tfnl(k)=TFNL;
    tfnr(k)=TFNR;
    tfnf(k)=TFNF;

    

end
    TP_left_class = sum(tp_left);
    TN_left_class = sum(tn_left);
    FP_left_class = sum(fp_left);
    FN_left_class = sum(fn_left);
    
    LFNR = sum(lfnr);
    LFNF = sum(lfnf);
    LFNT = sum(lfnt);

    
    TP_right_class = sum(tp_right);
    TN_right_class = sum(tn_right);
    FP_right_class = sum(fp_right);
    FN_right_class = sum(fn_right);
    
    RFNL = sum(rfnl);
    RFNF = sum(rfnf);
    RFNT = sum(rfnt);
   
    
    TP_foot_class = sum(tp_foot);
    TN_foot_class = sum(tn_foot);
    FP_foot_class = sum(fp_foot);
    FN_foot_class = sum(fn_foot);
    
    FFNL = sum(ffnl);
    FFNR = sum(ffnr);
    FFNT = sum(ffnt);
    
    
    TP_tongue_class = sum(tp_tongue);
    TN_tongue_class = sum(tn_tongue);
    FP_tongue_class = sum(fp_tongue);
    FN_tongue_class = sum(fn_tongue);
    
    
    TFNL = sum(tfnl);
    TFNR = sum(tfnr);
    TFNF = sum(tfnf);
    
cv_acc_avg = mean(cv_acc); % averaged accuracy
cv_sensitivity_left_avg = mean(cv_sensitivity_left);  % averaged sensitivity
cv_specificity_left_avg = mean(cv_specificity_left);  % averaged specificity

cv_sensitivity_right_avg = mean(cv_sensitivity_right);  % averaged sensitivity
cv_specificity_right_avg = mean(cv_specificity_right);  % averaged specificity

cv_sensitivity_foot_avg = mean(cv_sensitivity_foot);  % averaged sensitivity
cv_specificity_foot_avg = mean(cv_specificity_foot);  % averaged specificity

cv_sensitivity_tongue_avg = mean(cv_sensitivity_tongue);  % averaged sensitivity
cv_specificity_tongue_avg = mean(cv_specificity_tongue);  % averaged specificity

LDA_CON_MATRIX= [TP_left_class  LFNR LFNF LFNT;  RFNL TP_right_class RFNF RFNT;  FFNL FFNR TP_foot_class FFNT; TFNL TFNR TFNF TP_tongue_class];
class_labels = ["Left Class";"Right Class";"Foot Class";"Tongue Class"];

figure()
confusionplot=confusionchart(LDA_CON_MATRIX,class_labels,'Title','LDA Confusion Matrix 4 Classes');
M= 10;
tp_left_svm=zeros(size(M));
tn_left_svm=zeros(size(M));
fp_left_svm=zeros(size(M));
fn_left_svm=zeros(size(M));
tp_right_svm=zeros(size(M));
tn_right_svm=zeros(size(M));
fp_right_svm=zeros(size(M));
fn_right_svm=zeros(size(M));
tp_foot_svm=zeros(size(M));
tn_foot_svm=zeros(size(M));
fp_foot_svm=zeros(size(M));
fn_foot_svm=zeros(size(M));
tp_tongue_svm=zeros(size(M));
tn_tongue_svm=zeros(size(M));
fp_tongue_svm=zeros(size(M));
fn_tongue_svm=zeros(size(M));

lfnr_svm=zeros(size(M));
lfnf_svm=zeros(size(M));
lfnt_svm=zeros(size(M));

rfnl_svm=zeros(size(M));
rfnf_svm=zeros(size(M));
rfnt_svm=zeros(size(M));

ffnl_svm=zeros(size(M));
ffnr_svm=zeros(size(M));
ffnt_svm=zeros(size(M));

tfnl_svm=zeros(size(M));
tfnr_svm=zeros(size(M));
tfnf_svm=zeros(size(M));
indices_svm=crossvalind('Kfold',all_labels,M);
for m=1:M
    cv_test_idx_svm = find(indices_svm == m); % indices for test samples in one trial of validation
    cv_train_idx_svm = find(indices_svm ~= m); % indices for training samples in one trial of validation
    Mdl = fitcecoc(all_samples(cv_train_idx_svm,:),all_labels(cv_train_idx_svm));
    cv_classout_svm = predict(Mdl,all_samples(cv_test_idx_svm,:));
    cv_acc_svm(m) = mean(cv_classout_svm==all_labels(cv_test_idx_svm));
    TP_left_class_svm = sum((cv_classout_svm==all_labels(cv_test_idx_svm))&(cv_classout_svm==1));
    TN_left_class_svm = sum((cv_classout_svm==all_labels(cv_test_idx_svm))&(cv_classout_svm~=1));
    FN_left_class_svm = sum((cv_classout_svm~=1)&(all_labels(cv_test_idx_svm)==1));
    FP_left_class_svm = sum((cv_classout_svm==1)&(all_labels(cv_test_idx_svm)~=1));
    
   
    LFNR_svm = sum((cv_classout_svm~=all_labels(cv_test_idx_svm))&(cv_classout_svm==1)&(all_labels(cv_test_idx_svm)==2));
    LFNF_svm = sum((cv_classout_svm~=all_labels(cv_test_idx_svm))&(cv_classout_svm==1)&(all_labels(cv_test_idx_svm)==3));
    LFNT_svm = sum((cv_classout_svm~=all_labels(cv_test_idx_svm))&(cv_classout_svm==1)&(all_labels(cv_test_idx_svm)==4));
    

    TP_right_class_svm = sum((cv_classout_svm==all_labels(cv_test_idx_svm))&(cv_classout_svm==2));
    TN_right_class_svm = sum((cv_classout_svm==all_labels(cv_test_idx_svm))&(cv_classout_svm~=2));
    FN_right_class_svm = sum((cv_classout_svm~=2)&(all_labels(cv_test_idx_svm)==2));
    FP_right_class_svm = sum((cv_classout_svm==2)&(all_labels(cv_test_idx_svm)~=2));
    
    RFNL_svm = sum((cv_classout_svm~=all_labels(cv_test_idx_svm))&(cv_classout_svm==2)&(all_labels(cv_test_idx_svm)==1));
    RFNF_svm = sum((cv_classout_svm~=all_labels(cv_test_idx_svm))&(cv_classout_svm==2)&(all_labels(cv_test_idx_svm)==3));
    RFNT_svm = sum((cv_classout_svm~=all_labels(cv_test_idx_svm))&(cv_classout_svm==2)&(all_labels(cv_test_idx_svm)==4));
    

    TP_foot_class_svm = sum((cv_classout_svm==all_labels(cv_test_idx_svm))&(cv_classout_svm==3));
    TN_foot_class_svm = sum((cv_classout_svm==all_labels(cv_test_idx_svm))&(cv_classout_svm~=3));
    FN_foot_class_svm = sum((cv_classout_svm~=3)&(all_labels(cv_test_idx_svm)==3));
    FP_foot_class_svm = sum((cv_classout_svm==3)&(all_labels(cv_test_idx_svm)~=3));
    
    FFNL_svm = sum((cv_classout_svm~=all_labels(cv_test_idx_svm))&(cv_classout_svm==3)&(all_labels(cv_test_idx_svm)==1));
    FFNR_svm = sum((cv_classout_svm~=all_labels(cv_test_idx_svm))&(cv_classout_svm==3)&(all_labels(cv_test_idx_svm)==2));
    FFNT_svm = sum((cv_classout_svm~=all_labels(cv_test_idx_svm))&(cv_classout_svm==3)&(all_labels(cv_test_idx_svm)==4));


    TP_tongue_class_svm = sum((cv_classout_svm==all_labels(cv_test_idx_svm))&(cv_classout_svm==4));
    TN_tongue_class_svm = sum((cv_classout_svm==all_labels(cv_test_idx_svm))&(cv_classout_svm~=4));
    FN_tongue_class_svm = sum((cv_classout_svm~=4)&(all_labels(cv_test_idx_svm)==4));
    FP_tongue_class_svm = sum((cv_classout_svm==4)&(all_labels(cv_test_idx_svm)~=4));
    
    TFNL_svm = sum((cv_classout_svm~=all_labels(cv_test_idx_svm))&(cv_classout_svm==4)&(all_labels(cv_test_idx_svm)==1));
    TFNR_svm = sum((cv_classout_svm~=all_labels(cv_test_idx_svm))&(cv_classout_svm==4)&(all_labels(cv_test_idx_svm)==2));
    TFNF_svm = sum((cv_classout_svm~=all_labels(cv_test_idx_svm))&(cv_classout_svm==4)&(all_labels(cv_test_idx_svm)==3));

    cv_sensitivity_left_svm(m) = TP_left_class_svm/(TP_left_class_svm+FN_left_class_svm);%Tells what fraction of predictions as left class were actually left
    cv_specificity_left_svm(m) = TN_left_class_svm/(TN_left_class_svm+FP_left_class_svm);%Tells what fraction of all not left samples are correctly predicted as not left by the classifier.

    cv_sensitivity_right_svm(m) = TP_right_class_svm/(TP_right_class_svm+FN_right_class_svm);%Tells what fraction of predictions as right class were actually right
    cv_specificity_right_svm(m) = TN_right_class_svm/(TN_right_class_svm+FP_right_class_svm);%Tells what fraction of all not right samples are correctly predicted as not right by the classifier.

    cv_sensitivity_foot_svm(m) = TP_foot_class_svm/(TP_foot_class_svm+FN_foot_class_svm);%Tells what fraction of predictions as foot class were actually foot
    cv_specificity_foot_svm(m) = TN_foot_class_svm/(TN_foot_class_svm+FP_foot_class_svm);%Tells what fraction of all not foot samples are correctly predicted as not foot by the classifier.

    cv_sensitivity_tongue_svm(m) = TP_tongue_class_svm/(TP_tongue_class_svm+FN_tongue_class_svm);%Tells what fraction of predictions as tongue class were actually tongue 
    cv_specificity_tongue_svm(m) = TN_tongue_class_svm/(TN_tongue_class_svm+FP_tongue_class_svm);%Tells what fraction of all not tongue samples are correctly predicted as not tongue by the classifier.
    tp_left_svm(m)=TP_left_class_svm;
    tn_left_svm(m)=TN_left_class_svm;
    fp_left_svm(m)=FP_left_class_svm;
    fn_left_svm(m)=FN_left_class_svm;
    
    lfnr_svm(m)=LFNR_svm;
    lfnf_svm(m)=LFNF_svm;
    lfnt_svm(m)=LFNT_svm;
    
    
    
    tp_right_svm(m)=TP_right_class_svm;
    tn_right_svm(m)=TN_right_class_svm;
    fp_right_svm(m)=FP_right_class_svm;
    fn_right_svm(m)=FN_right_class_svm;
    
    rfnl_svm(m)=RFNL_svm;
    rfnf_svm(m)=RFNF_svm;
    rfnt_svm(m)=RFNT_svm;
    
    tp_foot_svm(m)=TP_foot_class_svm;
    tn_foot_svm(m)=TN_foot_class_svm;
    fp_foot_svm(m)=FP_foot_class_svm;
    fn_foot_svm(m)=FN_foot_class_svm;
    
    ffnl_svm(m)=FFNL_svm;
    ffnr_svm(m)=FFNR_svm;
    ffnt_svm(m)=FFNT_svm;
    
    tp_tongue_svm(m)=TP_tongue_class_svm;
    tn_tongue_svm(m)=TN_tongue_class_svm;
    fp_tongue_svm(m)=FP_tongue_class_svm;
    fn_tongue_svm(m)=FN_tongue_class_svm;
    
    tfnl_svm(m)=TFNL_svm;
    tfnr_svm(m)=TFNR_svm;
    tfnf_svm(m)=TFNF_svm;

    
end
    TP_left_class_svm = sum(tp_left_svm);
    TN_left_class_svm = sum(tn_left_svm);
    FP_left_class_svm = sum(fp_left_svm);
    FN_left_class_svm = sum(fn_left_svm);
    
    LFNR_svm = sum(lfnr_svm);
    LFNF_svm = sum(lfnf_svm);
    LFNT_svm = sum(lfnt_svm);

    
    TP_right_class_svm = sum(tp_right_svm);
    TN_right_class_svm = sum(tn_right_svm);
    FP_right_class_svm = sum(fp_right_svm);
    FN_right_class_svm = sum(fn_right_svm);
    
    RFNL_svm = sum(rfnl_svm);
    RFNF_svm = sum(rfnf_svm);
    RFNT_svm = sum(rfnt_svm);
   
    
    TP_foot_class_svm = sum(tp_foot_svm);
    TN_foot_class_svm = sum(tn_foot_svm);
    FP_foot_class_svm = sum(fp_foot_svm);
    FN_foot_class_svm = sum(fn_foot_svm);
    
    FFNL_svm = sum(ffnl_svm);
    FFNR_svm = sum(ffnr_svm);
    FFNT_svm = sum(ffnt_svm);
    
    
    TP_tongue_class_svm = sum(tp_tongue_svm);
    TN_tongue_class_svm = sum(tn_tongue_svm);
    FP_tongue_class_svm = sum(fp_tongue_svm);
    FN_tongue_class_svm = sum(fn_tongue_svm);
    
    
    TFNL_svm = sum(tfnl_svm);
    TFNR_svm = sum(tfnr_svm);
    TFNF_svm = sum(tfnf_svm);
cv_acc_avg_svm = mean(cv_acc_svm);
cv_sensitivity_left_avg_svm = mean(cv_sensitivity_left_svm);  % averaged sensitivity
cv_specificity_left_avg_svm = mean(cv_specificity_left_svm);  % averaged specificity

cv_sensitivity_right_avg_svm = mean(cv_sensitivity_right_svm);  % averaged sensitivity
cv_specificity_right_avg_svm = mean(cv_specificity_right_svm);  % averaged specificity

cv_sensitivity_foot_avg_svm = mean(cv_sensitivity_foot_svm);  % averaged sensitivity
cv_specificity_foot_avg_svm = mean(cv_specificity_foot_svm);  % averaged specificity

cv_sensitivity_tongue_avg_svm = mean(cv_sensitivity_tongue_svm);  % averaged sensitivity
cv_specificity_tongue_avg_svm = mean(cv_specificity_tongue_svm);  % averaged specificity

SVM_CON_MATRIX= [TP_left_class_svm  LFNR_svm LFNF_svm LFNT_svm;  RFNL_svm TP_right_class_svm RFNF_svm RFNT_svm;  FFNL_svm FFNR_svm TP_foot_class_svm FFNT_svm; TFNL_svm TFNR_svm TFNF_svm TP_tongue_class_svm];
class_labels_svm = ["Left Class";"Right Class";"Foot Class";"Tongue Class"];

figure()
confusionplot_svm=confusionchart(SVM_CON_MATRIX,class_labels_svm,'Title','SVM Confusion Matrix 4 Classes');

J= 10;
tp_left_knn=zeros(size(J));
tn_left_knn=zeros(size(J));
fp_left_knn=zeros(size(J));
fn_left_knn=zeros(size(J));
tp_right_knn=zeros(size(J));
tn_right_knn=zeros(size(J));
fp_right_knn=zeros(size(J));
fn_right_knn=zeros(size(J));
tp_foot_knn=zeros(size(J));
tn_foot_knn=zeros(size(J));
fp_foot_knn=zeros(size(J));
fn_foot_knn=zeros(size(J));
tp_tongue_knn=zeros(size(J));
tn_tongue_knn=zeros(size(J));
fp_tongue_knn=zeros(size(J));
fn_tongue_knn=zeros(size(J));

lfnr_knn=zeros(size(J));
lfnf_knn=zeros(size(J));
lfnt_knn=zeros(size(J));

rfnl_knn=zeros(size(J));
rfnf_knn=zeros(size(J));
rfnt_knn=zeros(size(J));

ffnl_knn=zeros(size(J));
ffnr_knn=zeros(size(J));
ffnt_knn=zeros(size(J));

tfnl_knn=zeros(size(J));
tfnr_knn=zeros(size(J));
tfnf_knn=zeros(size(J));
indices_knn=crossvalind('Kfold',all_labels,M);
for j=1:J
    cv_test_idx_knn = find(indices_knn == j); % indices for test samples in one trial of validation
    cv_train_idx_knn = find(indices_knn ~= j); % indices for training samples in one trial of validation
    Mdl = fitcknn(all_samples(cv_train_idx_knn,:),all_labels(cv_train_idx_knn),'NumNeighbors',5,'NSMethod','exhaustive','Distance','minkowski');
    cv_classout_knn = predict(Mdl,all_samples(cv_test_idx_knn,:));
    cv_acc_knn(j) = mean(cv_classout_knn==all_labels(cv_test_idx_knn));
    TP_left_class_knn = sum((cv_classout_knn==all_labels(cv_test_idx_knn))&(cv_classout_knn==1));
    TN_left_class_knn = sum((cv_classout_knn==all_labels(cv_test_idx_knn))&(cv_classout_knn~=1));
    FN_left_class_knn = sum((cv_classout_knn~=1)&(all_labels(cv_test_idx_knn)==1));
    FP_left_class_knn = sum((cv_classout_knn==1)&(all_labels(cv_test_idx_knn)~=1));
    
   
    LFNR_knn = sum((cv_classout_knn~=all_labels(cv_test_idx_knn))&(cv_classout_knn==1)&(all_labels(cv_test_idx_knn)==2));
    LFNF_knn = sum((cv_classout_knn~=all_labels(cv_test_idx_knn))&(cv_classout_knn==1)&(all_labels(cv_test_idx_knn)==3));
    LFNT_knn = sum((cv_classout_knn~=all_labels(cv_test_idx_knn))&(cv_classout_knn==1)&(all_labels(cv_test_idx_knn)==4));
    

    TP_right_class_knn = sum((cv_classout_knn==all_labels(cv_test_idx_knn))&(cv_classout_knn==2));
    TN_right_class_knn = sum((cv_classout_knn==all_labels(cv_test_idx_knn))&(cv_classout_knn~=2));
    FN_right_class_knn = sum((cv_classout_knn~=2)&(all_labels(cv_test_idx_knn)==2));
    FP_right_class_knn = sum((cv_classout_knn==2)&(all_labels(cv_test_idx_knn)~=2));
    
    RFNL_knn = sum((cv_classout_knn~=all_labels(cv_test_idx_knn))&(cv_classout_knn==2)&(all_labels(cv_test_idx_knn)==1));
    RFNF_knn = sum((cv_classout_knn~=all_labels(cv_test_idx_knn))&(cv_classout_knn==2)&(all_labels(cv_test_idx_knn)==3));
    RFNT_knn = sum((cv_classout_knn~=all_labels(cv_test_idx_knn))&(cv_classout_knn==2)&(all_labels(cv_test_idx_knn)==4));
    

    TP_foot_class_knn = sum((cv_classout_knn==all_labels(cv_test_idx_knn))&(cv_classout_knn==3));
    TN_foot_class_knn = sum((cv_classout_knn==all_labels(cv_test_idx_knn))&(cv_classout_knn~=3));
    FN_foot_class_knn = sum((cv_classout_knn~=3)&(all_labels(cv_test_idx_knn)==3));
    FP_foot_class_knn = sum((cv_classout_knn==3)&(all_labels(cv_test_idx_knn)~=3));
    
    FFNL_knn = sum((cv_classout_knn~=all_labels(cv_test_idx_knn))&(cv_classout_knn==3)&(all_labels(cv_test_idx_knn)==1));
    FFNR_knn = sum((cv_classout_knn~=all_labels(cv_test_idx_knn))&(cv_classout_knn==3)&(all_labels(cv_test_idx_knn)==2));
    FFNT_knn = sum((cv_classout_knn~=all_labels(cv_test_idx_knn))&(cv_classout_knn==3)&(all_labels(cv_test_idx_knn)==4));


    TP_tongue_class_knn = sum((cv_classout_knn==all_labels(cv_test_idx_knn))&(cv_classout_knn==4));
    TN_tongue_class_knn = sum((cv_classout_knn==all_labels(cv_test_idx_knn))&(cv_classout_knn~=4));
    FN_tongue_class_knn = sum((cv_classout_knn~=4)&(all_labels(cv_test_idx_knn)==4));
    FP_tongue_class_knn = sum((cv_classout_knn==4)&(all_labels(cv_test_idx_knn)~=4));
    
    TFNL_knn = sum((cv_classout_knn~=all_labels(cv_test_idx_knn))&(cv_classout_knn==4)&(all_labels(cv_test_idx_knn)==1));
    TFNR_knn = sum((cv_classout_knn~=all_labels(cv_test_idx_knn))&(cv_classout_knn==4)&(all_labels(cv_test_idx_knn)==2));
    TFNF_knn = sum((cv_classout_knn~=all_labels(cv_test_idx_knn))&(cv_classout_knn==4)&(all_labels(cv_test_idx_knn)==3));

    cv_sensitivity_left_knn(j) = TP_left_class_knn/(TP_left_class_knn+FN_left_class_knn);%Tells what fraction of predictions as left class were actually left
    cv_specificity_left_knn(j) = TN_left_class_knn/(TN_left_class_knn+FP_left_class_knn);%Tells what fraction of all not left samples are correctly predicted as not left by the classifier.

    cv_sensitivity_right_knn(j) = TP_right_class_knn/(TP_right_class_knn+FN_right_class_knn);%Tells what fraction of predictions as right class were actually right
    cv_specificity_right_knn(j) = TN_right_class_knn/(TN_right_class_knn+FP_right_class_knn);%Tells what fraction of all not right samples are correctly predicted as not right by the classifier.

    cv_sensitivity_foot_knn(j) = TP_foot_class_knn/(TP_foot_class_knn+FN_foot_class_knn);%Tells what fraction of predictions as foot class were actually foot
    cv_specificity_foot_knn(j) = TN_foot_class_knn/(TN_foot_class_knn+FP_foot_class_knn);%Tells what fraction of all not foot samples are correctly predicted as not foot by the classifier.

    cv_sensitivity_tongue_knn(j) = TP_tongue_class_knn/(TP_tongue_class_knn+FN_tongue_class_knn);%Tells what fraction of predictions as tongue class were actually tongue 
    cv_specificity_tongue_knn(j) = TN_tongue_class_knn/(TN_tongue_class_knn+FP_tongue_class_knn);%Tells what fraction of all not tongue samples are correctly predicted as not tongue by the classifier.
    tp_left_knn(j)=TP_left_class_knn;
    tn_left_knn(j)=TN_left_class_knn;
    fp_left_knn(j)=FP_left_class_knn;
    fn_left_knn(j)=FN_left_class_knn;
    
    lfnr_knn(j)=LFNR_knn;
    lfnf_knn(j)=LFNF_knn;
    lfnt_knn(j)=LFNT_knn;
    
    
    
    tp_right_knn(j)=TP_right_class_knn;
    tn_right_knn(j)=TN_right_class_knn;
    fp_right_knn(j)=FP_right_class_knn;
    fn_right_knn(j)=FN_right_class_knn;
    
    rfnl_knn(j)=RFNL_knn;
    rfnf_knn(j)=RFNF_knn;
    rfnt_knn(j)=RFNT_knn;
    
    tp_foot_knn(j)=TP_foot_class_knn;
    tn_foot_knn(j)=TN_foot_class_knn;
    fp_foot_knn(j)=FP_foot_class_knn;
    fn_foot_knn(j)=FN_foot_class_knn;
    
    ffnl_knn(j)=FFNL_knn;
    ffnr_knn(j)=FFNR_knn;
    ffnt_knn(j)=FFNT_knn;
    
    tp_tongue_knn(j)=TP_tongue_class_knn;
    tn_tongue_knn(j)=TN_tongue_class_knn;
    fp_tongue_knn(j)=FP_tongue_class_knn;
    fn_tongue_knn(j)=FN_tongue_class_knn;
    
    tfnl_knn(j)=TFNL_knn;
    tfnr_knn(j)=TFNR_knn;
    tfnf_knn(j)=TFNF_knn;

end
    TP_left_class_knn = sum(tp_left_knn);
    TN_left_class_knn = sum(tn_left_knn);
    FP_left_class_knn = sum(fp_left_knn);
    FN_left_class_knn = sum(fn_left_knn);
    
    LFNR_knn = sum(lfnr_knn);
    LFNF_knn = sum(lfnf_knn);
    LFNT_knn = sum(lfnt_knn);

    
    TP_right_class_knn = sum(tp_right_knn);
    TN_right_class_knn = sum(tn_right_knn);
    FP_right_class_knn = sum(fp_right_knn);
    FN_right_class_knn = sum(fn_right_knn);
    
    RFNL_knn = sum(rfnl_knn);
    RFNF_knn = sum(rfnf_knn);
    RFNT_knn = sum(rfnt_knn);
   
    
    TP_foot_class_knn = sum(tp_foot_knn);
    TN_foot_class_knn = sum(tn_foot_knn);
    FP_foot_class_knn = sum(fp_foot_knn);
    FN_foot_class_knn = sum(fn_foot_knn);
    
    FFNL_knn = sum(ffnl_knn);
    FFNR_knn = sum(ffnr_knn);
    FFNT_knn = sum(ffnt_knn);
    
    
    TP_tongue_class_knn = sum(tp_tongue_knn);
    TN_tongue_class_knn = sum(tn_tongue_knn);
    FP_tongue_class_knn = sum(fp_tongue_knn);
    FN_tongue_class_knn = sum(fn_tongue_knn);
    
    
    TFNL_knn = sum(tfnl_knn);
    TFNR_knn = sum(tfnr_knn);
    TFNF_knn = sum(tfnf_knn);
cv_acc_avg_knn = mean(cv_acc_knn);
cv_sensitivity_left_avg_knn = mean(cv_sensitivity_left_knn);  % averaged sensitivity
cv_specificity_left_avg_knn = mean(cv_specificity_left_knn);  % averaged specificity

cv_sensitivity_right_avg_knn = mean(cv_sensitivity_right_knn);  % averaged sensitivity
cv_specificity_right_avg_knn = mean(cv_specificity_right_knn);  % averaged specificity

cv_sensitivity_foot_avg_knn = mean(cv_sensitivity_foot_knn);  % averaged sensitivity
cv_specificity_foot_avg_knn = mean(cv_specificity_foot_knn);  % averaged specificity

cv_sensitivity_tongue_avg_knn = mean(cv_sensitivity_tongue_knn);  % averaged sensitivity
cv_specificity_tongue_avg_knn = mean(cv_specificity_tongue_knn);  % averaged specificity

KNN_CON_MATRIX= [TP_left_class_knn  LFNR_knn LFNF_knn LFNT_knn;  RFNL_knn TP_right_class_knn RFNF_knn RFNT_knn;  FFNL_knn FFNR_knn TP_foot_class_knn FFNT_knn; TFNL_knn TFNR_knn TFNF_knn TP_tongue_class_knn];
class_labels_knn = ["Left Class";"Right Class";"Foot Class";"Tongue Class"];

figure()
confusionplot_knn=confusionchart(KNN_CON_MATRIX,class_labels_knn,'Title','KNN Confusion Matrix 4 Classes');
