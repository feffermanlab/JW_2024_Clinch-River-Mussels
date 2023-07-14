clear all

load('mature_Dem5.mat')
load('State_Sub5_Prop.mat')
%Set the number of beds
nBeds = 5;
% Set number of age classes
nAge = 3;
%Set initial population states
Istate = [zeros(2,nBeds); 0 0 0 0 .1*State_Sub5(3,5)];
Sstate = State_Sub5 - Istate;
% Sstate = [zeros(1,nBeds);
%     12.7072   30.8097   55.6399   75.9820   69.9554;
%    42.6071  103.3045  186.5601  254.7670  234.5599];

% Set projection matrices
reduction = .9;
pS_surv = 0.971; %survival prob for sus pop
pI_surv = pS_surv-reduction*.971; %survival prob for inf pop
fecS = 10^6;
fecI = (1-reduction)*fecS;
MS_r = [pS_surv*(11/12) 0;  pS_surv*1/12 pS_surv*1]; %aging during reproduction for sus pop
MI_r = [pI_surv*(11/12) 0;  pI_surv*1/12 pI_surv*1]; %aging during reproduction for inf pop
MS_a = [0 0 0; .21*(1/fecS) pS_surv*(11/12) 0; 0 pS_surv*1/12 pS_surv*1]; %aging matrix for sus pop
MI_a = [0 0 0; .21*(1-reduction)*(1/fecS) pI_surv*(11/12) 0; 0 pI_surv*1/12 pI_surv*1]; %aging matrix for inf pop

%Set infection probability
p_inf = 0.01;

%Bed distribution proportions
b_0 = 0.3;
b_m1 = 0.2;
b_p1 = 0.3;
b_p2 = 0.2;

steps = 20;
mature = Sstate(3,:)+Istate(3,:);
mature = [mature ;zeros(steps,nBeds)];

for step_counter = 1:steps
    Istate_copy = Istate;
    Sstate_copy = Sstate;
    if mod(step_counter,4) ==0 
        %Reproduction season

        for bed_counter = 1:nBeds
            numI = sum(Istate_copy(:,bed_counter));
            theta = numI ~= 0;
            Sstate(2:3,bed_counter) = (1-theta*p_inf)*MS_r*Sstate_copy(2:3,bed_counter); %age
            Istate(2:3,bed_counter) = MI_r*Istate_copy(2:3,bed_counter)+theta*p_inf*MS_r*Sstate_copy(2:3,bed_counter);

            Sstate(1,bed_counter) = Sstate(1,bed_counter)+b_0*(1-theta*p_inf)*fecS*Sstate_copy(3,bed_counter); %same bed
            Sstate(1,max(1,bed_counter-1)) = Sstate(1,max(1,bed_counter-1)) + (bed_counter-1>0)*(b_m1*(1-theta*p_inf)*fecS*Sstate_copy(3,bed_counter)); %goch sent upstream
            Sstate(1,min(nBeds,bed_counter+1)) = Sstate(1,min(nBeds,bed_counter+1)) + (bed_counter+1<=nBeds)*(b_p1*(1-theta*p_inf)*fecS*Sstate_copy(3,bed_counter)); %goch 1 downstream
            Sstate(1,min(nBeds,bed_counter+2)) = Sstate(1,min(nBeds,bed_counter+2)) + (bed_counter+2<=nBeds)*(b_p2*(1-theta*p_inf)*fecS*Sstate_copy(3,bed_counter)); %goch 2 downstream

            Istate(1,bed_counter) = Istate(1,bed_counter)+b_0*theta*p_inf*fecS*Sstate_copy(3,bed_counter)+b_0*fecI*Istate_copy(3,bed_counter); %same bed
            Istate(1,max(1,bed_counter-1)) = Istate(1,max(1,bed_counter-1)) + (bed_counter-1>0)*((b_m1*theta*p_inf*fecS*Sstate_copy(3,bed_counter))+b_m1*fecI*Istate_copy(3,bed_counter)); %goch sent upstream
            Istate(1,min(nBeds,bed_counter+1)) = Istate(1,min(nBeds,bed_counter+1)) + (bed_counter+1<=nBeds)*((b_p1*theta*p_inf*fecS*Sstate_copy(3,bed_counter))+b_p1*fecI*Istate_copy(3,bed_counter)); %goch 1 downstream
            Istate(1,min(nBeds,bed_counter+2)) = Istate(1,min(nBeds,bed_counter+2)) + (bed_counter+2<=nBeds)*((b_p2*theta*p_inf*fecS*Sstate_copy(3,bed_counter))+b_p2*fecI*Istate_copy(3,bed_counter)); %goch 2 downstream

        end
    else
        %No reproduction
        for bed_counter = 1:nBeds
        numI = sum(Istate(:,bed_counter));
        theta = numI ~= 0;
        Sstate(:,bed_counter) = (1-theta*p_inf)*MS_a*Sstate_copy(:,bed_counter);
        Istate(:,bed_counter) = theta*p_inf*MS_a*Sstate_copy(:,bed_counter)+MI_a*Istate_copy(:,bed_counter);
        end
    end
    mature(step_counter+1,:) = Sstate(3,:)+Istate(3,:);
end
% env_cont_level
% bed_cont_level
% state
MarkerIndex = 1:4:steps;
MarkerShapes = {'^','o','s','d','v'};
h=figure('Position',.6*get(0,'ScreenSize'));
hold on
for k = 1:nBeds
    hold on
    plot(0:steps,mature(:,k),'LineWidth',2,'Marker',MarkerShapes(k),'MarkerIndices',MarkerIndex,'MarkerSize',12);
end
ylim padded;
legend('1','2','3','4','5')
title("Pathogen")
set(gca,'FontSize',32)
xlabel('Time');
ylabel('Mature Pop Size')

h2=figure('Position',.6*get(0,'ScreenSize'));
hold on
for k = 1:nBeds
    hold on
    plot(0:steps,abs(mature_Dem5(:,k)-mature(:,k))./mature_Dem5(:,k),'LineWidth',2,'Marker',MarkerShapes(k),'MarkerSize',12);
end
ylim padded;
legend('1','2','3','4','5')
title("Pathogen")
set(gca,'FontSize',32)
xlabel('Time');
ylabel('Percent Change')

