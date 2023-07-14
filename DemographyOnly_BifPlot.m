steps = 1000;
% Setting up for plot
h2=figure('Position',.6*get(0,'ScreenSize'));
hold on
nDots = 4;
EndIndex = steps-(nDots-1)*4:4:steps;

Colors = [      0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880];
%%%%%%%%%%%%%%%%%%%%

for delta = -50:2:50
%Set the number of beds
nBeds = 5;
% Set number of age classes
nAge = 3;
%Set initial population states
%state = [zeros(1,nBeds);100*ones(nAge-2,nBeds);100*ones(1,nBeds)];
% state = [27750291.4022342	67282963.1439813	121507881.915295	165931553.242432	152770508.374037;
% 10.5236939322237	25.5155991131068	46.0792191357980	62.9259294156950	57.9348896400269;
% 35.2858416778074	85.5535513330852	154.503165520762	210.990017987569	194.255110903730]; %steady state
% state = [zeros(1,nBeds);
%     100   30.8097   55.6399   75.9820   69.9554;
%    100  103.3045  186.5601  254.7670  234.5599]; %recovers to steady
%    state
state = [100 100 100 100 100;
50 50 50 50 50;
100 100 100 100 100]+delta*ones(3,5); 

% Set projection matrices
p_surv = 0.971;
fec = 10^6;
M_r = [p_surv*(11/12) 0;  p_surv*1/12 p_surv*1]; %aging during reproduction
M_a = [0 0 0; .20024*(1/fec) p_surv*(11/12) 0; 0 p_surv*1/12 p_surv*1]; %aging matrix


%Bed distribution proportions
b_0 = 0.3;
b_m1 = 0.2;
b_p1 = 0.3;
b_p2 = 0.2;

bed_pop = sum(state(2:3,:));

mature = state(3,:);
mature = [mature; zeros(steps,nBeds)];

for step_counter = 1:steps
    %update bed population
    if mod(step_counter,4) ==0
        %Reproduction season
        state_copy = state;
        for bed_counter = 1:nBeds
            state(2:3,bed_counter) = M_r*state_copy(2:3,bed_counter); %age
            state(1,bed_counter) = state(1,bed_counter)+b_0*fec*state_copy(3,bed_counter); %same bed
            state(1,max(1,bed_counter-1)) = state(1,max(1,bed_counter-1)) + (bed_counter-1>0)*(b_m1*fec*state_copy(3,bed_counter)); %goch sent upstream
            state(1,min(nBeds,bed_counter+1)) = state(1,min(nBeds,bed_counter+1)) + (bed_counter+1<=nBeds)*(b_p1*fec*state_copy(3,bed_counter)); %goch 1 downstream
            state(1,min(nBeds,bed_counter+2)) = state(1,min(nBeds,bed_counter+2)) + (bed_counter+2<=nBeds)*(b_p2*fec*state_copy(3,bed_counter)); %goch 2 downstream

        end
    else
        %No reproduction
        state_copy = state;
        for bed_counter = 1:nBeds
        state(:,bed_counter) = M_a*state_copy(:,bed_counter);
        end
    end
    mature(step_counter+1,:) = state(3,:);
end

time = 0:steps;
MarkerIndex = 1:50:steps;
MarkerShapes = {'^','o','s','d','v'};
% 
% open("BaseFigure.fig");

% h=figure('Position',.6*get(0,'ScreenSize'));
% hold on
% for k = 1:nBeds
%     hold on
%     plot(time,mature(:,k),'LineWidth',1,'Color',Colors(k,:));
% end
% ylim padded;
% legend('1','2','3','4','5')
% title("Base Model")
% set(gca,'FontSize',24)
% xlabel('Time');
% ylabel('Mature Population Size')




for k = 1:nBeds
    hold on
    scatter(delta*ones(nDots,1),mature(EndIndex,k),'LineWidth',2,'Marker',MarkerShapes(k),'MarkerEdgeColor',Colors(k,:));
end
ylim padded;
set(gca,'FontSize',24)
xlabel('$\delta$-Value','Interpreter','latex');
ylabel('Mature Population Size')
end

