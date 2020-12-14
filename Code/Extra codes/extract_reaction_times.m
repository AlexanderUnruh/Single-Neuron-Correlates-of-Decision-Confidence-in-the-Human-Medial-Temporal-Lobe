function [RT,endResponse] = extract_reaction_times(events)

% This function extracts the real reaction times of satiation,
% i.e., the time between the likert is shown to the first
% press buttom
%
% Input: events of satiation. Events contains the timestamps of each event
% in the first column and the event_ID in the second
% Output: reaction times

ind1=find(events(:,2)==50);


if all(events(ind1-1,2)==44) &&  length(ind1)==120
    RT{2}=events(ind1(1:60)+1,1) - events(ind1(1:60),1);
    RT{5}=events(ind1(61:end)+1,1) - events(ind1(61:end),1);
    
else
    error ('Event 44 should be followed by event 50 in rating trials and 120 trials, some error here!')
end


ind3=find(events(:,2)==51);

if all(events(ind3-1,2) == 43) && length(ind3)== 380
    RT{3}=events(ind3(1:190)+1,1) - events(ind3(1:190),1);  % difference between time of likert shown (51) and the next ttl
    RT{6}=events(ind3(191:end)+1,1) - events(ind3(191:end),1);
    
else
    error ('Event 43 should be followed by event 51 in 2AFC trials and 380 trials, some error here!')
end




ind5 = find(events(:,2)==50);


if all(events(ind5-1,2)==44) &&  length(ind5)==120
    for kk=1:120
        if kk<60
            mm=find(events(ind5(kk):end,2)==2,1,'first')-2; events(ind5(kk)+mm:ind5(kk)+mm+1,2);
            if events(ind5(kk)+mm,2) > 15
                endResponse{2}(kk,1)= (events(ind5(kk)+mm,1) - events(ind5(kk),1)) - (0.3 *1e3);
            else
                error('This event should be bigger tha 15 to be a likert event')
            end
        elseif kk==60
            mm=find(events(ind5(kk):end,2)==1,1,'first')-2; events(ind5(kk)+mm:ind5(kk)+mm+1,2);
            
            if events(ind5(kk)+mm,2) > 15
                endResponse{2}(kk,1)= (events(ind5(kk)+mm,1) - events(ind5(kk),1)) - (0.3 *1e3);
            else
                error('This event should be bigger tha 15 to be a likert event')
            end
        elseif kk>60 && kk<120
            mm=find(events(ind5(kk):end,2)==2,1,'first')-2; events(ind5(kk)+mm:ind5(kk)+mm+1,2);
            if events(ind5(kk)+mm,2) > 15
                endResponse{5}(kk-60,1)= (events(ind5(kk)+mm,1) - events(ind5(kk),1)) - (0.3 *1e3);
            else
                error('This event should be bigger tha 15 to be a likert event')
            end
        elseif kk==120
            mm=find(events(ind5(kk):end,2)==1,1,'first')-2; events(ind5(kk)+mm:ind5(kk)+mm+1,2);
            if events(ind5(kk)+mm,2) > 15
                endResponse{5}(kk-60,1)= (events(ind5(kk)+mm,1) - events(ind5(kk),1)) - (0.3 *1e3);
            else
                error('This event should be bigger tha 15 to be a likert event')
            end
        end
    end
    
else
    error ('Event 44 should be followed by event 50 in rating trials and 120 trials, some error here!')
end






%% (0.3 *1e3) are deleted because aftr the ok press, a WaitSec(0.3) is run. Thus, to account to the exact time of the ok press.


ind6 = find(events(:,2)==51)-1;

if all(events(ind6,2) == 43) && length(ind6)== 380
    
    for kk=1:380
        if kk <=190
            mm = find(events(ind6(kk):end,2)==32,1,'first')-2; events(ind6(kk)+mm:ind6(kk)+mm+1,2);
            if events(ind6(kk)+mm,2) > 15
                endResponse{3}(kk,1) = (events(ind6(kk)+mm,1) - events(ind6(kk),1)) - (0.3 *1e3);
            else
                error('This event should be bigger tha 15 to be a likert event')
            end
            
        elseif kk>190
            mm = find(events(ind6(kk):end,2)==32,1,'first')-2; events(ind6(kk)+mm:ind6(kk)+mm+1,2);
            if events(ind6(kk)+mm,2) > 15
                endResponse{6}(kk-190,1)= (events(ind6(kk)+mm,1) - events(ind6(kk),1)) - (0.3 *1e3);
            else
                error('This event should be bigger tha 15 to be a likert event')
            end
            
        end
    end
    
 
    
else
    error ('Event 43 should have 380 trials, some error here!')
end