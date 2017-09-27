function [tauX1,tauX2,tauY1,tauY2,tauPer1,tauPer2] = ...
    FindTaus(responsedata,threshold,responseFinalValue,responseInitialValue)
    %% Calculate tau using initial slope method
    
    % initialize variables
    indSt = 1;
    
    % cut off the initial flat data
    for i = 1:size(responsedata,1)
        % This condition will be met after the initial index, therefore if this
        % condition is met, the loop will exit
        if (responsedata(i,1) > threshold)
            indSt = i;
            break;
        end
    end
    indEnd = indSt + 4;
    
    % polyfit the initial slope segment data to find the equation of a line
    % then use the equation values and final value to find the intersection
    p = polyfit(responsedata(indSt:indEnd,1),responsedata(indSt:indEnd,2),1);
    tauX1 = (responseFinalValue - p(2))/p(1); % x = (y - y(0))/m
    for i = indEnd:size(responsedata,1) % we know the tau must be beyond indEnd
        if (responsedata(i,1) >= tauX1)
            tauY1 = responsedata(i,2);
            break;
        end
    end

    % check for the tau percentage
    tauPer1 = (tauY1-responseInitialValue)/(responseFinalValue-responseInitialValue);

    %% Calculate tau using the 63.2% method
    tauY2 = (responseFinalValue-responseInitialValue)*0.632 + responseInitialValue;
    for i = indEnd:size(responsedata,1)
        if (responsedata(i,2) >= tauY2)
            tauX2 = responsedata(i,1);
            break;
        end
    end

    % check for the tau percentage
    tauPer2 = (tauY2-responseInitialValue)/(responseFinalValue-responseInitialValue);
    
end % end of function