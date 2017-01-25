function [ intersect, closest_linept, second_closest_linept ] = pt2line_intersect( line, point )
%This function takes an XY point that is not on a line and finds where it
%would intersect the line by way of the shortest distance between line and
%point. Both point and line need to be X=column1, Y=column2
% by Rachel Kay 9/10/12
[linesize,~] = size(line);
distances= zeros(linesize,1);
for j = 1: linesize
        distances(j,1) = sqrt((point(1,1)-line(j,1))^2+((point(1,2)-line(j,2))^2));
end
    [~, closestindex] = sort(distances(:,1));
    [X1]= line(closestindex(1),1);% closest two line points
    [X2]= line(closestindex(2),1);
    [Y1]= line(closestindex(1),2);
    [Y2]= line(closestindex(2),2);
    [XD]= point(1,1); %  point
    [YD]= point(1,2);
    u = ((XD-X1)*(X2-X1)+(YD-Y1)*(Y2-Y1))/((norm([X1;Y1]-[X2;Y2]))^2);
    xintersct = X1 + u*(X2-X1);
    yintersct = Y1 + u*(Y2-Y1);
    intersect = [xintersct, yintersct];
    closest_linept = [line(closestindex(1),:),closestindex(1)];
    second_closest_linept = [line(closestindex(2),:),closestindex(2)];
    

end

