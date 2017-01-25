function [] = single_section_analysis(clockwise, output_file, section )
% example function call:
% single_section_analysis('example.xlsx',1,'example_matlabed.xlsx')
addpath('C:\Users\Rachel\Desktop\AONanalysis') %to change input folder change this path

%for an individual section from neurolucida, data must be imported into
%matlab, and must include AON boundaries and a lateralis-dorsalis boundary,
%cells and background. This script puts all the data into usable info for
%matlab
% by Rachel Kay 9/14/2012

section.textdata.Sheet2 = section.textdata.Sheet2(2:end,:);
section.textdata.Sheet3 = section.textdata.Sheet3(2:end,:);


%%
%1st this script finds the boundaries in the data and makes a binary array
%to show if they are present in the following order
bndstrg = {'AON outer boundary';'AON inner boundary';'externa boundary';'lateralis-dorsalis boundary';'dorsalis-medialis boundary';'lateralis-vp boundary';'med-vp boundary'};
boundaries = zeros(7,1);
for i = 1:7
    boundaries(i,1)= sum(sum(strncmp(bndstrg(i),section.textdata.Sheet3,11)));
end

%%

%put the X Y coordinates with the boundaries
if boundaries(1,1) == 0
    error('No AON outer boundary')
else outer= section.data.Sheet1((min(strmatch(bndstrg(1),section.textdata.Sheet2))+7):(max(strmatch(bndstrg(1),section.textdata.Sheet2))),26:27);
end
if boundaries(2,1) == 0
    error('No AON inner boundary')
    else inner= section.data.Sheet1((min(strmatch(bndstrg(2),section.textdata.Sheet2))+7):(max(strmatch(bndstrg(2),section.textdata.Sheet2))),26:27);
end
if boundaries (3,1) == 1
    externa= section.data.Sheet1((min(strmatch(bndstrg(3),section.textdata.Sheet2))+7):(max(strmatch(bndstrg(3),section.textdata.Sheet2))),26:27);
else
    externa = zeros(2,2);
end
if boundaries (4,1) == 1
    lat_dors= section.data.Sheet1((min(strmatch(bndstrg(4),section.textdata.Sheet2))+7):(max(strmatch(bndstrg(4),section.textdata.Sheet2))),26:27);
else
    lat_dors = zeros(2,2);
end
if boundaries (5,1) == 1
    dors_med= section.data.Sheet1((min(strmatch(bndstrg(5),section.textdata.Sheet2))+7):(max(strmatch(bndstrg(5),section.textdata.Sheet2))),26:27);
else
    dors_med = zeros(2,2);
end
if boundaries (6,1) == 1
    lat_vp= section.data.Sheet1((min(strmatch(bndstrg(6),section.textdata.Sheet2))+7):(max(strmatch(bndstrg(6),section.textdata.Sheet2))),26:27);
else
    lat_vp= zeros(2,2);
end
if boundaries (7,1) == 1
    med_vp= section.data.Sheet1((min(strmatch(bndstrg(7),section.textdata.Sheet2))+7):(max(strmatch(bndstrg(7),section.textdata.Sheet2))),26:27);
else
    med_vp = zeros(2,2);
end
%%

%get the x y coordinates and brightness for the CELLS

xmlcells = strmatch('CELL',section.textdata.Sheet2); %each cells takes up four rows in the xml, we need to narrow it down to just the center XY pts
cellindexml = xmlcells(3:4:end);
clipcells = strmatch('CELL',section.textdata.Sheet3);

cells1(:,1)= section.data.Sheet1((cellindexml),26); % x coords
cells1(:,2)= section.data.Sheet1((cellindexml),27); % y coords
cells1(:,3)= section.data.Sheet3((clipcells),5); % brightness
if isempty(cells1(:,3))
    error('No brightness measurements for cells, check countor measurements table')
end

%%

%get the background mean and sd
section.textdata.Sheet3 = section.textdata.Sheet3(2:end, :);
bkgrdindex = strmatch('BACKGROUND', section.textdata.Sheet3);
bkgrd = section.data.Sheet3((bkgrdindex),5);
if isempty(bkgrd)
    error('No BACKGROUND measurements')
end
threshold = mean(bkgrd)-(1.5*std(bkgrd)); %%%%%% RIGHT HERE IS WHERE YOU CAN CHANGE 
%THRESHOLD. IT IS THE MEAN MINUS 1.5X STANDARD DEVIATION or something else
%you might want, maybe?
%use the threshold to cut off cells
threshed_cells_index = cells1(:,3)< threshold;
cells = cells1(threshed_cells_index,1:2);

if clockwise ~=1
    inner = flipud(inner);
    outer = flipud(outer);
    externa = flipud(externa);
end
%%

%Finding the Midline- thie midline is found by taking every point on the
%outer boundary finding the closest point on the inner boundary to it, then
%placing a point halfway between them
[innersize,~] = size(inner);
[outersize,~] = size(outer);
distances = zeros (3,(outersize*innersize));
%we will be calculating all the distances between each outer and every
%inner to find the lowest
closestpairXs = zeros(outersize,2);
closestpairYs = zeros(outersize,2);
midpoints= zeros(outersize,2);
for i = 1:outersize
    for j = 1:innersize
        distances(((i-1)*innersize+j),3) = sqrt((inner(j,1)-outer(i,1))^2+((inner(j,2)-outer(i,2))^2));
        distances(((i-1)*innersize+j),1)= i;
        distances(((i-1)*innersize+j),2)= j;
    end
    [~, innermatch] = sort(distances(((i-1)*innersize+1):((i-1)*innersize+innersize),3));
    closestpairXs(i,:)= [outer(i,1),inner(innermatch(1),1)];
    closestpairYs(i,:)= [outer(i,2),inner(innermatch(1),2)];
    midpoints(i,1) = (sum(closestpairXs(i,:)))/2;
    midpoints(i,2) = (sum(closestpairYs(i,:)))/2;
end
%%
%now we have to put the midpoints on a line instead of XY coords, the line
%is based on taking the distance between the first point on the midline and
%the second, and making pt1 = 0, pt2= distance bt 1 and
%2, and pt3 will be added based on its distance to pt2, and so on
point_distances = zeros(outersize,1);
midline_unrolled =zeros(outersize,2);
for i = 1:(outersize-1)
   point_distances((i+1),1)= sqrt((midpoints((i+1),1)-midpoints(i,1))^2+((midpoints((i+1),2)-midpoints(i,2))^2));
end
for i = 2:outersize
    midline_unrolled(i,1) = point_distances(i,1)+midline_unrolled((i-1),1);
end
%%
%if externa present in section, separate out externa cells, otherwise put
%cells on line
if externa ==0
    externa_unrolled = 0;
    pEcells_on_line = zeros(1,3);
    pEcells = zeros(1,2);
    pPcells_on_line = zeros(length(cells),1);
    pPcells = cells;
    for i = 1: length(cells)
        [out_intersect, ~, ~] = pt2line_intersect(outer, cells(i,:));
        [mid_intersect, p1mid, p2mid] = pt2line_intersect(midpoints, cells(i,:));
        distmidint2outer = sqrt((out_intersect(1,1)-mid_intersect(1,1))^2+(out_intersect(1,2)-mid_intersect(1,2))^2);
        distcell2outer = sqrt((cells(i,1)-out_intersect(1,1))^2+(cells(i,2)-out_intersect(1,2))^2);
        distcell2mid = sqrt((cells(i,1)-mid_intersect(1,1))^2+(cells(i,2)-mid_intersect(1,2))^2);
        distP1mid2int = sqrt((p1mid(1,1)-mid_intersect(1,1))^2+(p1mid(1,2)-mid_intersect(1,2))^2);
        distP1mid2outer = sqrt((p1mid(1,1)-(outer(p1mid(1,3),1)))^2+(p1mid(1,2)-(outer(p1mid(1,3),2)))^2);
        if p1mid(1,3)< p2mid(1,3)
            pPcells_on_line(i,1) =midline_unrolled(p1mid(1,3),1)+distP1mid2int;
        else
            pPcells_on_line(i,1) =midline_unrolled(p1mid(1,3),1)-distP1mid2int;
        end
        if distmidint2outer > distcell2outer
            pPcells_on_line(i,2) = distcell2mid;
        else
            pPcells_on_line(i,2) = distcell2mid *-1;
        end
        pPcells_on_line(i,3) = pPcells_on_line(i,2)/distP1mid2outer;
        if pPcells_on_line(i,3)>1
            pPcells_on_line(i,3) =1;
        end
        if pPcells_on_line(i,3)< -1
            pPcells_on_line(i,3)= -1;
        end
    end
else %externa is present
    %find externa midline
    [externasize,~] = size(externa);
    extdistances = zeros (3,(outersize*externasize));
    extclosestpairXs = zeros(externasize,2);
    extclosestpairYs = zeros(externasize,2);
    extmidpoints= zeros(externasize,2);
    for i = 1:externasize  
        for j = 1:outersize
        extdistances(((i-1)*outersize+j),3) = sqrt((outer(j,1)-externa(i,1))^2+((outer(j,2)-externa(i,2))^2));
        extdistances(((i-1)*outersize+j),1)= i;
        extdistances(((i-1)*outersize+j),2)= j;
        end
    [~, outermatch] = sort(extdistances(((i-1)*outersize+1):((i-1)*outersize+outersize),3));
    extclosestpairXs(i,:)= [externa(i,1),outer(outermatch(1),1)];
    extclosestpairYs(i,:)= [externa(i,2),outer(outermatch(1),2)];
    extmidpoints(i,1) = (sum(extclosestpairXs(i,:)))/2;
    extmidpoints(i,2) = (sum(extclosestpairYs(i,:)))/2;
    externa_midline = extmidpoints;
    end
    extpt_distances = zeros(externasize,1);
    externa_unrolled =zeros(externasize,2);
    for j = 1:(externasize-1)
        extpt_distances((j+1),1)= sqrt((externa_midline((j+1),1)-externa_midline(j,1))^2+((externa_midline((j+1),2)-externa_midline(j,2))^2));
    end
    for j = 2:externasize
    externa_unrolled(j,1) = extpt_distances(j,1)+externa_unrolled((j-1),1);
    end
    %take externa cells off the list of cells, and put them on externa line
    if isempty(cells)
        pEcells = 0;
        pEcells_on_line = 0;
        pPcells = 0;
        pPcells_on_line = 0;
    else
        pPcells = zeros(length(cells),2);
        pPcells_on_line = zeros(length(cells),3);
        pEcells = zeros(length(cells),2);
        pEcells_on_line = zeros(length(cells),3);%column1 is place on flat line colm 2 is real distance + or - from flat line, colm 3 is scaled to -1 to +1
        for i = 1:length(cells)
            [emid_intersect, p1emid, p2emid] = pt2line_intersect(externa_midline, cells(i,:));
            [out_intersect, p1out, ~] = pt2line_intersect(outer, cells(i,:));
            [ext_intersect, p1ext, ~] = pt2line_intersect(externa, cells(i,:));
            [mid_intersect, p1mid, p2mid] = pt2line_intersect(midpoints, cells(i,:));
            
            distcell2mid = sqrt((cells(i,1)-mid_intersect(1,1))^2+(cells(i,2)-mid_intersect(1,2))^2);
            distP1out2ext = sqrt((p1out(1,1)-p1ext(1,1))^2+(p1out(1,2)-p1ext(1,2))^2);
            distP1emid2int = sqrt((p1emid(1,1)-emid_intersect(1,1))^2+(p1emid(1,2)-emid_intersect(1,2))^2);
            distcell2ext = sqrt((cells(i,1)-ext_intersect(1,1))^2+(cells(i,2)-ext_intersect(1,2))^2);
            distemid2ext = sqrt((p1emid(1,1)-p1ext(1,1))^2+(p1emid(1,2)-p1ext(1,2))^2);
            distemid2cell = sqrt((emid_intersect(1,1)- cells(i,1))^2+(emid_intersect(1,2)-cells(i,2))^2);
            if (distcell2ext < distP1out2ext) && (distemid2cell< distcell2mid) %externa cells have a shorter distance to externa border, than the closest outer pt has to that same externa boarder
                                            %and closer to externa midline than pP midline 
               pEcells(i,:) = cells(i,:);
               pPcells(i,:)= [0,0];
               if p1emid(1,3) <p2emid(1,3)
                   pEcells_on_line(i,1)= externa_unrolled((p1emid(1,3)),1)+ distP1emid2int;% does it go before or after closest midline pt?
               else %closestindex(1)>closestindex(2)
                   pEcells_on_line(i,1)= externa_unrolled((p1emid(1,3)),1)- distP1emid2int; 
               end
               if distemid2ext> distcell2ext %if the midline point has a greater distance to the outer boundary
        %than the distance from the data point to the outer point- then
        %it's positive (ie, Superficial cells are positive)
                pEcells_on_line(i,2)= distemid2cell;
               else
                   pEcells_on_line(i,2)= distemid2cell*-1;
               end
               pEcells_on_line(i,3)= pEcells_on_line(i,2)/distemid2ext;
               if pEcells_on_line(i,3) > 1
                   pEcells_on_line(i,3) = 1;
               end
               if pEcells_on_line(i,3) < -1
                   pEcells_on_line(i,3) = -1;
               end
               
            else % NOT EXTERNA CELL
                  pPcells(i,:) = cells(i,:);
                  pEcells(i,:)= [0,0];
                 distcell2outer = sqrt((cells(i,1)-out_intersect(1,1))^2+(cells(i,2)-out_intersect(1,2))^2);
                 distP1mid2int = sqrt((p1mid(1,1)-mid_intersect(1,1))^2+(p1mid(1,2)-mid_intersect(1,2))^2);
                 distP1mid2outer = sqrt((p1mid(1,1)-(outer(p1mid(1,3),1)))^2+(p1mid(1,2)-(outer(p1mid(1,3),2)))^2);
                distmidint2outer = sqrt((out_intersect(1,1)-mid_intersect(1,1))^2+(out_intersect(1,2)-mid_intersect(1,2))^2);
                 if p1mid(1,3)< p2mid(1,3)
                     pPcells_on_line(i,1) =midline_unrolled(p1mid(1,3),1)+distP1mid2int;
                 else
                     pPcells_on_line(i,1) =midline_unrolled(p1mid(1,3),1)-distP1mid2int;
                 end
                 if distmidint2outer > distcell2outer
                     pPcells_on_line(i,2) = distcell2mid;
                 else
                     pPcells_on_line(i,2) = distcell2mid *-1;
                 end
                 pPcells_on_line(i,3) = pPcells_on_line(i,2)/distP1mid2outer;
                 if pPcells_on_line(i,3)>1
                     pPcells_on_line(i,3) =1;
                 end
                 if pPcells_on_line(i,3)< -1
                     pPcells_on_line(i,3)= -1;
                 end
            end
        end
        %get rid of all the zero place holders
        c1=pPcells(:,1);
        c2=pPcells(:,2);
        e1=pEcells(:,1);
        e2=pEcells(:,2);
        col1=pPcells_on_line(:,1);
        col2=pPcells_on_line(:,2);
        col3=pPcells_on_line(:,3);
        ecl1=pEcells_on_line(:,1);
        ecl2=pEcells_on_line(:,2);
        ecl3=pEcells_on_line(:,3);
        c1(c1==0)=[];
        c2(c2==0)=[];
        e1(e1==0)=[];
        e2(e2==0)=[];
        col1(col1==0)=[];
        col2(col2==0)=[];
        col3(col3==0)=[];
        ecl1(ecl1==0)=[];
        ecl2(ecl2==0)=[];
        ecl3(ecl3==0)=[];

        pPcells=[c1 c2];
        pEcells=[e1 e2];
        pPcells_on_line =[col1 col2 col3];
        pEcells_on_line = [ecl1 ecl2 ecl3];
    end
end
%%
%find boarder intersections with midlinefind closest points to boarder points on the midline
if length(lat_dors)>2|| length(dors_med)>2 || length(lat_vp)>2 || length(med_vp)>2
    error('lines drawn for boarders between subregions (ie, lateralis-dorsalis boundary)cannot be more than two points')
else
    boarders = [lat_dors;dors_med;lat_vp;med_vp];
    boarders_on_line = zeros(4,1);
    for i = 1:4
        if boarders((i*2-1),1) == 0
            boarders_on_line(i) = 0;
        else
            [BX1]= boarders((i*2-1),1);
            [BX2]= boarders((i*2),1);
            [BY1]= boarders((i*2-1),2);
            [BY2]= boarders((i*2),2);
            BX = [BX1;BX2];
            BY = [BY1;BY2];
            [Xmp, ~, ~, closeto] = intersections(BX, BY, midpoints(:,1), midpoints(:,2));
            if (i== 1) && isempty(Xmp)
            boarders_on_line(1)=0;% if lateralis_dorsalis boarder doesn't cross, put in first midpoints and first externa midpoints
            elseif (i == 4) && isempty(Xmp)
                boarders_on_line(4)= max(midline_unrolled); % if med_vp doesn't cross, put in last midpoints
            else 
                on_line1 = midline_unrolled(floor(closeto));
                on_line2 = midline_unrolled(ceil(closeto));
                boarders_on_line(i) = on_line1+((on_line2-on_line1)*mod(closeto,1));
            end
        end
    end
end
if (med_vp(1,1) == 0) && (dors_med(1,1) == 0) && (lat_vp(1,1) == 0) % only lat_dors boarder
    [Lidx]= find(pPcells_on_line(:,1) > boarders_on_line(1));
    [Didx] = find(pPcells_on_line(:,1)< boarders_on_line(1));
    [Midx] = [];
    [Vidx] = [];
elseif (med_vp(1,1) == 0) && (lat_vp(1,1) == 0) %lat_dors + dors_med boarder
    [Lidx]= find(pPcells_on_line(:,1) > boarders_on_line(1));
    [Didx]= find(pPcells_on_line(:,1)< boarders_on_line(1)& pPcells_on_line(:,1)> boarders_on_line(2));
    [Midx] = find(pPcells_on_line(:,1)<boarders_on_line(2));
    [Vidx] =[];
elseif (med_vp(1,1) == 0) && (dors_med(1,1) == 0) %lat_dors + lat_vp boarder
    [Lidx]= find(pPcells_on_line(:,1) > boarders_on_line(1)& pPcells_on_line(:,1)<boarders_on_line(3));
    [Didx]= find(pPcells_on_line(:,1)< boarders_on_line(1));
    [Vidx] = find(pPcells_on_line(:,1)>boarders_on_line(3));
    [Midx] = [];
else %full circle
    [Lidx]= find(pPcells_on_line(:,1) > boarders_on_line(1)& pPcells_on_line(:,1)<boarders_on_line(3));
    [Didx]= find(pPcells_on_line(:,1)< boarders_on_line(1)& pPcells_on_line(:,1)> boarders_on_line(2));
    [Vidx] = find(pPcells_on_line(:,1)>boarders_on_line(3));
    [Midx] = find(pPcells_on_line(:,1)<boarders_on_line(2));
end
pPM_cells = pPcells(Midx,:);
pPM_cells_on_line = pPcells_on_line(Midx,:);
pPD_cells = pPcells(Didx,:);
pPD_cells_on_line = pPcells_on_line(Didx,:);
pPL_cells = pPcells(Lidx,:);
pPL_cells_on_line = pPcells_on_line(Lidx,:);
pPV_cells = pPcells(Vidx,:);
pPV_cells_on_line = pPcells_on_line(Vidx,:);


pPLsupdep = pPL_cells_on_line(:,3);
pPDsupdep = pPD_cells_on_line(:,3);
pPMsupdep = pPM_cells_on_line(:,3);
pPVsupdep = pPV_cells_on_line(:,3);
pEsupdep = pEcells_on_line(:,3);

super_lat = length(pPLsupdep(pPLsupdep>0));
deep_lat = length(pPLsupdep(pPLsupdep<0));
super_dors = length(pPDsupdep(pPDsupdep>0));
deep_dors = length(pPDsupdep(pPDsupdep<0));
super_med = length(pPMsupdep(pPMsupdep>0));
deep_med = length(pPMsupdep(pPMsupdep<0));
super_vp = length(pPVsupdep(pPVsupdep>0));
deep_vp = length(pPVsupdep(pPVsupdep<0));
super_PE =length(pEsupdep(pEsupdep>0));
deep_PE = length(pEsupdep(pEsupdep<0));
%%
%write to excel file Sheet1 = just totals
xlswrite(output_file, {'lateralis','dorsalis','medialis','vp','externa'},'Sheet1','B1');
xlswrite(output_file, {'superficial';'deep';'total'},'Sheet1','A2');
xlswrite(output_file,{super_lat,super_dors,super_med,super_vp,super_PE;deep_lat,deep_dors,deep_med,deep_vp,deep_PE; (super_lat+deep_lat), (super_dors+deep_dors),(super_med+deep_med),(super_vp+deep_vp),(super_PE+deep_PE)},'Sheet1','B2');

%write to excel file Sheet2 line coordinates, all set so lat-dors boarder=0
lat_adjust = boarders_on_line(1);
xlswrite(output_file,{'midline','lat cells on line','lat cells offset','lat dors boarder','dors cells on line','dors cells offset', 'dors med boarder','med cells on line', 'med cells offset','lat vp boarder','vp cells on line', 'vp cells offset','med vp boarder','externa on line','externa offset'},'Sheet2','A1');
xlswrite(output_file,(0-lat_adjust),'Sheet2','A2');
xlswrite(output_file,((max(midline_unrolled(:,1))-lat_adjust)),'Sheet2','A3');
if ~isempty(pPL_cells_on_line) && pPL_cells_on_line(1,1) ~= 0
    xlswrite(output_file,(pPL_cells_on_line(:,1)-lat_adjust),'Sheet2','B2');
    xlswrite(output_file,pPL_cells_on_line(:,3),'Sheet2','C2');
end
xlswrite(output_file,boarders_on_line(1)-lat_adjust,'Sheet2','D2');
if ~isempty(pPD_cells_on_line) && pPD_cells_on_line(1,1) ~= 0
    xlswrite(output_file,(pPD_cells_on_line(:,1)-lat_adjust),'Sheet2','E2');
    xlswrite(output_file,pPD_cells_on_line(:,3),'Sheet2','F2');
end
if boarders_on_line(2) ~= 0
    xlswrite(output_file,boarders_on_line(2)-lat_adjust,'Sheet2','G2');
end
if ~isempty(pPM_cells_on_line) && pPM_cells_on_line(1,1) ~= 0
    xlswrite(output_file,(pPM_cells_on_line(:,1)-lat_adjust),'Sheet2','H2');
    xlswrite(output_file,pPM_cells_on_line(:,3),'Sheet2','I2');
end
if boarders_on_line(3) ~= 0
    xlswrite(output_file,boarders_on_line(3)-lat_adjust,'Sheet2','J2');
end
if ~isempty(pPV_cells_on_line) && pPV_cells_on_line(1,1) ~= 0
    xlswrite(output_file,(pPV_cells_on_line(:,1)-lat_adjust),'Sheet2','K2');
    xlswrite(output_file,pPV_cells_on_line(:,3),'Sheet2','L2');
end
if boarders_on_line(4) ~= 0
    xlswrite(output_file,boarders_on_line(1)-lat_adjust,'Sheet2','M2');
end
if ~isempty(pEcells_on_line) && pEcells_on_line(1,1) ~= 0
    xlswrite(output_file,(pEcells_on_line(:,1)-lat_adjust),'Sheet2','N2');
    xlswrite(output_file,pEcells_on_line(:,3),'Sheet2','O2');
end

%write to excel file Sheet3 raw XY coordinates
xlswrite(output_file,{'midlineX','midlineY','outerX','outerY','innerX','innerY', 'externaX', 'externaY','lateralisX','lateralisY', 'dorsalisX','dorsalisY','medialisX','medialisY', 'ventrpX','ventrpY', 'boardersX','boardersY'},'Sheet3','A1');
xlswrite(output_file,midpoints(:,1),'Sheet3','A2');
xlswrite(output_file,midpoints(:,2),'Sheet3','B2');
xlswrite(output_file,outer(:,1),'Sheet3','C2');
xlswrite(output_file,outer(:,2),'Sheet3','D2');
xlswrite(output_file,inner(:,1),'Sheet3','E2');
xlswrite(output_file,inner(:,2),'Sheet3','F2');
if ~isempty(externa) && externa(1,1) ~=0
xlswrite(output_file,externa(:,1),'Sheet3','G2');
xlswrite(output_file,externa(:,2),'Sheet3','H2');
end
if ~isempty(pPL_cells) && pPL_cells(1,1) ~=0
    xlswrite(output_file,pPL_cells(:,1),'Sheet3','I2');
    xlswrite(output_file,pPL_cells(:,2),'Sheet3','J2');
end
if ~isempty(pPD_cells) && pPD_cells(1,1) ~=0
    xlswrite(output_file,pPD_cells(:,1),'Sheet3','K2');
    xlswrite(output_file,pPD_cells(:,2),'Sheet3','L2');
end
if ~isempty(pPM_cells) && pPM_cells(1,1) ~=0
    xlswrite(output_file,pPM_cells(:,1),'Sheet3','M2');
    xlswrite(output_file,pPM_cells(:,2),'Sheet3','N2');
end
if ~isempty(pPV_cells) && pPV_cells(1,1) ~=0
    xlswrite(output_file,pPV_cells(:,1),'Sheet3','O2');
    xlswrite(output_file,pPV_cells(:,2),'Sheet3','P2');
end
if ~isempty(boarders) && boarders(1,1) ~=0
xlswrite(output_file,boarders(:,1),'Sheet3','Q2');
xlswrite(output_file,boarders(:,2),'Sheet3','R2');
end

figure
plot(outer(:,1),outer(:,2),'k');
hold on
plot(inner(:,1),inner(:,2),'k');
plot(midpoints(:,1),midpoints(:,2),'--','Color','k');
if ~isempty(externa)&& externa(1,1) ~=0
    plot(externa(:,1),externa(:,2),'--','Color',[0, .4,0]);
end
plot(lat_dors(:,1),lat_dors(:,2), 'y');
plot(dors_med(:,1),dors_med(:,2), 'y');
plot(lat_vp(:,1),lat_vp(:,2),'y');
plot(med_vp(:,1),med_vp(:,2),'y');
if ~isempty(pPL_cells)&& pPL_cells(1,1) ~=0
    plot(pPL_cells(:,1),pPL_cells(:,2),'s','Color','m');
end
if ~isempty(pPD_cells)&& pPD_cells(1,1) ~=0
    plot(pPD_cells(:,1),pPD_cells(:,2),'d','Color','b');
end
if ~isempty(pPM_cells)&& pPM_cells(1,1) ~=0
    plot(pPM_cells(:,1),pPM_cells(:,2),'s','Color','r');
end
if ~isempty(pPV_cells)&& pPV_cells(1,1) ~=0
    plot(pPV_cells(:,1),pPV_cells(:,2),'d','Color','c');
end
if ~isempty(pEcells)&& pEcells(1,1) ~=0
    plot(pEcells(:,1),pEcells(:,2),'*','Color','g');
end
%%

figure
plot([0;0],[(-lat_adjust);(max(midline_unrolled(:,1))-lat_adjust)],'k');
hold on
plot([-1;1],[boarders_on_line(1,1)-lat_adjust;boarders_on_line(1,1)-lat_adjust],'y');
plot([-1;1],[boarders_on_line(2,1)-lat_adjust;boarders_on_line(2,1)-lat_adjust],'y');
plot([-1;1],[boarders_on_line(3,1)-lat_adjust;boarders_on_line(3,1)-lat_adjust],'y');
plot([-1;1],[boarders_on_line(4,1)-lat_adjust;boarders_on_line(4,1)-lat_adjust],'y');
if ~isempty(pPL_cells_on_line) && pPL_cells_on_line(1,1) ~= 0
    plot(pPL_cells_on_line(:,3),pPL_cells_on_line(:,1)-lat_adjust,'s', 'Color','m')
end
if ~isempty(pPD_cells_on_line) && pPD_cells_on_line(1,1) ~= 0
    plot(pPD_cells_on_line(:,3),pPD_cells_on_line(:,1)-lat_adjust,'d', 'Color','b');
end
if ~isempty(pPM_cells_on_line) && pPM_cells_on_line(1,1) ~= 0
    plot(pPM_cells_on_line(:,3),pPM_cells_on_line(:,1)-lat_adjust,'s', 'Color','r');
end
if ~isempty(pPV_cells_on_line) && pPV_cells_on_line(1,1) ~= 0
    plot(pPV_cells_on_line(:,3),pPV_cells_on_line(:,1)-lat_adjust,'d', 'Color','c');
end

figure
plot([0;0],[(-lat_adjust);(max(midline_unrolled(:,1))-lat_adjust)],'k');
hold on
plot([-.1;.1],[boarders_on_line(1,1)-lat_adjust;boarders_on_line(1,1)-lat_adjust],'y');
plot([-.1;.1],[boarders_on_line(2,1)-lat_adjust;boarders_on_line(2,1)-lat_adjust],'y');
plot([-.1;.1],[boarders_on_line(3,1)-lat_adjust;boarders_on_line(3,1)-lat_adjust],'y');
plot([-.1;.1],[boarders_on_line(4,1)-lat_adjust;boarders_on_line(4,1)-lat_adjust],'y');
if ~isempty(pPL_cells_on_line) && pPL_cells_on_line(1,1) ~= 0
    plot((zeros(length(pPL_cells_on_line(:,3)))),pPL_cells_on_line(:,1)-lat_adjust,'.', 'Color','m')
end
if ~isempty(pPD_cells_on_line) && pPD_cells_on_line(1,1) ~= 0
    plot((zeros(length(pPD_cells_on_line(:,3)))),pPD_cells_on_line(:,1)-lat_adjust,'.', 'Color','b');
end
if ~isempty(pPM_cells_on_line) && pPM_cells_on_line(1,1) ~= 0
    plot((zeros(length(pPM_cells_on_line(:,3)))),pPM_cells_on_line(:,1)-lat_adjust,'.', 'Color','r');
end
if ~isempty(pPV_cells_on_line) && pPV_cells_on_line(1,1) ~= 0
    plot((zeros(length(pPV_cells_on_line(:,3)))),pPV_cells_on_line(:,1)-lat_adjust,'.', 'Color','c');
end

    