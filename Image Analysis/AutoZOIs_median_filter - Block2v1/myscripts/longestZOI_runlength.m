function zoi_roi = ZOI_runlength(norm_int,roi);

sz = numel(norm_int);
temp_roi = reshape((roi&norm_int<0),1,sz); %logical vector
temp_int = reshape(norm_int,1,sz);
roi_pos = find(~[temp_roi 0]); %all regions beyond edge and <0, padded the end;
[~, grpidx] = max(diff(roi_pos)); %index of largest group of roi idxs
zoi_roi = false(size(roi));
zoi_roi(roi_pos(grpidx)+1:roi_pos(grpidx+1)-1)=true;

