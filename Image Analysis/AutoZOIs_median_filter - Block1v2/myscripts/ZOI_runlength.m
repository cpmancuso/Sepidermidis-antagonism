function zoi_roi = ZOI_runlength(norm_int,roi);
runlen = 5;
sz = numel(norm_int);
temp_roi = reshape((roi&norm_int<0),1,sz); %logical vector
temp_int = reshape(norm_int,1,sz);
roi_pos = find(~[temp_roi 0]); %all regions beyond edge and <0, padded the end;

% ZOI_subarea(i,n) = -trapz(norm_int3(i,n,zoi_roi))

[~, grpidxs] = find(diff(roi_pos)>runlen); %indexes of of roi idxs with runlength>5

if grpidxs
    for n=1:numel(grpidxs) %loop through run lengths
        tempidx=grpidxs(n);
        ZOI_subarea(n)=-trapz(temp_int(roi_pos(tempidx)+1:roi_pos(tempidx+1)-1));
    end
    [~,max_run_len] = max(ZOI_subarea);
    max_run_len_idx = grpidxs(max_run_len);
    zoi_roi = false(size(roi));
    zoi_roi(roi_pos(max_run_len_idx)+1:roi_pos(max_run_len_idx+1)-1)=true;

else
    zoi_roi = false(size(roi));
end



