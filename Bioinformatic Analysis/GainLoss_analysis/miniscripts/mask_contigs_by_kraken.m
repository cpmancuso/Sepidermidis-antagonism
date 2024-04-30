clc
clear all
close all

fid = fopen('sepi_clade_49_krakSeq.txt');
tline = fgetl(fid);
taxes = [];
num=0;
while ischar(tline)
    num=num+1;
    test = strsplit(tline);
    taxes(num)=str2num(test{3});
    tline = fgetl(fid);
end
fclose(fid);
%%
unique_taxes = unique(taxes);
unique_taxes = unique_taxes(unique_taxes>2);
for n=1:numel(unique_taxes)
    fileID = fopen('out.xml','w');
    baseURL = ['http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=' num2str(unique_taxes(n)) '&retmode=xml'];
    searchReport = webread(baseURL);
    fprintf(fileID,searchReport);
    ncbi = readstruct('out.xml');
    disp(ncbi.Taxon.LineageEx.Taxon)
    fclose(fileID)
end
