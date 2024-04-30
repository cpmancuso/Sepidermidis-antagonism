function genbank_to_CDS(dir_genome,path_fasta,path_gb)

%% Initialize

% Load fasta
fastaFile = fastaread(path_fasta);
genbankAll = genbankread(path_gb); 

[~, ~, ~, ScafNames]= get_genome_stats(path_fasta);

%%
fprintf(1,'Reading in gb file...\n');

clear CDS

for i=1:length(ScafNames)
    %genome fasta must have as many contigs as genbank file
    if numel(genbankAll(:)) > length(ScafNames)
        error('scaffold number must equal or exceed the number of scaffolds in the genbank file. ')
    end
    %skip contigs without annotations
    if numel(genbankAll(:))<i
        fprintf(1,'Warning: Genbank file may have fewer contigs than exist in genome fasta.\n');
        continue
    else
        genbankContig = genbankAll(i);
    end
    % load genbank
    if ~isfield(genbankContig,'Sequence')
        fprintf(1,'Warning: Genbank file may have no CDS/Feature or first line of the genbank file may not be at least 79 characters including whitespaces.\n');
        genbankContig.Sequence = lower(fastaFile(i).Sequence);
    elseif isempty(genbankContig.Sequence)% check if sequence is empty
        genbankContig.Sequence = lower(fastaFile(i).Sequence);
    end
    % get CDS
    if isfield(genbankContig,'CDS') || isfield(genbankContig,'Features')
        genes = locustag_from_text(genbankContig.CDS) ;
        genes = div_add_in_nonCDS(genes, genbankContig.Features);
        
        if isfield(genbankContig, 'Sequence')
            CDS{i} = parse_all_locations_gb(genes, genbankContig.Sequence);
        else
            CDS{i} = parse_all_locations_gb(genes, char([fastaFile.Sequence]+32)) ;  %also reverses strands in this function, sorts by position
        end
        
        %sort by position on genome
        if ~isempty(CDS{i})
            [~,sortedpositions]=sort([CDS{i}.loc1]);
            CDS{i}=CDS{i}(sortedpositions);
        end
    end
end

%% Save CDS

save([dir_genome '/cds_sorted.mat'],'CDS')

end
