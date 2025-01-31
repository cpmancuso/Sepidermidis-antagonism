function make_clade_coverage_plot(SampleNames, this_clade_name, contig_start_indices, contig_lengths, median_copy_by_contig, norm_median_cov_by_contig,mask_copynum,mask_spades,mask_length,mask_custom,mask,copynum_cutoff,spades_cutoff,length_cutoff, spades_contig_coverage, dir_save_coverage)
    %% Plot Coverage By Sample and By Contig
    f=figure('Position',[0 0 1200 800],'visible','off');

    subplot(3,3,1:3)
    mksz=20;
    alpha=0.5;
    num_samples = numel(SampleNames);
    num_contigs = numel(contig_start_indices);

    for i=1:num_samples
        aplot=stairs(contig_start_indices,norm_median_cov_by_contig(i,:));
%         adata = dataTipTextRow('Sample',repmat({SampleNames{i}},size(aplot.XData)));
%         aplot.DataTipTemplate.DataTipRows(end+1) = adata;
        hold on
    end
    ylim([-2 4])
    if numel(SampleNames)<16
        legend(SampleNames,'location','eastoutside','Interpreter', 'none')
    end
    xlabel('Genome Position')
    ylabel('Normalized Coverage (standard deviations)')
    
    %% Plot spades coverage by length
    subplot(3,3,[4 7])
    scatter(contig_lengths(mask),spades_contig_coverage(mask),mksz,'k','filled','MarkerFaceAlpha',alpha)
    hold on
    % low kmer coverage, kept
    scatter(contig_lengths((mask_copynum&~mask_spades)&mask_length),spades_contig_coverage((mask_copynum&~mask_spades)&mask_length),mksz,'b','filled','MarkerFaceAlpha',alpha)
    % small contigs, filtered
    scatter(contig_lengths(~mask_length),spades_contig_coverage(~mask_length),mksz,'y','filled','MarkerFaceAlpha',alpha)
    % low both, filtered
    scatter(contig_lengths(mask_length&~mask),spades_contig_coverage(mask_length&~mask),mksz,'r','filled','MarkerFaceAlpha',alpha)

    plot([100 1000000],[spades_cutoff,spades_cutoff],'r--')
    plot([length_cutoff length_cutoff],[0.01 100000],'r--')
    ylim([0.01 100000])
    xlim([100 1000000])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    legend('Included Contig','Included Contig, Coverage Disagree','Masked Contig, Short','Masked Contig, Low Coverage','location','northoutside')
    xlabel('Contig Length')
    ylabel('Contig SPAdes k-mer coverage')
    
    %% Plot copy number by length
    subplot(3,3,[5 8])

    copy_number_to_plot=max(median_copy_by_contig)+0.001;
    scatter(contig_lengths(mask),copy_number_to_plot(mask),mksz,'k','filled','MarkerFaceAlpha',alpha)
    hold on
    % low kmer coverage, kept
    scatter(contig_lengths((mask_copynum&~mask_spades)&mask_length),copy_number_to_plot((mask_copynum&~mask_spades)&mask_length),mksz,'b','filled','MarkerFaceAlpha',alpha)
    % small contigs, filtered
    scatter(contig_lengths(~mask_length),copy_number_to_plot(~mask_length),mksz,'y','filled','MarkerFaceAlpha',alpha)
    % low both, filtered
    scatter(contig_lengths(mask_length&~mask),copy_number_to_plot(mask_length&~mask),mksz,'r','filled','MarkerFaceAlpha',alpha)

    plot([100 1000000],[copynum_cutoff,copynum_cutoff],'r--')
    plot([length_cutoff length_cutoff],[0.001 1000],'r--')
    ylim([0.001 1000])
    xlim([100 1000000])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    legend('Included Contig','Included Contig, Coverage Disagree','Masked Contig, Short','Masked Contig, Low Coverage','location','northoutside')
    xlabel('Contig Length')
    ylabel('Maximum Copy Number in any Sample')

    %% Plot copy number by spades coverage
    subplot(3,3,[6 9])
    scatter(spades_contig_coverage(mask),copy_number_to_plot(mask),mksz,'k','filled','MarkerFaceAlpha',alpha)
    hold on
    % low kmer coverage but high copy, kept
    scatter(spades_contig_coverage((mask_copynum&~mask_spades)&mask_length),copy_number_to_plot((mask_copynum&~mask_spades)&mask_length),mksz,'b','filled','MarkerFaceAlpha',alpha)
    % small contigs, filtered
    scatter(spades_contig_coverage(~mask_length),copy_number_to_plot(~mask_length),mksz,'y','filled','MarkerFaceAlpha',alpha)
    % low both, filtered
    scatter(spades_contig_coverage(mask_length&~mask),copy_number_to_plot(mask_length&~mask),mksz,'r','filled','MarkerFaceAlpha',alpha)
    plot([0.01 100000],[copynum_cutoff,copynum_cutoff],'r--')
    plot([spades_cutoff,spades_cutoff],[0.001 1000],'r--')
    ylim([0.001 1000])
    xlim([0.01 100000])
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    legend('Included Contig','Included Contig, Coverage Disagree','Masked Contig, Short','Masked Contig, Low Coverage','location','northoutside')
    xlabel('Contig SPAdes k-mer coverage')
    ylabel('Maximum Copy Number in any Sample')
    
    old_length = round(sum(contig_lengths)./1000000,1);
    new_length = round(sum(contig_lengths(mask))./1000000,1);
    
    temp_title = [this_clade_name '; Samples ' num2str(num_samples) '; included ' num2str(sum(mask)) '/' num2str(num_contigs) ' contigs; Length ' num2str(new_length) '/' num2str(old_length) 'Mb'];
    if verLessThan('matlab','9.5')
        suptitle(strrep(temp_title,'_','-'))
    else
        sgtitle(temp_title,'Interpreter', 'none')
    end    
%     saveas(f,[dir_save_coverage '/coverage_' this_clade_name '.fig'])
    print([dir_save_coverage '/coverage_' this_clade_name '.png'],'-dpng')

end