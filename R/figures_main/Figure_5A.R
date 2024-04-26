
.Figure_5A <-function(all_cutoff_passed_tumor12=TRUE, sel_study='na') {
  
  point_cols <-c(pal_simpsons(palette = c("springfield"), alpha = 1)(10),pal_simpsons(palette = c("springfield"), alpha = 0.3)(10))
  
  require(ggthemes)
  
  theme_base1 <-
    theme_bw() +
    theme(plot.title = element_text(size = 12),
          plot.margin = unit(c(0,0.5,0,0.5),"lines"),
          panel.margin = unit(0.3, "lines"),
          axis.text.x = element_text(size=12,angle=50,colour = 'black', face = 'italic',hjust=1,vjust=1),
          axis.text.x.top = element_text(size=12,angle=50,colour = 'black',  hjust = 0, vjust = 0, face = 'bold.italic'),
          axis.text.x.bottom = element_text(size=12,angle=50,colour = 'black',  hjust = 1, vjust = 1, face = 'bold.italic'),
          axis.text.y = element_text(size=12,angle=0,colour = 'black', face = 'italic',hjust=1,vjust=0.5),
          axis.title.x = element_text(size=13,colour = 'black', face = 'bold'),
          axis.title.y = element_text(size=13,colour = 'black', face = 'bold'),
          legend.title  = element_text(size=13,face="bold.italic"),
          legend.text = element_text(size=12),
          legend.key.size = unit(1, "line"),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.background = element_blank(),
          axis.ticks.x =element_blank(),
          axis.ticks.y =element_blank(),
          axis.ticks.length = unit(0,"null"),
          axis.ticks=element_line(colour="white"),
          strip.background = element_rect(fill="white",colour=NA),
          strip.text.x = element_text(size = 13, colour="black",face="bold",angle = 0, vjust=0,hjust=0.5),
          strip.text.y = element_text(size = 13, colour="black",face="italic",angle = 0,hjust=0))
  
  stopifnot(all_cutoff_passed_tumor12)
  
  #' Load data
  data <-readRDS(file.path(output_basedir, 'nat_genet_github_1st',
                           sprintf('data_aaSuite_somatic_ms_kataegis_amplicons_overlap-%s.rds','all_cutoff_passed_tumor12')))
  
  
  output_basename <-file.path(output_basedir, 'nat_genet_github_1st',
                              sprintf('Figure_5A-%s', 'all_cutoff_passed_tumor12'))
  
  dt_panecdnams_amplicon_similarity <-data$dt_panecdnams_amplicon_similarity
  dt_panecdnams_amplicons <-data$dt_panecdnams_amplicons
  panecdnams_kataegis_calls <-data$panecdnams_kataegis_calls
  
  sel_samples             <-data$sel_samples
  sel_patients            <-data$sel_patients
  
  d_amp1    <-data$d_amp1
  d_amp2    <-data$d_amp2
  
  d_amp1_kataegis_vafs    <-data$d_amp1_kataegis_vafs
  d_amp2_kataegis_vafs    <-data$d_amp2_kataegis_vafs
  
  d_snvs_sample1_with_amplicons <-data$d_snvs_sample1_with_amplicons
  d_snvs_sample2_with_amplicons <-data$d_snvs_sample2_with_amplicons
  
  
  #' #############################################################
  #' 
  #' Plot VAFs of clustered and unclustered mutations in shared ecDNAs and private ecDNAs
  #' 
  #' #############################################################
  #' ======================================================
  #' Note that shared amplicons should have shared mutations only, and private amplicons have private mutations only
  #' ======================================================
  .make_data_tmp <-function(A_amplicon="2",B_amplicon="1") {
    stopifnot(A_amplicon!=B_amplicon)
    stopifnot(identical(sort(c(A_amplicon, B_amplicon)), sort(c('1','2'))))
    
    if (A_amplicon=='2') {
      
      d_amp_kataegis_vafs           <-d_amp2_kataegis_vafs
      d_snvs_sample_with_amplicons  <-d_snvs_sample2_with_amplicons
      d_amp                         <-d_amp2 %>% 
        dplyr::select(-one_of('Intervals')) %>%
        left_join(dt_panecdnams_amplicons %>% dplyr::select(amplicon_barcode, Intervals), by='amplicon_barcode')
      
      d_snvs_sample_with_amplicons %>% group_by(aa_barcode) %>% summarise(n=n()) %>% ungroup() %>% summarise(mean_snv=mean(n), median_snv=median(n)) %>% as.data.frame()
      
      stopifnot(identical(sort(as.character(unique(d_amp$amplicon2_has_similar_amplicon1))), sort(c('Yes','No'))))
      d_amp <-d_amp %>% mutate(is_shared_amplicon=ifelse(amplicon2_has_similar_amplicon1=='Yes','Shared_amplicon','Private_amplicon'))
      stopifnot(all(!is.element(filter(d_amp, is_shared_amplicon=='Private_amplicon')$amplicon_barcode, 
                                filter(dt_panecdnams_amplicon_similarity, 
                                       SimScorePvalue<0.05, 
                                       amplicon_classification1==amplicon_classification2,
                                       amplicon_classification1!='No amp/Invalid',
                                       amplicon_classification2!='No amp/Invalid')$Amp2)))
      stopifnot(all(is.element(filter(d_amp, is_shared_amplicon=='Shared_amplicon')$amplicon_barcode, 
                               filter(dt_panecdnams_amplicon_similarity, 
                                      SimScorePvalue<0.05, 
                                      amplicon_classification1==amplicon_classification2,
                                      amplicon_classification1!='No amp/Invalid',
                                      amplicon_classification2!='No amp/Invalid')$Amp2)))
      
      line_color <-'#E08214'
      #   mean_snv median_snv
      # 1 25374.84      12406
    } else if (A_amplicon=='1') {
      
      d_amp_kataegis_vafs           <-d_amp1_kataegis_vafs
      d_snvs_sample_with_amplicons  <-d_snvs_sample1_with_amplicons
      d_amp                         <-d_amp1 %>%
        dplyr::select(-one_of('Intervals')) %>%
        left_join(dt_panecdnams_amplicons %>% dplyr::select(amplicon_barcode, Intervals), by='amplicon_barcode')
      
      d_snvs_sample_with_amplicons %>% group_by(aa_barcode) %>% summarise(n=n()) %>% ungroup() %>% summarise(mean_snv=mean(n), median_snv=median(n)) %>% as.data.frame()
      
      stopifnot(identical(sort(as.character(unique(d_amp$amplicon1_has_similar_amplicon2))), sort(c('Yes','No'))))
      d_amp <-d_amp %>% mutate(is_shared_amplicon=ifelse(amplicon1_has_similar_amplicon2=='Yes','Shared_amplicon','Private_amplicon'))
      stopifnot(all(!is.element(filter(d_amp, is_shared_amplicon=='Private_amplicon')$amplicon_barcode, 
                                filter(dt_panecdnams_amplicon_similarity, 
                                       SimScorePvalue<0.05, 
                                       amplicon_classification1==amplicon_classification2,
                                       amplicon_classification1!='No amp/Invalid',
                                       amplicon_classification2!='No amp/Invalid')$Amp1)))
      stopifnot(all(is.element(filter(d_amp, is_shared_amplicon=='Shared_amplicon')$amplicon_barcode, 
                               filter(dt_panecdnams_amplicon_similarity, 
                                      SimScorePvalue<0.05, 
                                      amplicon_classification1==amplicon_classification2,
                                      amplicon_classification1!='No amp/Invalid',
                                      amplicon_classification2!='No amp/Invalid')$Amp1)))
      
      line_color <-'#8073AC'
      
    }
    
    #' ======================================================
    #' Check input file
    #' ======================================================
    stopifnot(all(!duplicated(d_amp_kataegis_vafs$tag)))
    stopifnot(all(!duplicated(d_snvs_sample_with_amplicons$tag)))
    stopifnot(all(is.element(d_amp_kataegis_vafs$tag, d_snvs_sample_with_amplicons$tag)))
    stopifnot(identical(sort(unique(d_snvs_sample_with_amplicons$is_shared)), sort(c("Private","Shared"))))
    stopifnot(all(!is.na(d_amp_kataegis_vafs$kataegis_chr)))
    stopifnot(all(is.element(d_amp$aa_barcode, d_snvs_sample_with_amplicons$aa_barcode)))
    
    #' Group mutations into 
    #' clustered mutations
    #' non clustered mutations
    stopifnot(all(!duplicated(d_amp_kataegis_vafs$tag)))
    stopifnot(all(!duplicated(d_snvs_sample_with_amplicons$tag)))
    stopifnot(all(is.element(d_amp_kataegis_vafs$tag, d_snvs_sample_with_amplicons$tag)))
    
    stopifnot(all(!duplicated(d_amp$amplicon_barcode)))
    stopifnot(all(is.element(d_amp_kataegis_vafs$amplicon_barcode, d_amp$amplicon_barcode)))
    
    d_clustered <-d_amp_kataegis_vafs %>%
      dplyr::select(-one_of(c('amplicon_has_oncogene','is_shared_amplicon'))) %>%
      left_join(d_amp %>% dplyr::select(amplicon_barcode, amplicon_has_oncogene, is_shared_amplicon), by='amplicon_barcode') %>%
      dplyr::select(aa_barcode, amplicon_barcode, amplicon_classification, source, study, is_shared, snv_vaf, amplicon_has_oncogene, mut_tag, tag, is_shared_amplicon) %>% 
      dplyr::mutate(is_clustered_mutations='clustered_mutations')
    d_nonclustered <-d_snvs_sample_with_amplicons %>% 
      dplyr::filter(!is.element(tag, d_clustered$tag)) %>%
      dplyr::mutate(snv_vaf=vaf) %>%
      dplyr::select(aa_barcode, source, study, is_shared, snv_vaf, mut_tag, tag, chrom, start, end) %>% 
      dplyr::mutate(is_clustered_mutations='non_clustered_mutations')
    
    #' Find non clustered mutations overlapping with ecDNA
    d_nonclustered_overlapping_w_ecdna_intervals <-ddply(filter(d_amp, amplicon_classification=='ecDNA'), c('aa_barcode','amplicon_barcode','Intervals','amplicon_classification','is_shared_amplicon'), function(x) {
      xx <-unlist(strsplit(x$Intervals, split=","))
      y  <-data.frame(loci=xx) %>%
        dplyr::mutate(chr=gsub("([0-9XYMT]+):([0-9]+)-([0-9]+)","\\1",loci),
                      start=gsub("([0-9XYMT]+):([0-9]+)-([0-9]+)","\\2",loci),
                      end=gsub("([0-9XYMT]+):([0-9]+)-([0-9]+)","\\3",loci)) %>%
        dplyr::mutate(chk=paste0(chr,':', start, '-', end)==loci)
      stopifnot(all(y$chk))
      Xamplicon_intervals <-y %>% dplyr::select(-one_of('chk'))
      
      #' Get non clustered mutations from this aa barcode
      Xd_nonclustered <-d_nonclustered %>% filter(aa_barcode==x$aa_barcode)
      
      #' Intersect Xd_nonclustered and Xamplicon_intervals
      require(GenomicRanges)
      Xamplicon_intervals_gr  <-makeGRangesFromDataFrame(Xamplicon_intervals, keep.extra.columns = T)
      Xd_nonclustered_gr      <-makeGRangesFromDataFrame(Xd_nonclustered, keep.extra.columns = T)
      ovs <-findOverlaps(Xd_nonclustered_gr, Xamplicon_intervals_gr)
      Xd_nonclustered_gr <-Xd_nonclustered_gr[queryHits(ovs)]
      as.data.frame(Xd_nonclustered_gr)
    })
    
    #' Make d_clustered_nonclustered_mutations_overlapping_w_ecdna_intervals ##########
    #' Add amplicon_has_oncogene
    d <-bind_rows(d_clustered %>%
                    dplyr::select(aa_barcode, amplicon_barcode, amplicon_classification, source, study, is_shared, snv_vaf, mut_tag, tag, is_clustered_mutations, is_shared_amplicon) %>%
                    dplyr::filter(amplicon_classification=='ecDNA'),
                  d_nonclustered_overlapping_w_ecdna_intervals %>%
                    dplyr::select(aa_barcode, amplicon_barcode, amplicon_classification, source, study, is_shared, snv_vaf, mut_tag, tag, is_clustered_mutations, is_shared_amplicon) %>%
                    dplyr::filter(amplicon_classification=='ecDNA'))
    stopifnot(all(!is.na(d$amplicon_barcode)))
    stopifnot(all(is.element(d$amplicon_barcode, d_amp$amplicon_barcode)))
    stopifnot(all(!duplicated(d_amp$amplicon_barcode)))
    stopifnot(all(!is.na(d_amp$amplicon_has_oncogene)))
    stopifnot(identical(unique(as.character(d$amplicon_classification)),'ecDNA'))
    stopifnot(identical(sort(unique(d_amp$amplicon_has_oncogene)), sort(c('Has_oncogene','No_oncogene'))))
    d <-d %>%
      dplyr::select(-one_of('amplicon_has_oncogene')) %>%
      left_join(d_amp %>% dplyr::select(amplicon_barcode, amplicon_has_oncogene), by='amplicon_barcode')
    
    stopifnot(identical(sort(unique(as.character(d$amplicon_has_oncogene))), sort(c('Has_oncogene','No_oncogene'))))
    d <-d %>% mutate(amplicon_has_oncogene=factor(amplicon_has_oncogene, levels=c('Has_oncogene','No_oncogene')))
    
    stopifnot(identical(sort(unique(as.character(d$is_shared))), sort(c('Shared','Private'))))
    d <-d %>% mutate(is_shared=factor(is_shared, levels=c('Shared','Private')))
    
    stopifnot(identical(sort(unique(as.character(d$is_clustered_mutations))), sort(c('clustered_mutations','non_clustered_mutations'))))
    d <-d %>% mutate(is_clustered_mutations=factor(is_clustered_mutations, levels=c('clustered_mutations','non_clustered_mutations')))
    
    stopifnot(identical(sort(unique(as.character(d$is_shared_amplicon))), sort(c('Shared_amplicon','Private_amplicon'))))
    d <-d %>% mutate(is_shared_amplicon=factor(is_shared_amplicon, levels=c('Shared_amplicon','Private_amplicon')))
    
    table(d[,c('is_shared','is_shared_amplicon')])
    #                   is_shared_amplicon
    # is_shared Shared_amplicon Private_amplicon
    #   Shared             1363              473
    #   Private            5743             3988
    
    d <-d %>% 
      mutate(A_amplicon=A_amplicon, 
             B_amplicon=B_amplicon)
    
    return(d)
  }
  
  d12 <-.make_data_tmp(A_amplicon="1",B_amplicon="2")
  d21 <-.make_data_tmp(A_amplicon="2",B_amplicon="1")
  
  stopifnot(all(!duplicated(d12$tag)))
  stopifnot(all(!duplicated(d21$tag)))
  
  stopifnot(all(!duplicated(dt_panecdnams_amplicon_similarity$Amp1)))
  stopifnot(all(!duplicated(dt_panecdnams_amplicon_similarity$Amp2)))
  stopifnot(all(is.element(filter(d12, is_shared_amplicon=='Shared_amplicon')$amplicon_barcode, dt_panecdnams_amplicon_similarity$Amp1)))
  stopifnot(all(is.element(filter(d21, is_shared_amplicon=='Shared_amplicon')$amplicon_barcode, dt_panecdnams_amplicon_similarity$Amp2)))
  stopifnot(all(!duplicated(dt_panecdnams_amplicon_similarity$amp_pair)))
  
  filtered_sim_amp_pairs <-filter(dt_panecdnams_amplicon_similarity, 
                                  SimScorePvalue<0.05, 
                                  amplicon_classification1==amplicon_classification2,
                                  amplicon_classification1!='No amp/Invalid',
                                  amplicon_classification2!='No amp/Invalid')
  dt_shared_amp_snvs <-ddply(filtered_sim_amp_pairs,'amp_pair', function(Xamp_pair) {
    # Xamp_pair <-filter(dt_panecdnams_amplicon_similarity, is.element(Amp1, d12$amplicon_barcode))[1,]
    # Xamp_pair <-filter(dt_panecdnams_amplicon_similarity,
    #                    amp_pair=='0bcfd127-838e-4670-8d4a-b028dce79774__01c342ec-b7a7-47ff-9892-35d32e026a64-amplicon1--e6bd6240-68a5-4f89-b7a0-f86353725bfa__01c342ec-b7a7-47ff-9892-35d32e026a64-amplicon1')
    # Xamp_pair <-filter(dt_panecdnams_amplicon_similarity,amp_pair=='CPCT02010894T__CPCT02010894R-amplicon5--CPCT02010894TII__CPCT02010894R-amplicon5')
    
    cat('\n',Xamp_pair$amp_pair,'\n')
    
    Xd12_shared_amp_snvs <-d12 %>% filter(is.element(amplicon_barcode, Xamp_pair$Amp1))
    Xd21_shared_amp_snvs <-d21 %>% filter(is.element(amplicon_barcode, Xamp_pair$Amp2))
    
    if (nrow(Xd12_shared_amp_snvs)==0) {
      return(NULL)
    }
    if (nrow(Xd21_shared_amp_snvs)==0) {
      return(NULL)
    }
    
    shared_mut_tags <-intersect(Xd12_shared_amp_snvs$mut_tag, Xd21_shared_amp_snvs$mut_tag)
    Xd12_shared_amp_snvs <-Xd12_shared_amp_snvs %>% filter(is.element(mut_tag, shared_mut_tags))
    Xd21_shared_amp_snvs <-Xd21_shared_amp_snvs %>% filter(is.element(mut_tag, shared_mut_tags))
    
    if (length(shared_mut_tags)>0) {
      stopifnot(identical(unique(as.character(unique(Xd12_shared_amp_snvs$is_shared))), 'Shared'))
      stopifnot(identical(unique(as.character(unique(Xd21_shared_amp_snvs$is_shared))), 'Shared'))
      
    }
    
    bind_rows(Xd12_shared_amp_snvs, Xd21_shared_amp_snvs)
    
  })
  stopifnot(identical(as.character(unique(unique(dt_shared_amp_snvs$is_shared))),'Shared'))
  stopifnot(identical(as.character(unique(unique(dt_shared_amp_snvs$is_shared_amplicon))),'Shared_amplicon'))
  chk <-dt_shared_amp_snvs %>%
    group_by(A_amplicon, B_amplicon, tag) %>%
    summarise(n=n()) %>%
    ungroup()
  stopifnot(identical(as.numeric(unique(chk$n)), as.numeric(1)))
  
  
  dt_private_amp_snvs<-bind_rows(filter(d12, is_shared_amplicon=='Private_amplicon', is_shared=='Private'),
                                 filter(d21, is_shared_amplicon=='Private_amplicon', is_shared=='Private'))
  
  dt <-bind_rows(dt_shared_amp_snvs,
                 dt_private_amp_snvs) %>%
    mutate(AB_amplicon_pair=paste0(A_amplicon,'_vs_', B_amplicon))
  stopifnot(identical(as.character(unique(sort(dt$AB_amplicon_pair))), sort(c('1_vs_2','2_vs_1'))))
  dt <-dt %>% mutate(AB_amplicon_pair=factor(AB_amplicon_pair, levels=c('1_vs_2','2_vs_1')))
  chk <-dt %>%
    group_by(AB_amplicon_pair, tag) %>%
    summarise(n=n()) %>%
    ungroup()
  stopifnot(identical(as.numeric(unique(chk$n)), as.numeric(1)))
  
  dt <-dt %>%
    group_by(AB_amplicon_pair,tag) %>%
    mutate(cl=paste0(is_clustered_mutations,'-', is_shared_amplicon)) %>%
    ungroup()
  
  table(dt[,c('is_shared_amplicon','is_shared')])
  #                   is_shared
  # is_shared_amplicon Shared Private
  #   Shared_amplicon    2708       0
  #   Private_amplicon      0    5453
  
  table(dt[,c('is_clustered_mutations','is_shared_amplicon')])
  # is_shared_amplicon
  # is_clustered_mutations    Shared_amplicon Private_amplicon
  # clustered_mutations                 413              600
  # non_clustered_mutations            2295             4853
  
  stopifnot(identical(sort(as.character(unique(dt$cl))), 
                      sort(c('clustered_mutations-Private_amplicon', 
                             'non_clustered_mutations-Private_amplicon',
                             'clustered_mutations-Shared_amplicon',
                             'non_clustered_mutations-Shared_amplicon'))))
  dt <-dt %>%
    mutate(cl=factor(cl, levels=c('clustered_mutations-Private_amplicon', 'non_clustered_mutations-Private_amplicon','clustered_mutations-Shared_amplicon','non_clustered_mutations-Shared_amplicon')))
  
  dt_summary_count <-dt %>%
    group_by(AB_amplicon_pair, amplicon_classification, amplicon_has_oncogene, is_clustered_mutations, is_shared_amplicon) %>%
    summarise(n_aa_barcode=length(unique(aa_barcode)),
              n_mut=n(),
              mean_VAF=mean(snv_vaf),
              median_VAF=median(snv_vaf)) %>%
    ungroup()
  dt_summary_count2 <-dt %>%
    group_by(AB_amplicon_pair, amplicon_classification, cl) %>%
    summarise(n_aa_barcode=length(unique(aa_barcode)),
              n_mut=n(),
              mean_VAF=mean(snv_vaf),
              median_VAF=median(snv_vaf)) %>%
    ungroup()
  
  my_comparisons <- list( c('clustered_mutations-Shared_amplicon', 'clustered_mutations-Private_amplicon'), 
                          c('clustered_mutations-Shared_amplicon', 'non_clustered_mutations-Shared_amplicon'), 
                          c('clustered_mutations-Shared_amplicon', 'non_clustered_mutations-Private_amplicon'))
  
  my_comparisons <-NULL
  n<-0
  for (ix in 1:length(levels(dt$cl))) {
    for (iy in 1:length(levels(dt$cl))) {
      if (ix>iy) {
        n <-n+1
        my_comparisons[[n]]<-c(levels(dt$cl)[ix],levels(dt$cl)[iy])
      }
    }
  }
  
  p <-
    ggboxplot(dt,
              x = "cl", y = "snv_vaf",
              color = "AB_amplicon_pair",
              palette = 'npg',
              add.params = list(alpha = 0.3),
              short.panel.labs = T) +
    facet_grid(.~AB_amplicon_pair) +
    stat_compare_means(label = "p.signif",method = "wilcox.test", comparisons=my_comparisons) +
    theme_base1 +
    scale_y_continuous(
      limits = c(0,1.3),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0.00","0.25","0.50","0.75","1.00")) + 
    ylab('mutation VAFs') +
    ggtitle(sprintf('%s:Using samples having amplicons with overlapping clustered mutations in tumor 1 2',ifelse(all_cutoff_passed_tumor12,'all_cutoff_passed_tumor12','unfiltered')))
  
  p2 <-
    ggboxplot(dt,
              x = "cl", y = "snv_vaf",
              color = 'AB_amplicon_pair',
              palette = 'npg',
              add.params = list(alpha = 0.3),
              short.panel.labs = T) +
    facet_grid(.~AB_amplicon_pair) +
    stat_compare_means(label = "p.format",method = "wilcox.test", comparisons=my_comparisons) +
    theme_base1 +
    scale_y_continuous(
      limits = c(0,1.3),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0.00","0.25","0.50","0.75","1.00")) + 
    ylab('mutation VAFs') +
    ggtitle(sprintf('%s:Using samples having amplicons with overlapping clustered mutations in tumor 1 2',ifelse(all_cutoff_passed_tumor12,'all_cutoff_passed_tumor12','unfiltered')))
  
  ggsave(p,
         file=sprintf('%s-boxplot-1.pdf', output_basename),
         width=10, height=7)
  ggsave(p2,
         file=sprintf('%s-boxplot-2.pdf', output_basename),
         width=10, height=7)
  write_tsv(dt_summary_count,
            sprintf('%s-d_summary_count.tsv', output_basename))
  write_tsv(dt_summary_count2,
            sprintf('%s-d_summary_count2.tsv', output_basename))
  saveRDS(dt,
          file=sprintf('%s-data.rds', output_basename))
  
  
  #' plot density plot #########
  ggdensity(dt%>%filter(is.element(cl, c('clustered_mutations-Shared_amplicon','non_clustered_mutations-Shared_amplicon'))),
            x = "snv_vaf",
            fill = "cl",
            color = "cl",
            palette = 'npg',
            add.params = list(alpha = 0.3),
            short.panel.labs = T) +
    facet_grid(.~AB_amplicon_pair)    
  
  
}
.Figure_5A(all_cutoff_passed_tumor12=TRUE, sel_study='na')
