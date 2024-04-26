.Figure_4AB <-function(all_cutoff_passed_tumor12=TRUE, sel_study='na') {
  
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
                              sprintf('Figure_4AB-%s', 'all_cutoff_passed_tumor12'))
  
  
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
  
  
  d_amp1                         <-d_amp1 %>%
    dplyr::select(-one_of('Intervals')) %>%
    left_join(dt_panecdnams_amplicons %>% dplyr::select(amplicon_barcode, Intervals), by='amplicon_barcode')
  
  d_amp2                         <-d_amp2 %>%
    dplyr::select(-one_of('Intervals')) %>%
    left_join(dt_panecdnams_amplicons %>% dplyr::select(amplicon_barcode, Intervals), by='amplicon_barcode')
  
  
  #' Create d_kataegis_calls from panecdnams_kataegis_calls #######
  #' And then add is_shared_in_kataegis to d_amp_kataegis_vafs.
  #' is_shared_in_kataegis= shared mutations present in kataegis calls
  d_kataegis_calls    <-panecdnams_kataegis_calls %>% 
    filter(is.element(aa_barcode, sel_samples$aa_barcode)) %>% 
    filter(is.element(patient_barcode, sel_patients$patient_barcode)) %>%
    dplyr::mutate(mut_tag=paste0(chr, ':',start,'_',REF,'/',ALT)) %>%
    mutate(patient_mut_tag=paste0(patient_barcode,'-',mut_tag))
  x_shared <-d_kataegis_calls %>%
    group_by(patient_mut_tag) %>%
    summarise(n=n()) %>%
    ungroup()
  x_shared <-x_shared %>% filter(n==2)
  d_kataegis_calls <-d_kataegis_calls %>%
    mutate(is_shared=ifelse(is.element(patient_mut_tag, x_shared$patient_mut_tag),'Shared','Private')) %>%
    dplyr::select(-patient_mut_tag)
  
  #' Check how many shared mutations each kataegis_group_barcode has
  #' Filter kataegis_group_barcode with n_is_shared_muttag>=2
  d_kataegis_group_barcode_shared <-d_kataegis_calls %>%
    group_by(aa_barcode, patient_barcode, tumor12, kataegis_group_barcode) %>%
    summarise(n_is_shared_muttag=sum(is_shared=='Shared')) %>%
    ungroup() %>%
    filter(n_is_shared_muttag>=2)
  
  #' Revise is_shared using d_kataegis_group_barcode_shared
  d_kataegis_calls <-d_kataegis_calls %>%
    group_by(kataegis_barcode) %>%
    mutate(is_shared_rev=is_shared,
           is_shared_rev=ifelse(is.element(kataegis_group_barcode, d_kataegis_group_barcode_shared$kataegis_group_barcode),is_shared_rev,'Private')) %>%
    ungroup() %>%
    mutate(is_shared_bak=is_shared,
           is_shared=is_shared_rev)
  
  d_amp1_kataegis_vafs <-d_amp1_kataegis_vafs %>%
    dplyr::select(-one_of(c('is_shared_in_kataegis'))) %>%
    left_join(d_kataegis_calls %>% dplyr::select(kataegis_tag, is_shared) %>% dplyr::rename(is_shared_in_kataegis=is_shared), by='kataegis_tag')
  d_amp2_kataegis_vafs <-d_amp2_kataegis_vafs %>%
    dplyr::select(-one_of(c('is_shared_in_kataegis'))) %>%
    left_join(d_kataegis_calls %>% dplyr::select(kataegis_tag, is_shared) %>% dplyr::rename(is_shared_in_kataegis=is_shared), by='kataegis_tag')
  
  #' =============================
  #' 
  #' plot Figure 4A
  #' 
  #' =============================
  .plot_4A <-function() {
    
    d_amp1                         <-d_amp1 %>% 
      dplyr::select(-one_of('Intervals')) %>%
      left_join(dt_panecdnams_amplicons %>% dplyr::select(amplicon_barcode, Intervals), by='amplicon_barcode') %>%
      mutate(amplicon1_has_kataegis=factor(amplicon1_has_kataegis, levels=c('No','Yes')))
    
    d_lab <-d_amp1 %>%
      group_by(amplicon_classification, amplicon1_has_similar_amplicon2, amplicon1_has_kataegis) %>%
      summarise(n=n()) %>%
      ungroup() %>%
      left_join(d_amp1 %>%
                  group_by(amplicon_classification, amplicon1_has_similar_amplicon2) %>%
                  summarise(N=n()) %>%
                  ungroup(),
                by=c('amplicon_classification', 'amplicon1_has_similar_amplicon2')) %>%
      dplyr::mutate(f=n/N) %>%
      dplyr::filter(amplicon1_has_kataegis=='Yes') %>%
      mutate(lab=paste0(sprintf('%.f',f*100),'%\n(',n,'/',N,')'),
             y=f,
             y=ifelse(y<0.1, 0.15, y))
    #' fisher.test
    pvals <-ddply(d_amp1, 'amplicon_classification', function(xd) {
      
      xx <-table(xd[,c('amplicon1_has_similar_amplicon2','amplicon1_has_kataegis')], useNA="ifany")
      p_chisq  <-chisq.test(as.matrix(xx)) # p_chisq$p.value
      p_fisher <-fisher.test(as.matrix(xx))
      
      p_null  <-nrow(filter(xd, amplicon1_has_similar_amplicon2=='No', amplicon1_has_kataegis=='Yes'))/nrow(filter(xd, amplicon1_has_similar_amplicon2=='No'))
      x_nb    <-nrow(filter(xd, amplicon1_has_similar_amplicon2=='Yes', amplicon1_has_kataegis=='Yes'))
      n_nb    <-nrow(filter(xd, amplicon1_has_similar_amplicon2=='Yes'))
      p_binom <-binom.test(x_nb, n_nb, p_null)
      
      y<-data.frame(amplicon='amplicon1',
                    p_chisq=p_chisq$p.value,
                    p_fisher=p_fisher$p.value,
                    odds_ratio=p_fisher$estimate['odds ratio'],
                    conf_int=paste0(p_fisher$conf.int, collapse = ','),
                    p_binom=p_binom$p.value,
                    binom_x=x_nb,
                    binom_n=n_nb,
                    binom_p_null=p_null)
      xx_df<-xd%>%group_by(amplicon1_has_kataegis,amplicon1_has_similar_amplicon2) %>% summarise(n=n()) %>% ungroup()
      cbind(y,as.data.frame(xx_df))
    })
    write_tsv(pvals, 
              sprintf('%s-plot_4A-pvalue.tsv', output_basename))
    p <-
      ggplot() +
      geom_bar(data=d_amp1, aes(x=amplicon1_has_similar_amplicon2, fill=amplicon1_has_kataegis), colour='black',position = position_fill()) +
      facet_grid(.~amplicon_classification) +
      scale_fill_manual(breaks=c('Yes','No'), values=c('#8073AC','#D8DAEB'), labels=c('Yes','No')) +
      theme_base1 +
      scale_y_continuous(labels = scales::percent) +
      ylab('% of amplicons in tumor 1') +
      guides(fill = guide_legend(title = "Amplicons in tumor 1 have\n overlapping clustered mutations")) +
      geom_text(d=d_lab, aes(x=amplicon1_has_similar_amplicon2, y=y, label=lab), angle=90) +
      xlab('Amplicons in tumor 1 have\nsimilar amplicons in tumor 2') +
      ggtitle(sprintf('%s\nsamples with amplicons',ifelse(all_cutoff_passed_tumor12,'all_cutoff_passed_tumor12','unfiltered')))
    ggsave(p,
           file=sprintf('%s-plot_4A-barplot.pdf', output_basename),
           width=7, height=5)
    
  }
  .plot_4A()
  
  #' =============================
  #' 
  #' plot Figure 4B
  #' 
  #' =============================
  .plot_4B <-function() {
    d_amp2 <-d_amp2 %>%
      mutate(amplicon2_has_kataegis=factor(amplicon2_has_kataegis, levels=c('No','Yes')))
    
    #' overlapping kataegis in shared vs private ecDNAs, shared vs private non-ecDNAs in GLASS+HMF
    d_lab <-d_amp2 %>%
      group_by(amplicon_classification, amplicon2_has_similar_amplicon1, amplicon2_has_kataegis) %>%
      summarise(n=n()) %>%
      ungroup() %>%
      left_join(d_amp2 %>%
                  group_by(amplicon_classification, amplicon2_has_similar_amplicon1) %>%
                  summarise(N=n()) %>%
                  ungroup(),
                by=c('amplicon_classification', 'amplicon2_has_similar_amplicon1')) %>%
      dplyr::mutate(f=n/N) %>%
      dplyr::filter(amplicon2_has_kataegis=='Yes') %>%
      mutate(lab=paste0(sprintf('%.f',f*100),'%\n(',n,'/',N,')'),
             y=f,
             y=ifelse(y<0.1, 0.15, y))
    #' fisher.test
    pvals <-ddply(d_amp2, 'amplicon_classification', function(xd) {
      xx <-table(xd[,c('amplicon2_has_similar_amplicon1','amplicon2_has_kataegis')], useNA="ifany")
      p_chisq <-chisq.test(as.matrix(xx)) # p_chisq$p.value
      p_fisher <-fisher.test(as.matrix(xx))
      
      p_null  <-nrow(filter(xd, amplicon2_has_similar_amplicon1=='No', amplicon2_has_kataegis=='Yes'))/nrow(filter(xd, amplicon2_has_similar_amplicon1=='No'))
      x_nb    <-nrow(filter(xd, amplicon2_has_similar_amplicon1=='Yes', amplicon2_has_kataegis=='Yes'))
      n_nb    <-nrow(filter(xd, amplicon2_has_similar_amplicon1=='Yes'))
      p_binom <-binom.test(x_nb, n_nb, p_null)
      
      y<-data.frame(amplicon='amplicon2',
                    p_chisq=p_chisq$p.value,
                    p_fisher=p_fisher$p.value,
                    odds_ratio=p_fisher$estimate['odds ratio'],
                    conf_int=paste0(p_fisher$conf.int, collapse = ','),
                    p_binom=p_binom$p.value,
                    binom_x=x_nb,
                    binom_n=n_nb,
                    binom_p_null=p_null)
      
      
      xx_df<-xd%>%group_by(amplicon2_has_kataegis,amplicon2_has_similar_amplicon1) %>% summarise(n=n()) %>% ungroup()
      cbind(y,as.data.frame(xx_df))
    })
    write_tsv(pvals, 
              sprintf('%s-plot_4B-pvalue.tsv', output_basename))
    p<-
      ggplot() +
      geom_bar(data=d_amp2, aes(x=amplicon2_has_similar_amplicon1, fill=amplicon2_has_kataegis), colour='black',position = position_fill()) +
      facet_grid(.~amplicon_classification) +
      scale_fill_manual(breaks=c('Yes','No'), values=c('#E08214','#FEE0B6'), labels=c('Yes','No')) +
      theme_base1 +
      scale_y_continuous(labels = scales::percent) +
      ylab('% of amplicons in tumor 2') +
      guides(fill = guide_legend(title = "Amplicons in tumor 2 have\n overlapping clustered mutations")) +
      geom_text(d=d_lab, aes(x=amplicon2_has_similar_amplicon1, y=y, label=lab), angle=90) +
      xlab('Amplicons in tumor 2 have\nsimilar amplicons in tumor 1') +
      ggtitle(sprintf('%s\nsamples with amplicons',ifelse(all_cutoff_passed_tumor12,'all_cutoff_passed_tumor12','unfiltered')))
    ggsave(p,
           file=sprintf('%s-plot_4B-barplot.pdf', output_basename),
           width=7, height=5)
    
    
  }
  .plot_4B()
  
  
}
.Figure_4AB(all_cutoff_passed_tumor12=TRUE, sel_study='na')

