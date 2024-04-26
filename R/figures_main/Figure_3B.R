.Figure_3B <-function(all_cutoff_passed_tumor12=TRUE) {
  
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
  
  output_basename <-file.path(output_basedir, 'nat_genet_github_1st',
                              sprintf('Figure_3B-%s', ifelse(all_cutoff_passed_tumor12,'all_cutoff_passed_tumor12','')))
  
  dt_samples <-panecdnams_samples %>%
    filter(is.element(aa_barcode, filter(panecdnams_platinum, all_cutoff_passed_tumor12=='passed')$aa_barcode))
  
  #' Get filtered amplicons
  dt_amplicons <-panecdnams_amplicons %>%
    filter(is.element(aa_barcode, dt_samples$aa_barcode)) %>%
    filter(amplicon_classification!="No amp/Invalid")
  
  #' Get significantly similar amplicon pairs ###
  #' Filter by platinum
  #' Exclude Invalid amps
  #' same amplicon classification in tumor 1 and 2
  #' pvalue<0.05
  dt_sig_amp_pairs <-panecdnams_amplicon_similarity %>%
    mutate(SimScorePvalue=as.numeric(SimScorePvalue)) %>%
    filter(SimScorePvalue<0.05) %>%
    filter(is.element(aa_barcode1, dt_samples$aa_barcode)) %>%
    filter(is.element(aa_barcode2, dt_samples$aa_barcode)) %>%
    filter(!is.element(amplicon_classification1, 'No amp/Invalid')) %>%
    filter(!is.element(amplicon_classification2, 'No amp/Invalid')) %>%
    filter(amplicon_classification1==amplicon_classification2)
  
  
  #' #########################################
  #' 
  #' Plot
  #' 
  #' #########################################
  d_amp1 <-dt_amplicons %>%
    filter(tumor12=='tumor1') %>%
    filter(amplicon_classification!="No amp/Invalid") %>%
    mutate(has_sim_amp=ifelse(is.element(amplicon_barcode, dt_sig_amp_pairs$Amp1),'Yes','No'))
  
  stopifnot(all(is.element(dt_sig_amp_pairs$Amp2, filter(dt_amplicons, tumor12=='tumor2')$amplicon_barcode)))
  d_amp2 <-dt_amplicons %>%
    filter(tumor12=='tumor2') %>%
    filter(amplicon_classification!="No amp/Invalid") %>%
    mutate(has_sim_amp=ifelse(is.element(amplicon_barcode, dt_sig_amp_pairs$Amp2),'Yes','No'))
  
  d_amp <-bind_rows(d_amp1, d_amp2) %>%
    mutate(tumor12=factor(tumor12, levels=c('tumor1','tumor2'))) %>%
    mutate(has_sim_amp=factor(has_sim_amp, levels=c('No','Yes')))
  
  stopifnot(identical(sort(unique(as.character(d_amp$amplicon_classification))), 
                      sort(c("ecDNA","ChrAmp"))))
  l <-.factor_amplicon_classification(d_amp$amplicon_classification, method = 'ecDNA')
  d_amp$amplicon_classification  <-l$cl
  aa_class_orders                 <-l$aa_class_orders
  
  #' summarize w.r.t. tumor 1 amplicons ########
  d_amp_summary <-d_amp %>%
    group_by(tumor12, amplicon_classification, has_sim_amp) %>%
    summarise(n=n()) %>%
    ungroup() %>%
    left_join(d_amp %>%
                group_by(tumor12, amplicon_classification) %>%
                summarise(N=n()) %>%
                ungroup(),
              by=c('tumor12','amplicon_classification')) %>%
    mutate(f=n/N) %>%
    mutate(lab=paste0(formatC(f*100, format = "f", digits = 1),'% (',n,' / ',N,')')) %>%
    mutate(y=f)
  
  #' chsq pval w.r.t. tumor 1 amplicons ########
  require(dplyr)
  require(knitr)
  require(pander)
  require(broom)
  d_amp <-d_amp %>% 
    mutate(amplicon_classification=factor(amplicon_classification, levels=rev(c("ecDNA","ChrAmp"))))
  
  d_amp_summary_chisq <-ddply(d_amp,c('tumor12'), function(xd) {
    p <-chisq.test(table(xd[,c('has_sim_amp','amplicon_classification')]))
    f1 <-glm(has_sim_amp~amplicon_classification, data=xd, family=binomial)
    or <- exp(coef(f1))
    ci <- exp(confint(f1))
    or_data <- data.frame(predictor = names(or), odds_ratio = or, lower_ci = ci[,1], upper_ci = ci[,2])
    xtidy<-tidy(f1)
    stopifnot(identical(sort(as.character(xtidy$term)),
                        sort(as.character(or_data$predictor))))
    stopifnot(all(!is.element(colnames(xtidy), colnames(or_data))))
    or_data <-or_data %>% left_join(xtidy%>% dplyr::rename(predictor=term), by='predictor') %>%
      mutate(is_sig=ifelse(p.value<0.05, 'sig','n.s.'),
             is_sig=factor(is_sig, levels=c('sig','n.s.'))) %>%
      filter(predictor=='amplicon_classificationecDNA')
    colnames(or_data) <-paste0('OR.', colnames(or_data))
    or_data <-or_data %>% mutate(pval_chisq=p$p.value)
    return(or_data)
  })
  
  d_amp <-d_amp %>% 
    mutate(amplicon_classification=factor(amplicon_classification, levels=c("ecDNA","ChrAmp")))
  
  p <-
    ggplot() +
    geom_bar(data=d_amp, aes(x=amplicon_classification, fill=has_sim_amp), position = position_fill()) +
    geom_text(data=d_amp_summary %>% filter(has_sim_amp=='Yes'), aes(x=amplicon_classification, y=y+0.05, label=lab), angle=90, colour='#0072B2', size=4.5, hjust=0) +
    scale_y_continuous(labels = scales::percent) +
    #scale_fill_manual(breaks=aa_class_orders$class, labels=aa_class_orders$class, values=aa_class_orders$color) +
    theme_base1 +
    facet_grid(.~tumor12) +
    scale_fill_manual(breaks = c('Yes','No'), values = c('#0072B2','#E1E7EA'), labels = c('Yes','No')) +
    xlab('Amplicon classification') +
    ylab('% of amplicons') +
    ggtitle(sprintf('%s, %s',paste(sort(unique(d_amp$study)), collapse="/"), ifelse(all_cutoff_passed_tumor12,'all_cutoff_passed_tumor12',''))) +
    labs(fill='T1 and T2 share similar amplicon')
  print(p)
  
  saveRDS(list(d_amp=d_amp,
               d_amp_summary=d_amp_summary,
               d_amp_summary_chisq=d_amp_summary_chisq), 
          sprintf('%s-data.rds', output_basename))
  
}

.Figure_3B(all_cutoff_passed_tumor12=TRUE)
