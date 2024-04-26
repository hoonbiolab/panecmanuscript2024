.Figure_2A <-function(matched_tissue=TRUE) {
  
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
  
  stopifnot(matched_tissue)
  
  output_basename <-file.path(output_basedir, 'nat_genet_github_1st',
                              sprintf('Figure_2A-%s', ifelse(matched_tissue,'matched_tissue','all_tissue')))
  
  require(survival)
  require(survminer)
  
  survdata <-panecdnass_survivals %>%
    filter(is.element(aa_barcode, filter(panecdnass_samples, tissue_matched_primary_advanced=='matched_tissue')$aa_barcode)) %>%
    dplyr::select(-one_of(c('source', 'study', 'primaryTumorLocation', 'sample_classification'))) %>%
    left_join(panecdnass_samples %>% filter(tissue_matched_primary_advanced=='matched_tissue') %>% dplyr::select(aa_barcode, source, study, primaryTumorLocation, sample_classification), by='aa_barcode') %>%
    mutate(primaryTumorLocation=ifelse(is.na(primaryTumorLocation),'null',primaryTumorLocation)) %>%
    group_by(patient_barcode) %>%
    mutate(survival_days=ifelse(vital_status=='DEAD',days_to_death, days_to_last_followup),
           death=ifelse(vital_status=='DEAD',1,0)) %>%
    ungroup() %>%
    filter(!is.na(survival_days))
  
  ll <-.factor_sample_classification(survdata$sample_classification, method = 'ecDNA')
  survdata$sample_classification <-ll$cl
  survdata <-survdata %>%
    mutate(source=factor(source, levels=c("Primary","Advanced"))) %>%
    mutate(sample_classification=factor(sample_classification, levels=c('ecDNA','ChrAmp','NoAmp')))
  
  fit <- survminer::surv_fit( Surv(survival_days, death) ~ sample_classification + source, data = survdata )
  z <-survminer::ggsurvplot(
    fit = fit,
    xlab = "Days",
    ylab = "Survival probability",
    legend.title="Classification",
    legend='top',
    #legend.labs = Xaa_class_orders$class,
    xlim=c(0,2000),
    ncensor.plot = FALSE,
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.y.text = FALSE,
    break.time.by = 500,
    font.title    = c(16, "bold", "black"),
    font.caption  = c(14, "plain", "black"),
    font.x        = c(14, "bold.italic", "black"),
    font.y        = c(14, "bold.italic", "black"),
    font.xtickslab = c(12, "plain", "black"),
    pval = FALSE,
    #palette=as.character(Xaa_class_orders$color),
    palette=c("#D7191C", '#CC79A7',"#6495ED", '#56B4E9',  '#999999',"gray"),
    data=survdata,
    title=sprintf('total=%d', nrow( survdata )  )
  )
  z$plot +theme_bw() + facet_grid(. ~ factor(source, levels=c('Primary','Advanced')))
  
  ggsave(z$plot +theme_bw() + facet_grid(. ~ factor(source, levels=c('Primary','Advanced'))),
         file=sprintf('%s-kmplot.pdf', output_basename),
         width=12, height=5)
  
  surv_diff_primary  <- survdiff(Surv(survival_days, death) ~ sample_classification, data = survdata %>% filter(source=='Primary'))
  surv_diff_advanced <- survdiff(Surv(survival_days, death) ~ sample_classification, data = survdata %>% filter(source=='Advanced'))
  print(surv_diff_primary)
  print(surv_diff_advanced)
  
  saveRDS(survdata, file=sprintf('%s-data.rds',output_basename))
}

.Figure_2A(matched_tissue=TRUE)
