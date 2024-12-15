# - - - - STATS Z-SCORES - - - - #
library(compareGroups)
library(ggstatsplot)
library(dplyr)
library(ggplot2)
library(tidyr)

#-----------------------------------
# - - - - - Prepare data - - - - - -
#-----------------------------------
# Set paths for data
# Z-scores. Rows: subjects. Columns: brain regions Z-scores
z.score <- ""
# ID, age, sex, scanner, sample, dcode, subID, timepoint
covariates.test <- ""
# Z-scores. Rows: subjects. Columns: brain regions features
regional.features <-""

# Read data
df_Z_scores <- read.csv(z.score)
df_cov_te <- read.csv(covariates.test)#[,-12]
df_features <- read.csv(regional.features)
df_complete <- cbind(df_cov_te,df_Z_scores) # covariates and Z-scores

clinical_samples <- unique(df_complete[df_complete$dcode==1,"sample"]) # get clinical samples names
# add all clinical sites at the end of the list to evaluate them together
clinical_samples <- c(clinical_samples,list(clinical_samples))

# create vector with region names
features <- setdiff(names(df_Z_scores), names(df_cov_te))

# Brief check of numbers within each group
table(df_complete$dcode) #diagnosis
table(df_complete$sex) #sex

#-----------------------------------
#-------- Only males ---------------
#-----------------------------------
df_Z_scores <- df_Z_scores %>%  filter(df_cov_te$sex == 1)
df_features <- df_features %>%  filter(df_cov_te$sex == 1)
df_complete <- df_complete %>%  filter(df_cov_te$sex == 1)
df_cov_te <- df_cov_te %>%  filter(df_cov_te$sex == 1)
#-----------------------------------
#-------- Only females ---------------
#-----------------------------------
df_Z_scores <- df_Z_scores %>%  filter(df_cov_te$sex == 2)
df_features <- df_features %>%  filter(df_cov_te$sex == 2)
df_complete <- df_complete %>%  filter(df_cov_te$sex == 2)
df_cov_te <- df_cov_te %>%  filter(df_cov_te$sex == 2)

#-----------------------------------
# - - Case-control global feature -  
#-----------------------------------
for (sample in clinical_samples) {
  df_feature <- df_features
  df_feature<- cbind(df_cov_te,df_feature)
  df_feature <- df_feature[df_feature$sample %in% sample,]
  tps <- unique(df_feature$timepoint)
  
  for (tp in tps) {
    # for all clinical samples, only look at baseline subjects
    # if you want to run it only on first timepoint subjects of all clinical samples uncomment the following line
    # if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    #filter df per timepoint
    df_feature_tp <- df_feature[df_feature$timepoint == tp, ]
    dcode <- df_feature_tp$dcode
    # create vector with region names
    features <- setdiff(names(df_feature_tp), names(df_cov_te))
    
    # prepare data for t-test
    global_feature <- rowMeans(df_feature_tp[, features])
    a <- as.data.frame(cbind(dcode, global_feature))
    a$dcode <- ifelse(a$dcode == 0, 1, ifelse(a$dcode == 1, 0, a$dcode))
    a$dcode <- factor(a$dcode, labels = c('Schizophrenia', 'Healthy'))
    #dcode <- a$dcode
    # t-test (automatically chooses between normal or non-normal distribution & parametric or non-parametric test)
    t_test_global <- compareGroups(dcode ~ ., method = NA, data = a, alpha = 0.05)
    print(paste0("P-value global feature (", sample_name,") timepoint ",tp,": ",t_test_global[["global_feature"]][["p.overall"]]))
    normality <- attributes(t_test_global[[1]])$method[2]
    if (normality=="normal"){param <- "p" #parametric test
    } else {#non-normal: non-parametric test
      param <- "np"}

    plot <- ggbetweenstats(data=a, x=dcode, y=global_feature,
                           ylab = "Global MS",
                           xlab="",
                           type=param)+
                           #caption=paste0(" ")+
                           scale_color_manual(values = c("#1F78B4","#33A02C"))
    sub_1 <- plot[["plot_env"]][["subtitle"]][[2]][3]
    sub_1 <- gsub("\\(|\\)", "", sub_1)
    sub_2 <- plot[["labels"]][["subtitle"]][[4]][3]
    sub_2 <- gsub("\\(|\\)", "", sub_2)
    sub_3 <- plot[["labels"]][["subtitle"]][[5]][3]
    sub_3 <- gsub("\\(|\\)", "", sub_3)
    sub_3 <- as.numeric(gsub("[^0-9.-]", "", sub_3))
    sub_4 <- plot[["labels"]][["subtitle"]][[6]][2]
    sub_4 <- gsub("\\(|\\)", "", sub_4)
    sub_5 <- plot[["plot_env"]][["subtitle"]][[2]][2]
    sub_5 <- as.numeric(gsub("[^0-9.]", "", sub_5))
    if (normality=="normal"){plot <- plot + 
      labs(subtitle = ggplot2::expr(paste(italic("t")["Welch"],"(",!!sub_5,")"," = ", !!sub_1 ,
                                          ", ", widehat(italic("d"))["Cohen"], " = ", !!sub_2 ,
                                          ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]")))
    } else{plot <- plot + 
      labs(subtitle = ggplot2::expr(paste(italic("W")["Mann-Whitney"]," = ", !!trunc(as.numeric(sub_1)) ,
                                          ", ", widehat(italic("r"))["biserial"]^"rank", " = ", !!sub_2 ,
                                          ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]"
      )))}
    plot <- plot + theme(axis.text = element_text(size = 15),  
                         axis.title = element_text(size = 16),  
                         plot.subtitle = element_text(size = 16, face="italic"),
                         plot.title = element_text(size = 18))
    if (length(tps) > 1) {
      png(file = paste0("global_feature_", sample_name, "_tp", tp, ".png"), width = 6.5, height = 4.5, units = "in", res=600)
    } else {
      png(file = paste0("global_feature_", sample_name, ".png"), width = 6.5, height = 4.5, units = "in", res=600)
    }
    print(plot)
    dev.off()
  }
}
rm(a,df_feature_tp,df_feature,plot,t_test_global,dcode,normality,param,sample,sample_name,tp,tps)

#-----------------------------------
# - - Case-control global Z-score  -
#-----------------------------------
for (sample in clinical_samples) {
  df <- df_Z_scores
  df<- cbind(df,df_cov_te)
  df <- df[df$sample %in% sample,]
  tps <- unique(df$timepoint)
  for (tp in tps) {
    # for all clinical samples, only look at baseline subjects
    # if you want to run it only on first timepoint subjects of all clinical samples uncomment the following line
    #if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    df_tp <- df[df$timepoint == tp, ]
    dcode <- df_tp$dcode
    
    global <- rowMeans(df_tp[,features])
    a <- as.data.frame(cbind(dcode,global))
    t_test_global <- compareGroups(dcode ~ .,method = NA, data = a,alpha = 0.05)
    print(paste0("P-value global Z-score (", sample_name,") timepoint ",tp,": ",t_test_global[["global"]][["p.overall"]]))
    normality <- attributes(t_test_global[[1]])$method[2]
    if (normality=="normal"){param <- "p" #parametric test
    } else {#non-normal: non-parametric test
      param <- "np"}
    a$dcode <- factor(a$dcode, labels = c('Controls', 'Patients'))
    plot <- ggbetweenstats(data=a, x=dcode, y=global,
                           ylab = "Global Z-score",
                           xlab="Group",
                           type=param, 
                           caption=paste0("Sample: ",sample_name,", timepoint ",tp))
    if (length(tps) > 1) {
      png(file = paste0("global_Z_score_",sample_name,"_tp",tp,".png"), width=6, height=4.5, units="in", res=600)
    } else {
      png(file = paste0("global_Z_score_",sample_name,".png"), width=6, height=4.5, units="in", res=600)
    }
    print(plot)
    dev.off()
  }
}
rm(a,df_tp,df,plot,t_test_global,dcode,global,global_feature,normality,param,sample,sample_name,tp,tps)

#-----------------------------------
# - - Binarised outlier matrix  - - 
#-----------------------------------

bin_outlier_function <- function(df, cov_te) {
  outlier_threshold_u <- 1.96 # 95% confidence interval
  outlier_threshold_l <- -1.96 
  
  # Feature names
  features <- setdiff(names(df), names(cov_te))
  
  # Only POS deviations
  df_pos <- as.data.frame(ifelse(df[, features] > outlier_threshold_u, 1, 0))
  # Only NEG deviations
  df_neg <- as.data.frame(ifelse(df[, features] < outlier_threshold_l, 1, 0))
  # All deviations
  df2 <- as.data.frame(ifelse(df[, features] > outlier_threshold_u | df[, features] < outlier_threshold_l, 1, 0))
  df2 <- df2 %>% rename_all(paste0, "_bin") # Rename all binarised columns to have the suffix "_bin"
  
  # Calculate total outlier score: total, positive and negative
  df$total_outlier_score <- rowSums(df2) # Sum all of binarised outliers
  df$total_outlier_score_POS <- rowSums(df_pos)
  df$total_outlier_score_NEG <- rowSums(df_neg)
  
  # Add outlier matrix to original df
  df <- cbind(df, df2)
  rm(outlier_threshold_u, outlier_threshold_l, df2, df_pos, df_neg)
  return(df)
}

#-----------------------------------
# - - Outlier subjects and regions - 
#-----------------------------------
library(RColorBrewer)

for (sample in clinical_samples) {
  df <- df_Z_scores
  cov_te <- df_cov_te
  df<- cbind(df_cov_te,df)
  df <- df[df$sample %in% sample,]
  tps <- unique(df$timepoint)
  for (tp in tps) {
    print(tp)
    # for all clinical samples, only look at baseline subjects
    # if you want to run it only on frist timepoint subjects of all clinical samples uncomment the following line
    #if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    df_tp <- df[df$timepoint == tp, ]
    df_tp <- bin_outlier_function(df_tp, cov_te)
    if (tp==1){tp_text="Baseline"}else{tp_text="Follow-up"}
    # - - % OUTLIER SUBJECTS PER REGION  - -
    
    # count outliers HC
    out_controls <- df_tp %>% filter(dcode == "0") %>%
      select(contains("_bin"))
    out_controls <- colSums(out_controls)/nrow(out_controls)*100
    # count outliers HC positive deviations only
    df3 <- as.data.frame(ifelse(df_tp[df_tp$dcode==0,features] > 1.96, 1, 0))
    perc_out_pos_con <- colSums(df3)/nrow(df3)*100
    # count outliers HC negative deviations only
    df3 <- as.data.frame(ifelse(df_tp[df_tp$dcode==0,features] < -1.96, 1, 0))
    perc_out_neg_cont <- colSums(df3)/nrow(df3)*100
    out_controls_df <- as.data.frame(cbind(features,out_controls,perc_out_neg_cont,perc_out_pos_con))
    colnames(out_controls_df) <- c("features","perc_outliers","perc_outliers_neg","perc_outliers_pos")
    write.csv(out_controls_df, file = paste0("",sample_name,"_tp",tp,".csv"),row.names=FALSE)
    
    # count outliers PATS
    out_pats <- df_tp %>% filter(dcode == "1") %>%
      select(contains("_bin"))
    out_pats <- colSums(out_pats)/nrow(out_pats)*100
    # count outliers PATS negative deviations only
    df3 <- as.data.frame(ifelse(df_tp[df_tp$dcode==1,features] < -1.96, 1, 0))
    perc_out_neg_pat <- colSums(df3)/nrow(df3)*100
    # count outliers PATS positive deviations only
    df3 <- as.data.frame(ifelse(df_tp[df_tp$dcode==1,features] > 1.96, 1, 0))
    perc_out_pos_pat <- colSums(df3)/nrow(df3)*100
    
    out_pats_df <- as.data.frame(cbind(features,out_pats,perc_out_neg_pat,perc_out_pos_pat))
    colnames(out_pats_df) <- c("features","perc_outliers","perc_outliers_neg","perc_outliers_pos")
    write.csv(out_pats_df, file = paste0("",sample_name,"_tp",tp,".csv"),row.names=FALSE)
    
    rm(df3,perc_out_neg_cont,perc_out_pos_con,perc_out_pos_pat,perc_out_neg_pat)
    
    # plots
    patients_df <- df_tp[df_tp$dcode== 1, ]
    controls_df <- df_tp[df_tp$dcode== 0, ]
    
    out_perc <- data.frame(region = 0:length(features), perc_pos_pat=NA, perc_neg_pat=NA, perc_pos_hc=NA, perc_neg_hc=NA)
    
    for (i in 1:nrow(out_perc)) {
      #Positive outliers patients
      count <- sum(patients_df$total_outlier_score_POS == out_perc$region[i])
      percentage_out <- count / nrow(patients_df) * 100
      out_perc$perc_pos_pat[i] <- percentage_out
      # Negative outliers patients
      count <- sum(patients_df$total_outlier_score_NEG == out_perc$region[i])
      percentage_out <- count / nrow(patients_df) * 100
      out_perc$perc_neg_pat[i] <- percentage_out
      # Positive outliers controls
      count <- sum(controls_df$total_outlier_score_POS == out_perc$region[i])
      percentage_out <- count / nrow(controls_df) * 100
      out_perc$perc_pos_hc[i] <- percentage_out
      # Negative outliers controls
      count <- sum(controls_df$total_outlier_score_NEG == out_perc$region[i])
      percentage_out <- count / nrow(controls_df) * 100
      out_perc$perc_neg_hc[i] <- percentage_out}
    
    out_perc_pos <- out_perc %>%
      gather(key = "type_dev", value = "perc_out", -region) %>%  
      filter(type_dev == "perc_pos_hc" | type_dev == "perc_pos_pat") %>%  
      filter(perc_out != 0)
    
    plot<-ggplot(out_perc_pos, aes(x = region, y = perc_out, fill = type_dev)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = "Extremely positive deviated regions", y = "% subjects", 
           title = "Supranormal deviations") +
      scale_fill_manual(values = c("#33A02C", "#1F78B4"),
      #scale_fill_manual(values = brewer.pal(3, "Dark2"),
                        labels = c('Healthy', 'Schizophrenia'), name = NULL) +
      theme_minimal()
    plot <- plot + theme(text = element_text(size = 14),
                         axis.text = element_text(size = 17),  
                         axis.title = element_text(size = 17),  
                         plot.subtitle = element_text(size = 11),
                         plot.title = element_text(size = 18)) 
    png(file=paste0("",sample_name,"_tp",tp,".png"),width=6.5, height=4.5, units="in",res=600)
    print(plot)
    dev.off()
    
    out_perc_neg <- out_perc %>%
      gather(key = "type_dev", value = "perc_out", -region) %>%  
      filter(type_dev == "perc_neg_hc" | type_dev == "perc_neg_pat") %>%  
      filter(perc_out != 0)
    plot<-ggplot(out_perc_neg, aes(x = region, y = perc_out, fill = type_dev)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = "Extremely negative deviated regions", y = "% subjects", 
           title = "Infranormal deviations") +
      scale_fill_manual(values = c("#33A02C", "#1F78B4"),
      #scale_fill_manual(values = brewer.pal(3, "Dark2"),
                        labels = c('Healthy', 'Schizophrenia'), name = NULL) +
      theme_minimal()
    plot <- plot + theme(text = element_text(size = 14),
                         axis.text = element_text(size = 17),  
                         axis.title = element_text(size = 17),  
                         plot.subtitle = element_text(size = 11),
                         plot.title = element_text(size = 18))  
    png(file=paste0("",sample_name,"_tp",tp,".png"),width=6.5, height=4.5, units="in",res=600)
    print(plot)
    dev.off()
    
    rm(controls_df,patients_df,out_perc,out_perc_neg,out_perc_pos)
    
    # - - OUTLIER REGIONS PER SUBJECT  - - -
    
    # Scatter plot to see outliers according to diagnosis 
    df_tp$dcode <- as.factor(df_tp$dcode)
    plot(df_tp$total_outlier_score, col=df_tp$dcode, pch=as.numeric(df_tp$dcode))
    
    # t-test total_outlier_score
    dcode <- factor(df_tp$dcode)
    total_outlier_score <- df_tp$total_outlier_score
    a <- as.data.frame(cbind(dcode,total_outlier_score))
    
    t_test_out <- compareGroups(dcode ~ .,method = NA, data = a,alpha = 0.05)
    cat(paste0("P-value total outlier score (site=", sample_name," tp ",tp,"): ", t_test_out[["total_outlier_score"]][["p.overall"]],"\n"))
    normality <- attributes(t_test_out[[1]])$method[2]
    if (normality=="normal"){param <- "p" #parametric test
    } else {#non-normal: non-parametric test
      param <- "np"}
    a$dcode <- factor(a$dcode, labels = c('Healthy', 'Schizophrenia'))
    #custom_palette <- c("blue" = brewer.pal(3, "Blues")[2], "green" = brewer.pal(3, "Greens")[2])
    plot <- ggbetweenstats(data=a, x=dcode, y=total_outlier_score,
                           ylab = "Total outlier regions",
                           xlab="",
                           type=param, 
                           caption=paste0(tp_text),
                           ggplot.component = list(ggplot2::scale_y_continuous(
                             breaks = seq(0,40,10),
                             limits = (c(0,40)))))+
                           scale_color_manual(values = c("#33A02C", "#1F78B4")) # brewer paired palette
    sub_1 <- plot[["plot_env"]][["subtitle"]][[2]][3]
    sub_1 <- gsub("\\(|\\)", "", sub_1)
    sub_2 <- plot[["labels"]][["subtitle"]][[4]][3]
    sub_2 <- gsub("\\(|\\)", "", sub_2)
    sub_3 <- plot[["labels"]][["subtitle"]][[5]][3]
    sub_3 <- gsub("\\(|\\)", "", sub_3)
    sub_3 <- as.numeric(gsub("[^0-9.-]", "", sub_3))
    sub_4 <- plot[["labels"]][["subtitle"]][[6]][2]
    sub_4 <- gsub("\\(|\\)", "", sub_4)
    sub_5 <- plot[["plot_env"]][["subtitle"]][[2]][2]
    sub_5 <- as.numeric(gsub("[^0-9.]", "", sub_5))
    if (normality=="normal"){plot <- plot + 
      labs(subtitle = ggplot2::expr(paste(italic("t")["Welch"],"(",!!sub_5,")"," = ", !!sub_1 ,
                                          ", ", widehat(italic("d"))["Cohen"], " = ", !!sub_2 ,
                                          ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]")))
    } else{plot <- plot + 
      labs(subtitle = ggplot2::expr(paste(italic("W")["Mann-Whitney"]," = ", !!trunc(as.numeric(sub_1)) ,
                                          ", ", widehat(italic("r"))["biserial"]^"rank", " = ", !!sub_2 ,
                                          ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]"
      )))}
    plot <- plot + theme(axis.text = element_text(size = 15),  
                         axis.title = element_text(size = 16),  
                         plot.subtitle = element_text(size = 15, face="italic"),
                         plot.title = element_text(size = 18))
    #plot <- plot + theme(text = element_text(size = 14),
    #                     axis.text = element_text(size = 17),  
    #                     axis.title = element_text(size = 17),  
    #                     plot.subtitle = element_text(size = 11),
    #                     plot.title = element_text(size = 18)) 
    png(file=paste0("",sample_name,"_tp",tp,".png"),width=6, height=4.5, units="in",res=600)
    print(plot)
    dev.off()
    
    # POSITIVE AND NEGATIVE DEVIATIONS
    
    # t-test total_outlier_score POS
    #dcode <- factor(df_tp$dcode)
    total_outlier_score <- df_tp$total_outlier_score_POS
    a <- as.data.frame(cbind(dcode,total_outlier_score))
    t_test_out <- compareGroups(dcode ~ .,method = NA, data = a,alpha = 0.05)
    cat(paste0("P-value total outlier score POS deviations (site=", sample_name," tp ",tp,"): ", t_test_out[["total_outlier_score"]][["p.overall"]],"\n"))
    normality <- attributes(t_test_out[[1]])$method[2]
    if (normality=="normal"){param <- "p" #parametric test
    } else {#non-normal: non-parametric test
      param <- "np"}
    a$dcode <- factor(a$dcode, labels = c('Healthy', 'Schizophrenia'))
    plot <- ggbetweenstats(data=a, x=dcode, y=total_outlier_score,
                           ylab = "Total positive outlier regions",
                           xlab="",
                           type=param, 
                           caption=paste0(tp_text))+
      scale_color_manual(values = c("#33A02C", "#1F78B4"))
    sub_1 <- plot[["plot_env"]][["subtitle"]][[2]][3]
    sub_1 <- gsub("\\(|\\)", "", sub_1)
    sub_2 <- plot[["labels"]][["subtitle"]][[4]][3]
    sub_2 <- gsub("\\(|\\)", "", sub_2)
    sub_3 <- plot[["labels"]][["subtitle"]][[5]][3]
    sub_3 <- gsub("\\(|\\)", "", sub_3)
    sub_3 <- as.numeric(gsub("[^0-9.-]", "", sub_3))
    sub_4 <- plot[["labels"]][["subtitle"]][[6]][2]
    sub_4 <- gsub("\\(|\\)", "", sub_4)
    sub_5 <- plot[["plot_env"]][["subtitle"]][[2]][2]
    sub_5 <- as.numeric(gsub("[^0-9.]", "", sub_5))
    if (normality=="normal"){plot <- plot + 
      labs(subtitle = ggplot2::expr(paste(italic("t")["Welch"],"(",!!sub_5,")"," = ", !!sub_1 ,
                                          ", ", widehat(italic("d"))["Cohen"], " = ", !!sub_2 ,
                                          ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]")))
    } else{plot <- plot + 
      labs(subtitle = ggplot2::expr(paste(italic("W")["Mann-Whitney"]," = ", !!trunc(as.numeric(sub_1)) ,
                                          ", ", widehat(italic("r"))["biserial"]^"rank", " = ", !!sub_2 ,
                                          ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]"
      )))}
    plot <- plot + theme(axis.text = element_text(size = 15),  
                         axis.title = element_text(size = 16),  
                         plot.subtitle = element_text(size = 15, face="italic"),
                         plot.title = element_text(size = 18))
    #plot <- plot + theme(text = element_text(size = 14),
    #                     axis.text = element_text(size = 17),  
    #                     axis.title = element_text(size = 17),  
    #                     plot.subtitle = element_text(size = 11),
    #                     plot.title = element_text(size = 18)) 
    png(file=paste0("",sample_name,"_tp",tp,".png"),width=6, height=4.5, units="in",res=600)
    print(plot)
    dev.off()
    
    # t-test total_outlier_score NEG
    #dcode <- factor(df_tp$dcode)
    total_outlier_score <- df_tp$total_outlier_score_NEG
    a <- as.data.frame(cbind(dcode,total_outlier_score))
    t_test_out <- compareGroups(dcode ~ .,method = NA, data = a,alpha = 0.05)
    cat(paste0("P-value total outlier score NEG deviations (site=", sample_name," tp ",tp,"): ", t_test_out[["total_outlier_score"]][["p.overall"]],"\n"))
    normality <- attributes(t_test_out[[1]])$method[2]
    if (normality=="normal"){param <- "p" #parametric test
    } else {#non-normal: non-parametric test
      param <- "np"}
    a$dcode <- factor(a$dcode, labels = c('Healthy', 'Schizophrenia'))
    plot <- ggbetweenstats(data=a, x=dcode, y=total_outlier_score,
                           ylab = "Total negative outlier regions",
                           xlab="",
                           type=param, 
                           caption=paste0(tp_text))+
      scale_color_manual(values = c("#33A02C", "#1F78B4"))
    sub_1 <- plot[["plot_env"]][["subtitle"]][[2]][3]
    sub_1 <- gsub("\\(|\\)", "", sub_1)
    sub_2 <- plot[["labels"]][["subtitle"]][[4]][3]
    sub_2 <- gsub("\\(|\\)", "", sub_2)
    sub_3 <- plot[["labels"]][["subtitle"]][[5]][3]
    sub_3 <- gsub("\\(|\\)", "", sub_3)
    sub_3 <- as.numeric(gsub("[^0-9.-]", "", sub_3))
    sub_4 <- plot[["labels"]][["subtitle"]][[6]][2]
    sub_4 <- gsub("\\(|\\)", "", sub_4)
    sub_5 <- plot[["plot_env"]][["subtitle"]][[2]][2]
    sub_5 <- as.numeric(gsub("[^0-9.]", "", sub_5))
    if (normality=="normal"){plot <- plot + 
      labs(subtitle = ggplot2::expr(paste(italic("t")["Welch"],"(",!!sub_5,")"," = ", !!sub_1 ,
                                          ", ", widehat(italic("d"))["Cohen"], " = ", !!sub_2 ,
                                          ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]")))
    } else{plot <- plot + 
      labs(subtitle = ggplot2::expr(paste(italic("W")["Mann-Whitney"]," = ", !!trunc(as.numeric(sub_1)) ,
                                          ", ", widehat(italic("r"))["biserial"]^"rank", " = ", !!sub_2 ,
                                          ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]"
      )))}
    plot <- plot + theme(axis.text = element_text(size = 15),  
                         axis.title = element_text(size = 16),  
                         plot.subtitle = element_text(size = 15, face="italic"),
                         plot.title = element_text(size = 18))
    #plot <- plot + theme(text = element_text(size = 14),
    #                     axis.text = element_text(size = 17),  
    #                     axis.title = element_text(size = 17),  
    #                     plot.subtitle = element_text(size = 11),
    #                     plot.title = element_text(size = 18)) 
    png(file=paste0("",sample_name,"_tp",tp,".png"),width=6, height=4.5, units="in",res=600)
    print(plot)
    dev.off()
  }
}
rm(a,df,df_tp,plot,t_test_out,count,dcode,i,normality,out_controls,out_pats,param,
   percentage_out,sample,sample_name,total_outlier_score,tp,tps)
rm(out_controls_df,out_pats_df)

aaa <- read.csv("")
bbb <- read.csv("")
ccc <- cbind(aaa,bbb)

ddd <- read.csv("")
eee <- read.csv("")
fff <- cbind(ddd,eee)
#-----------------------------------
# - Longitudinal total outlier count
#-----------------------------------
# do it for each sample separately, remove from the list "all clinical samples" that are analyzed together
clinical_samples_long <- clinical_samples[-length(clinical_samples)] 

for (sample in clinical_samples_long) {
  df <- read.csv(z.score)
  cov_te <- df_cov_te
  df<- cbind(cov_te,df)
  df <- df[df$sample %in% sample,]
  tps <- unique(df$timepoint)
  if (length(tps)>1){
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    df <- bin_outlier_function(df, cov_te)
    
    # Longitudinal within-subject difference boxplot (for both groups)-POS and NEG deviations
    p1 <- ggwithinstats(
      data = df[df$dcode==0,],
      x    = timepoint,
      subject.id = subID,
      y    = total_outlier_score,
      ylab = "Total outlier count controls",
      type = "np",
      package = "RColorBrewer",
      palette = "Paired",
      caption=paste0("Sample: ",sample))
    p2 <- ggwithinstats(
      data = df[df$dcode==1,],
      x    = timepoint,
      subject.id = subID,
      y    = total_outlier_score,
      ylab = "Total outlier count patients",
      type = "np",
      package = "RColorBrewer",
      palette = "RdYlBu",
      caption=paste0("Sample: ",sample))
    plot <- combine_plots(
      plotlist = list(p1, p2),
      plotgrid.args = list(nrow = 1))
    plot
    png(file=paste0("withinsubject_total_outlier_count_",sample,".png"),width=11, height=6, units="in",res=600)
    print(plot)
    dev.off()
    
    # Longitudinal within-subject difference boxplot (for both groups)-POS deviations
    p1 <- ggwithinstats(
      data = df[df$dcode==0,],
      x    = timepoint,
      subject.id = subID,
      y    = total_outlier_score_POS,
      ylab = "Total outlier count supra-deviant controls",
      type = "np",
      package = "RColorBrewer",
      palette = "Paired",
      caption=paste0("Sample: ",sample))
    p2 <- ggwithinstats(
      data = df[df$dcode==1,],
      x    = timepoint,
      subject.id = subID,
      y    = total_outlier_score_POS,
      ylab = "Total outlier count supra-deviant patients",
      type = "np",
      package = "RColorBrewer",
      palette = "RdYlBu",
      caption=paste0("Sample: ",sample))
    plot <- combine_plots(
      plotlist = list(p1, p2),
      plotgrid.args = list(nrow = 1))
    plot
    png(file=paste0("withinsubject_total_outlier_count_POS_",sample,".png"),width=11, height=6, units="in",res=600)
    print(plot)
    dev.off()
    
    # Longitudinal within-subject difference boxplot (for both groups)-NEG deviations
    p1 <- ggwithinstats(
      data = df[df$dcode==0,],
      x    = timepoint,
      subject.id = subID,
      y    = total_outlier_score_NEG,
      ylab = "Total outlier count infra-deviant controls",
      type = "np",
      package = "RColorBrewer",
      palette = "Paired",
      caption=paste0("Sample: ",sample))
    p2 <- ggwithinstats(
      data = df[df$dcode==1,],
      x    = timepoint,
      subject.id = subID,
      y    = total_outlier_score_NEG,
      ylab = "Total outlier count infra-deviant patients",
      type = "np",
      package = "RColorBrewer",
      palette = "RdYlBu",
      caption=paste0("Sample: ",sample))
    plot <- combine_plots(
      plotlist = list(p1, p2),
      plotgrid.args = list(nrow = 1))
    plot
    png(file=paste0("withinsubject_total_outlier_count_NEG_",sample,".png"),width=11, height=6, units="in",res=600)
    print(plot)
    dev.off()
    
    #LONGITUDINAL CHANGE IN TOTAL_OUT_SCORE (first minus second timepoints)
    df_tp_1 <- df[df$timepoint == tps[1], ]
    total_tp1 <- df_tp_1$total_outlier_score 
    df_tp_2 <- df[df$timepoint == tps[2], ]
    total_tp2 <- df_tp_2$total_outlier_score
    
    total_change <- as.data.frame(cbind(total_tp1-total_tp2,df_tp_1$dcode, df_tp_1$age))# get baseline age
    names(total_change) <- c("total_change","group","age")
    total_change$group <- factor(total_change$group, labels = c('Healthy', 'Schizophrenia'))
    t_test <- compareGroups(group ~ total_change,method = NA, data = total_change,alpha = 0.05)
    print(paste0("P-value longitudinal change in total out score (", sample,"): ",t_test[["total_change"]][["p.overall"]]))
    normality <- attributes(t_test[[1]])$method[2]
    if (normality=="normal"){param <- "p" #parametric test
    } else {#non-normal: non-parametric test
      param <- "np"}
    plot <- ggbetweenstats(data=total_change, x=group, y=total_change,
                           ylab = "Longitudinal total outlier regions",
                           xlab="",
                           type=param)+
      scale_color_manual(values = c("#33A02C", "#1F78B4"))
    sub_1 <- plot[["plot_env"]][["subtitle"]][[2]][3]
    sub_1 <- gsub("\\(|\\)", "", sub_1)
    sub_2 <- plot[["labels"]][["subtitle"]][[4]][3]
    sub_2 <- gsub("\\(|\\)", "", sub_2)
    sub_3 <- plot[["labels"]][["subtitle"]][[5]][3]
    sub_3 <- gsub("\\(|\\)", "", sub_3)
    sub_3 <- as.numeric(gsub("[^0-9.-]", "", sub_3))
    sub_4 <- plot[["labels"]][["subtitle"]][[6]][2]
    sub_4 <- gsub("\\(|\\)", "", sub_4)
    sub_5 <- plot[["plot_env"]][["subtitle"]][[2]][2]
    sub_5 <- as.numeric(gsub("[^0-9.]", "", sub_5))
    if (normality=="normal"){plot <- plot + 
      labs(subtitle = ggplot2::expr(paste(italic("t")["Welch"],"(",!!sub_5,")"," = ", !!sub_1 ,
                                          ", ", widehat(italic("d"))["Cohen"], " = ", !!sub_2 ,
                                          ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]")))
    } else{plot <- plot + 
      labs(subtitle = ggplot2::expr(paste(italic("W")["Mann-Whitney"]," = ", !!trunc(as.numeric(sub_1)) ,
                                          ", ", widehat(italic("r"))["biserial"]^"rank", " = ", !!sub_2 ,
                                          ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]"
      )))}
    plot <- plot + theme(axis.text = element_text(size = 15),  
                         axis.title = element_text(size = 16),  
                         plot.subtitle = element_text(size = 15, face="italic"),
                         plot.title = element_text(size = 18))
    #plot <- plot + theme(text = element_text(size = 14),
    #                     axis.text = element_text(size = 17),  
    #                     axis.title = element_text(size = 17),  
    #                     plot.subtitle = element_text(size = 11),
    #                     plot.title = element_text(size = 18)) 
    png(file=paste0("",sample_name,".png"),width=6, height=4, units="in",res=600)
    print(plot)
    dev.off()
    
    # plot change in total outlier count
    # calculate 2std from mean (total out change in controls)
    upper <- mean(total_change[total_change["group"]=="Healthy",1]) + 2*(sd(total_change[total_change["group"]=="Healthy",1]))
    lower <- mean(total_change[total_change["group"]=="Healthy",1]) - 2*(sd(total_change[total_change["group"]=="Healthy",1]))
    total_change$upper <- upper
    total_change$lower <- lower
    plot <- ggplot(total_change, aes(age,total_change))+
      geom_point(aes(color=group))+
      geom_line(aes(age, upper))+
      geom_line(aes(age, lower))+
      #scale_color_brewer(palette = "Dark2")+
      theme_minimal()+
      labs(y = "Change in total outlier regions")+
      scale_color_manual(values = c("#33A02C", "#1F78B4"))
    plot <- plot + theme(text = element_text(size = 14),
                         axis.text = element_text(size = 16),  
                         axis.title = element_text(size = 17),  
                         plot.subtitle = element_text(size = 11),
                         plot.title = element_text(size = 18),
                         legend.text = element_text(size = 16),
                         legend.title=element_blank(),
                         axis.title.x = element_text(face = "bold"),
                         axis.title.y = element_text(face = "bold")) 
    png(file=paste0("longitudinal_total_outlier_regions_",sample,".png"),width=6.5, height=4.1, units="in",res=600)
    print(plot)
    dev.off()
  }
}
rm(clinical_samples_long,df,df_tp_1,df_tp_2,p1,p2,plot,total_change,lower,sample,total_tp1,total_tp2,tps,upper)

#-----------------------------------
# - - Proportion test reg outliers -
#-----------------------------------
#http://www.sthda.com/english/wiki/two-proportions-z-test-in-r

## warnings appear when running the proportion test function whenever there are 
## regions that both HC and PATS have 0 outliers. The output will be NAN.

prop_test_function <- function(out_HC, out_PAT, dataframe) {
  p_values <- x_squared <- numeric(length = length(out_HC))
  for (i in seq_along(out_HC)) {
    result <- prop.test(x = c(out_HC[i], out_PAT[i]),
                        n = c(length(which(dataframe$dcode == 0)), length(which(dataframe$dcode == 1))),
                        alternative = "two.sided")
    p_values[i] <- result$p.value
    x_squared[i] <- result$statistic[[1]]}
  return(data.frame(features=names(out_HC), p_values, x_squared))
}

for (sample in clinical_samples) {
  df <- df_Z_scores
  df<- cbind(df,df_cov_te)
  df <- df[df$sample %in% sample,]
  tps <- unique(df$timepoint)
  for (tp in tps) {
    # for all clinical samples, only look at baseline subjects
    # if you want to run it only on first timepoint subjects of all clinical samples uncomment the following line
    #if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    df_tp <- df[df$timepoint == tp, ]
    features <- setdiff(names(df_tp), names(df_cov_te))
    #Binarise outliers and create total_outlier score ACROSS total ROIs (score per participant)
    df3 <- as.data.frame(ifelse(df_tp[df_tp$dcode==0,1:62] > 1.96, 1, 0))
    out_pos_HC <- colSums(df3)
    df3 <- as.data.frame(ifelse(df_tp[df_tp$dcode==0,1:62] < -1.96, 1, 0))
    out_neg_HC <- colSums(df3)
    df3 <- as.data.frame(ifelse(df_tp[df_tp$dcode==1,1:62] > 1.96, 1, 0))
    out_pos_PAT <- colSums(df3)
    df3 <- as.data.frame(ifelse(df_tp[df_tp$dcode==1,1:62] < -1.96, 1, 0))
    out_neg_PAT <- colSums(df3)
    
    # Carry out proportions test for positive deviations HC and PATS
    res_stats_pos <- prop_test_function(out_pos_HC, out_pos_PAT, df_tp)
    res_stats_pos_signif <- res_stats_pos[res_stats_pos$p_values < 0.05 & complete.cases(res_stats_pos),]
    write.csv(res_stats_pos, file = paste0("Proportions_test_POS_",sample_name,"_tp",tp,".csv"),row.names=FALSE)
    if (nrow(res_stats_pos_signif) != 0){
      write.csv(res_stats_pos_signif, file = paste0("Proportions_test_POS_signif_",sample_name,"_tp",tp,".csv"),row.names=FALSE)}
    
    # Carry out proportions test for positive deviations HC and PATS
    res_stats_neg <- prop_test_function(out_neg_HC, out_neg_PAT, df_tp)
    res_stats_neg_signif <- res_stats_neg[res_stats_neg$p_values < 0.05 & complete.cases(res_stats_neg),]
    write.csv(res_stats_neg, file = paste0("Proportions_test_NEG_",sample_name,"_tp",tp,".csv"),row.names=FALSE)
    if (nrow(res_stats_neg_signif) != 0){
      write.csv(res_stats_neg_signif, file = paste0("Proportions_test_NEG_signif_",sample_name,"_tp",tp,".csv"),row.names=FALSE)}
  }
}
rm(df,df_tp,df3,res_stats_neg,res_stats_neg_signif,res_stats_pos,res_stats_pos_signif,
   out_neg_HC,out_neg_PAT,out_pos_HC,out_pos_PAT,sample,sample_name,tp,tps)

#-----------------------------------
# - - Two sample T-test z-scores - -
#-----------------------------------
for (sample in clinical_samples) {
  df <- df_Z_scores
  cov_te <- df_cov_te
  df<- cbind(cov_te,df)
  df <- df[df$sample %in% sample,]
  tps <- unique(df$timepoint)
  for (tp in tps) {
    # for all clinical samples, only look at baseline subjects
    # if you want to run it only on first timepoint subjects of all clinical samples uncomment the following line
    #if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    df_tp <- df[df$timepoint == tp, ]
    df_tp$dcode <- ifelse(df_tp$dcode == 0, 1, ifelse(df_tp$dcode == 1, 0, df_tp$dcode))
    dcode <- factor(df_tp$dcode, labels = c('Schizophrenia', 'Healthy'))
    a <- cbind(df_tp[,features],dcode)
    p_value_list <- list()
    t_statistic_list <- list()
    normality_test_list <- list()
    groups <- compareGroups(dcode ~ .,method = NA, data = a,alpha = 0.05)
    for (i in seq_along(features)) {
      normality <- attributes(groups[[i]])$method[2] # normal or non-normal
      if (normality == "normal") {
        # Parametric test (t-test)
        param <- "p" #parametric test
        t_test_result <- t.test(a[, i] ~ a$dcode)
        p_value <- t_test_result$p.value
        t_statistic <- t_test_result$statistic
      } else {
        # Non-parametric test (Man-Whitney U test)
        param <- "np" #non-parametric test
        wilcox_test_result <- wilcox.test(a[, i] ~ a$dcode)
        p_value <- wilcox_test_result$p.value
        t_statistic <- wilcox_test_result$statistic
      }
      p_value_list[[i]] <- p_value
      t_statistic_list[[i]] <- t_statistic
      dcode_column <- a$dcode
      reg_column <- a[features[i]]
      b <- as.data.frame(cbind(dcode_column,reg_column))
      b$dcode <- factor(b$dcode, labels = c('Schizophrenia', 'Healthy'))
      colnames(b) <- c("dcode_column","region")
      plot <- ggbetweenstats(data=b, x=dcode_column, y=region,
                             ylab = "Z-scores",
                             xlab="",
                             type=param, 
                             title=paste0("Region: ", features[i]," (sample: ", sample_name,", timepoint ",tp,")"))
    }
    fdr_p_value <- p.adjust(do.call(c, p_value_list), method = "fdr")
    p_value_list <- data.frame(feature = features, p_value = do.call(c, p_value_list),statistic = do.call(c, t_statistic_list),fdr_p = fdr_p_value)
    write.csv(p_value_list, file = paste0("Z_scores_p_t_",sample_name,"_tp",tp,".csv"),row.names=FALSE)
    p_value_list_signif <- p_value_list[p_value_list$fdr_p<0.05,]
    if (nrow(p_value_list_signif) != 0){
      write.csv(p_value_list_signif, file = paste0("Z_scores_p_t_signif_fdr_",sample_name,"_tp",tp,".csv"),row.names=FALSE)}
  }
} 

rm(a,b,df,df_tp,groups,p_value_list,p_value_list_signif,plot,reg_column,dcode,dcode_column,
   i,normality,p_value,param,sample,sample_name,tp,tps)

a1 <-read.csv(".")
a2<-read.csv(".")
a1fdr<-read.csv(".")
a2fdr<-read.csv(".")

# Do table
library(gridExtra)
library(grid)
library(gtable)
table <- tableGrob(a1[32:62,], theme=ttheme_minimal(), rows = NULL)
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 1, b = nrow(table), l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 1, l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 1, b=nrow(table),l = 4, r = 4)
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 1, b=nrow(table),l = 2, r = 2)
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 3, b=3 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 5, b=5 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 7, b=7 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 9, b=9 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 11, b=11 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 13, b=13 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 15, b=15 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 17, b=17 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 19, b=19 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 21, b=21 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 23, b=23 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 25, b=25 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 27, b=27 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 29, b=29 ,l = 1, r = ncol(table))
table <- gtable_add_grob(table,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 31, b=31 ,l = 1, r = ncol(table))
dev.off()
grid.newpage()
grid.draw(table)
dev.off()

# - Remove outliers and redo tests -
#-----------------------------------
df <- df_Z_scores
df3 <- as.data.frame(ifelse(df[,features] > 1.96 | df[,features] < -1.96, 1, 0))
results <- lapply(df3, function(column) {
  IDs <- df_cov_te$ID[which(column == 0)]
  return(IDs)})
names(results) <- names(df3)
df_filtrado <- list()
cov_te_filtrado <- list()
# remove outliers df
for (col in names(results)) {
  sujetos_quedar <- results[[col]]
  df_filtrado[[col]] <- df[df_cov_te$ID %in% sujetos_quedar, col]
}
# remove outliers cov_te
for (col in names(results)) {
  sujetos_quedar <- results[[col]]
  cov_te_filtrado[[col]] <- subset(df_cov_te, (ID %in% sujetos_quedar))
}

for (sample in clinical_samples) {
  filtered_data_frames <- list()
  filtered_value_lists <- list()
  for (i in seq_along(cov_te_filtrado)) {
    filtered_df <- filter(cov_te_filtrado[[i]], sample == "utrecht")
    filtered_values <- df_filtrado[[i]][filtered_df %>% row_number()]
    filtered_data_frames[[i]] <- filtered_df
    filtered_value_lists[[i]] <- filtered_values
  }
  tps <- unique(lapply(filtered_data_frames, function(df) unique(df$timepoint)))
  for (tp in tps[[1]]) {
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    filtered_data_frames_tp <- list()
    filtered_value_lists_tp <- list()
    for (i in seq_along(filtered_data_frames)) {
      filtered_df <- filter(filtered_data_frames[[i]], timepoint == tp)
      filtered_values <- filtered_value_lists[[i]][filtered_df %>% row_number()]
      filtered_data_frames_tp[[i]] <- filtered_df
      filtered_value_lists_tp[[i]] <- filtered_values
    }
    t_test_results <- list()
    for (col in 1:length(filtered_value_lists_tp)) {
      dcode_column <- filtered_data_frames_tp[[col]]$dcode
      Z_column <- filtered_value_lists_tp[[col]]
      a <- as.data.frame(cbind(dcode_column,Z_column))
      t_test <- compareGroups(dcode_column ~ Z_column, method = NA, data = a, alpha = 0.05)
      t_test_results[[col]] <- t_test[["Z_column"]][["p.overall"]]
    }
    fdr_p_value <- p.adjust(do.call(c, t_test_results), method = "fdr")
    t_test_results <- data.frame(features = features, p_value = do.call(c, t_test_results), fdr_p = fdr_p_value)
    write.csv(t_test_results, file = paste0("",sample_name,"_tp",tp,"_outliers_removed.csv"),row.names=FALSE)
    p_value_list_signif <- t_test_results[t_test_results$fdr_p<0.05,]
    if (nrow(p_value_list_signif) != 0){
      write.csv(p_value_list_signif, file = paste0("Z_scores_p_signif_fdr_",sample_name,"_tp",tp,"_outliers_removed.csv"),row.names=FALSE)}
    } 
  }
  
b1 <-read.csv("")
b2<-read.csv("")
b1fdr<-read.csv("")
b2fdr<-read.csv("")

#-----------------------------------
# - - - Linear regression - GAM  - -
#-----------------------------------
# age-quadratic, rest of covariates linear
# residualize for covariates using gam,
# get residuals and do a t-test on the residuals
# Read MS data
library(mgcv)

for (sample in clinical_samples) {
  df <- df_features
  cov_te <- df_cov_te
  df<- cbind(cov_te,df)
  df <- df[df$sample %in% sample,]
  tps <- unique(df$timepoint)
  for (tp in tps) {
    # for all clinical samples, only look at baseline subjects
    # if you want to run it only on first timepoint subjects of all clinical samples uncomment the following line
    #if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    df_tp <- df[df$timepoint == tp, ]
    features <- setdiff(names(df_tp), names(cov_te))
    
    # make scanner binary variables
    library(tidyverse)
    cov_te_scanners <- df_tp %>%
      mutate(value = 1) %>%
      pivot_wider(names_from=scanner, values_from=value, values_fill=0, names_prefix="scanner_")
    
    # Prepare covariates
    scanner_cols <- grep("^scanner", names(cov_te_scanners), value = TRUE)
    scanners <- cov_te_scanners[scanner_cols]
    colnames(scanners) <- paste0("scanner", 1:ncol(scanners))
    df_tp <- cbind(df_tp, scanners)
    df_tp <- subset(df_tp, select = -scanner) #remove scanner column
    
    columns_resid <- list()
    fitted_values_fit <- list()
    scanner_cols <- colnames(df_tp)[grepl("^scanner", colnames(df_tp))]
    formula <- paste("+ sex + euler_med +", paste(scanner_cols, collapse = " + "))
    fit <- lapply(df_tp[,features], function(x) gam(as.formula(paste("x ~ s(age, bs = 'cr', k = 5)", formula)), data = df_tp))
    
    for(i in seq_along(fit)){
      residuals<-fit[[i]]$residuals
      fitted_values_f<-fit[[i]]$fitted.values
      columns_resid[[i]]<-residuals
      fitted_values_fit[[i]]<-fitted_values_f
    }
    #plot(df$age,fitted_values_fit[[32]])
    
    residuals<-as.data.frame(columns_resid)
    names(residuals)<-paste0(features,"_resid")
    dcode<- df_tp$dcode
    residuals<- cbind(residuals,dcode)
    
    # using compare groups for t-test
    p_value_list <- list()
    groups <- compareGroups(dcode ~ .,method = NA, data = residuals,alpha = 0.05)
    for(i in seq_along(groups)){
      p_value<-groups[[i]]$p.overall
      p_value_list[[i]]<-p_value
    }
    fdr_p_value <- p.adjust(do.call(c, p_value_list), method = "fdr")
    p_value_list <- data.frame(feature = features, p_value = do.call(c, p_value_list),statistic = do.call(c, t_statistic_list),fdr_p = fdr_p_value)
    write.csv(p_value_list, file = paste0("GAM_p_",sample_name,"_tp",tp,".csv"),row.names=FALSE)
    p_value_list_signif <- p_value_list[p_value_list$fdr_p<0.05,]
    if (nrow(p_value_list_signif) != 0){
      write.csv(p_value_list_signif, file = paste0("GAM_p_signif_",sample_name,"_tp",tp,".csv"),row.names=FALSE)}
  }
}

rm(columns_resid,cov_te_scanners,df,df_tp,fit,fitted_values_fit,groups,p_value_list,
   p_value_list_signif,residuals,scanners,dcode,fitted_values_f,formula,i,p_value,
   sample,sample_name,scanner_cols,tp,tps)

#-----------------------------------
# - Person-Based Similarity Index  -
#-----------------------------------
# Person-Based Similarity Index (Janssen et al. 2021)


for (sample in clinical_samples) {
  # HC PBSI
  M_HC_df <- df_Z_scores[df_cov_te$sample %in% sample,]
  cov_te_HC <- df_cov_te[df_cov_te$sample %in% sample,]
  M_HC_df <- M_HC_df[cov_te_HC$dcode==0,]
  cov_te_HC <- cov_te_HC[cov_te_HC$dcode==0,]
  # PAT PBSI
  M_PT_df <- df_Z_scores[df_cov_te$sample %in% sample,]
  cov_te_PT <- df_cov_te[df_cov_te$sample %in% sample,]
  M_PT_df <- M_PT_df[cov_te_PT$dcode==1,]
  cov_te_PT <- cov_te_PT[cov_te_PT$dcode==1,]
  
  tps <- unique(cov_te_PT$timepoint)
  for (tp in tps) {
    # for all clinical samples, only look at baseline subjects
    # if you want to run it also on second timepoint subjects of all clinical samples comment the following line
    if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    M_HC <- M_HC_df[cov_te_HC$timepoint == tp, ]
    
    # # 1. Calculate PBSI_HC ---
    SIMS_HC <- cor(t(M_HC), use= "everything", method = "spearman")
    v_HC <- c()
    
    for (i in 1:ncol(SIMS_HC)){
      mat_HC <- sum(SIMS_HC[which(!is.na(SIMS_HC[,i])), i])
      v_HC <- c(v_HC, mat_HC)
    }
    
    PBSI_HC <- (v_HC-1)/((dim(SIMS_HC)[1])-1)
    PBSI_HC <- t(PBSI_HC)
    
    # # Check for outliers
    mHOMS_HC <- mean(PBSI_HC)
    sHOMS_HC <- sd(PBSI_HC)
    outliersN_HC <- mHOMS_HC - 2*sHOMS_HC # - 2SD
    outliersP_HC <- mHOMS_HC + 2*sHOMS_HC # + 2SD
    
    outliers_HC <- length(which(PBSI_HC < outliersN_HC)) +
      length(which(PBSI_HC > outliersP_HC))
    if (outliers_HC != 0){
      cat('\n Check your subjects, you seem to have', outliers_HC, 'HC outliers in', sample_name, 'timepoint', tp,'\n')
    }
    
    M_PT <- M_PT_df[cov_te_PT$timepoint == tp, ]
    
    # # 2. Calculate PBSI_PT ---
    SIMS_PT <- cor(t(M_PT), use= "everything", method = "spearman")
    v_PT <- c()
    
    for (j in 1:ncol(SIMS_PT)){
      mat_PT <- sum(SIMS_PT[which(!is.na(SIMS_PT[,j])), j])
      v_PT <- c(v_PT, mat_PT)
    }
    
    PBSI_PT <- (v_PT-1)/((dim(SIMS_PT)[1])-1)
    PBSI_PT <- t(PBSI_PT)
    
    # Check for outliers
    mHOMS_PT <- mean(PBSI_PT)
    sHOMS_PT <- sd(PBSI_PT)
    outliersN_PT <- mHOMS_PT - 2*sHOMS_PT # - 2SD
    outliersP_PT <- mHOMS_PT + 2*sHOMS_PT # + 2SD
    
    outliers_PT <- length(which(PBSI_PT < outliersN_PT)) +
      length(which(PBSI_PT > outliersP_PT))
    if (outliers_PT != 0){
      cat('\n Check your subjects, you seem to have', outliers_PT, 'PT outliers in', sample_name, 'timepoint', tp,'\n')
    }
    # # 3. Calculate PBSI_PT_HC ---
    PBSI_PT_HC <- SIMS_PT_HC <- c()
    for (k in 1:nrow(M_PT)){
      subj_PT <- M_PT[k,]
      SIMS_PT_HC1 <- cor(x = t(M_HC), y = t(subj_PT), use= "everything", method = "spearman")
      v_PT_HC <- c()
      SIMS_PT_HC <- cbind(SIMS_PT_HC, SIMS_PT_HC1)
      for (i in 1:ncol(SIMS_PT_HC1)){
        mat_PT_HC <- sum(SIMS_PT_HC1[which(!is.na(SIMS_PT_HC1[,i])), i])
        v_PT_HC <- c(v_PT_HC, mat_PT_HC)
      }
      
      PBSI_PT_HC <- c(PBSI_PT_HC, (v_PT_HC-1)/((dim(SIMS_PT_HC1)[1])-1))
    }
    PBSI_PT_HC <- t(PBSI_PT_HC)
    
    # Check for outliers
    mHOMS_PT_HC <- mean(PBSI_PT_HC)
    sHOMS_PT_HC <- sd(PBSI_PT_HC)
    outliersN_PT_HC <- mHOMS_PT_HC - 2*sHOMS_PT_HC # - 2SD
    outliersP_PT_HC <- mHOMS_PT_HC + 2*sHOMS_PT_HC # + 2SD
    
    outliers_PT_HC <- length(which(PBSI_PT_HC < outliersN_PT_HC)) +
      length(which(PBSI_PT_HC > outliersP_PT_HC))
    if (outliers_PT_HC != 0){
      cat('\n Check your subjects, you seem to have', outliers_PT_HC, 'PT-HC outliers in', sample_name, 'timepoint', tp,'\n')
    }
    
    PBSI_HC_2<-as.vector(PBSI_HC)
    PBSI_PT_2<-as.vector(PBSI_PT)
    PBSI_PT_HC_2<-as.vector(PBSI_PT_HC)
    
    # t-test pbsi hc- pbsi pats
    PBSI_HC <- as.data.frame(t(PBSI_HC))
    PBSI_HC$dcode <- 0
    PBSI_PT <- as.data.frame(t(PBSI_PT))
    PBSI_PT$dcode <- 1
    a <- rbind(PBSI_HC,PBSI_PT)
    groups <- compareGroups(dcode ~ .,method = NA, data = a,alpha = 0.05)
    cat(paste0("P-value PBSI_HC/PBSI_PAT (site=", sample_name, " tp ",tp,"): ", groups[["V1"]][["p.overall"]],"\n"))
  }
}
rm(a,cov_te_HC,cov_te_PT,groups,M_HC,M_HC_df,M_PT,M_PT_df,PBSI_HC,PBSI_PT,PBSI_PT_HC,
   SIMS_HC,SIMS_PT,SIMS_PT_HC,SIMS_PT_HC1,subj_PT,i,j,k,mat_HC,mat_PT,mat_PT_HC,
   mHOMS_HC,mHOMS_PT,mHOMS_PT_HC,outliers_HC,outliers_PT,outliers_PT_HC,outliersN_HC,
   outliersN_PT,outliersN_PT_HC,outliersP_HC,outliersP_PT,outliersP_PT_HC,PBSI_HC_2,
   PBSI_PT_2,PBSI_PT_HC_2,sample,sample_name,sHOMS_HC,sHOMS_PT,sHOMS_PT_HC,tp,tps,
   v_HC,v_PT,v_PT_HC)


#---------------------------------
# - - - - Hamming distance - - - -
#---------------------------------
library(e1071)

for (sample in clinical_samples) {
  df <- df_Z_scores
  cov_te <- df_cov_te
  df<- cbind(df_cov_te,df)
  df <- df[df$sample %in% sample,]
  tps <- unique(df$timepoint)
  for (tp in tps) {
    # for all clinical samples, only look at baseline subjects
    # if you want to run it only on first timepoint subjects of all clinical samples uncomment the following line
    #if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    df_tp <- df[df$timepoint == tp, ]
    df_tp <- bin_outlier_function(df_tp, cov_te)
    
    # - Outlier matrix plots
    ## Controls, using binerised z -scores
    controls_bin.mat <- as.matrix(df_tp %>% filter(dcode == 0) %>%
                                    select(contains("_bin")))
    
    ## SZ, using binerised z -scores
    pats_bin.mat <- as.matrix(df_tp %>% filter(dcode == 1) %>%
                                select(contains("_bin")))
    # Each cell is the Hamming distance between the binary (see above)
    # vectors of each subject, based on pairwise comparisons.
    ## Control distances
    #plot(hamming.distance(controls_bin.mat), border = NA)
    controls_hamming <- hamming.distance(controls_bin.mat)
    write.csv(controls_hamming, paste0("controls_hamming_",sample_name,"_tp_",tp,".csv"))
    
    ## PATS distances
    #plot (hamming.distance(MCI_bin.mat), border = NA)
    pats_hamming <- hamming.distance(pats_bin.mat)
    write.csv(pats_hamming, paste0("pats_hamming_",sample_name,"_tp_",tp,".csv"))
    
    # t-test controls/patients
    pats_hamming_mean <- as.data.frame(rowMeans(pats_hamming))
    names(pats_hamming_mean)<-"hamming"
    pats_hamming_mean$dcode <- 1
    controls_hamming_mean <- as.data.frame(rowMeans(controls_hamming))
    names(controls_hamming_mean)<-"hamming"
    controls_hamming_mean$dcode <- 0
    a <- rbind(controls_hamming_mean,pats_hamming_mean)
    groups <- compareGroups(dcode ~ .,method = NA, data = a,alpha = 0.05)
    cat(paste0("P-value hamming distance HC/PAT (site=", sample_name, " tp ",tp,"): ", groups[["hamming"]][["p.overall"]],"\n"))
    
    # Mean, median and standard deviation of Hamming distances
    ## Controls
    cat("Controls Hamming distance",sample_name,"timepoint",tp,"\nMean:",mean(hamming.distance(controls_bin.mat)),
        "\nMedian:",median(hamming.distance(controls_bin.mat)),"\nStandard Deviation:",sd(hamming.distance(controls_bin.mat)),"\n")
    
    ## Patients
    cat("Patients Hamming distance",sample_name,"timepoint",tp,"\nMean:",mean(hamming.distance(pats_bin.mat)),
        "\nMedian:",median(hamming.distance(pats_bin.mat)),"\nStandard Deviation:",sd(hamming.distance(pats_bin.mat)),"\n\n")
    
    # Plotting Hamming distance frequency density plot ---
    # Calculate density
    dc <- density(controls_hamming)
    dp <- density(pats_hamming)
    
    # Prepare data
    controls_hamming_g <- as.data.frame(controls_hamming)
    controls_hamming_g <-gather(controls_hamming_g)
    controls_hamming_g$a1 <-  c("Controls")
    
    pats_hamming_g <- as.data.frame(pats_hamming)
    pats_hamming_g <-gather(pats_hamming_g)
    pats_hamming_g$a1 <-  c("Patients")
    
    hamming_all <- rbind(controls_hamming_g,pats_hamming_g) #combine datasets
    rm(controls_hamming_g,pats_hamming_g) #remove dfs
    hamming_all <- hamming_all %>% dplyr::rename("Dissimilarity" = "value")
    hamming_all <- hamming_all %>% dplyr::rename("diagnosis" = "a1")
    
    hamming_all$diagnosis <- as.factor(hamming_all$diagnosis)
    
    plot <-  ggplot(hamming_all, aes(x=Dissimilarity, fill=diagnosis)) + #colour=diagnosis
      geom_density(adjust = 2,alpha=0.6)+
      scale_color_brewer(palette = "Dark2")+
      #scale_color_manual (values=c("Controls" =  "grey65" ,  "Patients" ="mediumturquoise"))+
      #scale_fill_manual (values=c("Controls" =  "grey65" ,  "Patients" ="mediumturquoise")) +
      labs(y= "Density", x= "Dissimilarity")+ 
      theme_light()
    png(file=paste0("Dissimilarity_",sample_name,"_tp",tp,".png"),width=9, height=6, units="in",res=600)
    print(plot)
    dev.off()
    rm(dc,dp)
  }
}

rm(a,controls_bin.mat,controls_hamming_mean,controls_hamming,df,df_tp,groups,hamming_all,
   pats_bin.mat,pats_hamming_mean,pats_hamming,plot,sample,sample_name,tp,tps)


#----------------------------------
# - - - - - - Clustering  - - - - -
#----------------------------------
library(NbClust)
library(reshape2)

for (sample in clinical_samples) {
  df <- df_Z_scores
  cov_te <- df_cov_te
  df<- cbind(df_cov_te,df)
  df <- df[df$sample %in% sample,]
  tps <- unique(df$timepoint)
  for (tp in tps) {
    # for all clinical samples, only look at baseline subjects
    # if you want to run it only on first timepoint subjects of all clinical samples uncomment the following line
    #if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    df_tp <- df[df$timepoint == tp, ]
    df_tp$dcode <- as.factor(df_tp$dcode)
    df_pats <- df_tp[df_tp$dcode== 1, features]
    set.seed(123)
    result <- NbClust(df_pats, diss = NULL,
                      min.nc = 2, max.nc = 10, method = "kmeans", index = "all")
    optimal_clusters <- length(unique(result$Best.partition))
    set.seed(123)
    kmeans_result <- kmeans(df_pats, centers = optimal_clusters)
    # Add cluster assignment to the original dataframe
    df_pats$cluster <- kmeans_result$cluster
    
    distance_matrix <- dist(df_pats[,features])
    distance_matrix <- as.matrix(distance_matrix)
    melted_cormat <- melt(as.matrix(distance_matrix))
    df_matrix <- distance_matrix %>% 
      as.matrix() %>% 
      as_tibble(rownames="A") %>% 
      pivot_longer(-A,names_to="B", values_to="distances")
    df_matrix$A_ordered <- factor(df_matrix$A)
    df_matrix$B_ordered <- factor(df_matrix$B)
    df_matrix[df_matrix == 0] <- NA
    plot1 <- df_matrix %>% ggplot(aes(x = A_ordered, y = B_ordered , fill = distances)) +
      geom_tile()+
      scale_fill_gradient(
        low = "purple",
        high = "yellow",
        na.value = "white",
        limits = c(0, max(df_matrix$distances)))+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank())
    #plot <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    #  geom_tile()
    png(file=paste0("Distance_matrix_",sample_name,"_tp",tp,".png"),width=11, height=6, units="in",res=600)
    print(plot1)
    dev.off()
    
    for (c_num in 1:optimal_clusters) {
      cluster_mean <- colMeans(df_pats[df_pats$cluster == c_num, ])[1:(length(features))]
      write.csv(cluster_mean, file = paste0("cluster", c_num, "_mean_", sample_name, "_tp_", tp, ".csv"))
    }
  }
}

rm(df,df_pats,df_tp,df_matrix,distance_matrix,kmeans_result,melted_cormat,plot1,result,
   cluster1_mean,cluster2_mean,optimal_clusters,sample,sample_name,tp,tps,c_num,cluster_mean)

#----------------------------------
# - Relationship clinical variables 
#----------------------------------
library(tidyr)
# my clinical data is in this file, but it should be in cov_te
clinical <- read.csv("")
#clinical <- clinical %>% filter(ID %in% df$ID)
clinical_samples="utrecht"
for (sample in clinical_samples) {
  df <- df_Z_scores
  cov_te <- df_cov_te
  df<- cbind(cov_te,df)
  df <- df[df$sample %in% sample,]
  df <- df[df$dcode==1,]
  clinical <- clinical[clinical$dcode==1,]
  df <- cbind(clinical[,3:100], df[,13:74])
  tps <- unique(df$timepoint)
  for (tp in tps) {
    print(tp)
    # for all clinical samples, only look at baseline subjects
    # if you want to run it also on second timepoint subjects of all clinical samples comment the following line
    if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    df_tp <- df[df$timepoint == tp, ]
    df_tp$dcode <- as.factor(df_tp$dcode)
    
    df_long <- gather(df_tp, region, Z.score, features[1]:features[length(features)], factor_key=TRUE)
    
    xvars <- c("IQ_total", "PANSS_total", "Illness_duration_at_scan", "Total_P",
               "Total_N", "Total_G")
    for (xvar in xvars){
      plot <- grouped_ggscatterstats(
        data = filter(df_long, region %in% features),
        x = !!sym(xvar),
        y = Z.score,
        grouping.var = region,
        xlab = xvar)
      for (i in c(1:length(features))) {
        corr <- extract_stats(plot[[i]])$subtitle_data$estimate
        if (abs(corr)> 0.3){
          print(paste0("Correlation region ",i," ",features[i]," ",xvar,": ",corr))
          print(plot[[i]])
          png(file=paste0("Relationship_Z_scores_",xvar,"_",features[i],"_",sample,"_tp",tp,".png"),width=11, height=6, units="in",res=600)
          dev.off()
        }
      }
    }
    
    # Change region and clinical variable that you desire to explore manually
    # without covariables
    #model <- lm(rh_precentral ~ Illness_duration_at_scan,data=df_tp)
    #summary(model)
    
    # with age as covariable
    #model <- lm(rh_precentral ~ Illness_duration_at_scan + age,data=df_tp)
    #summary(model)
  }
}

rm(clinical,df,df_long,df_tp,plot,corr,i,sample,sample_name,tp,tps,xvar,xvars)


#----------------------------------
# - - - - Atlases overlap - - - - -
#----------------------------------
# There are n number of classes in the atlases. 
# Each feature from the fsaverage_*_aparc belongs to 1 of these classes.
# We decide class-membership based on the largest overlap. 
# Look for each column-maximum in the file fsavergae_*_aparc and 
# then in the row find the class.

# Classes' numbers are not paired to original classes numbers, check class column to see original class name

# Change to desired atlas name. Economo must be named "economo", Yeo must be named "yeo"
atlas_name <- "economo" 
# Change to desired atlas
atlas <- read.csv("") 

# Prepare feature names
if(atlas_name!="yeo"){
  coln <- sub("^ctx\\.", "", colnames(atlas))  # Remove prefix "ctx." from feature names
  coln <- gsub("\\.", "_", coln)  # Replace dots by underscores to match feature names in my MS and Z-scores df
}else{
  atlas <- t(atlas)
  colnames(atlas)<-atlas[1,]
  atlas <- atlas[-1, ]
  atlas <- cbind(class = rownames(atlas), atlas)
  rownames(atlas) <-NULL
  atlas <-as.data.frame(atlas)
  coln <- sub("^ctx\\-", "", colnames(atlas))  # Remove prefix "ctx-"
  coln <- gsub("\\-", "_", coln)  # Replace "-" by "_" to match feature names "rh_insula", "lh_insula"
}
colnames(atlas) <- coln
names(atlas)

# Read dataframe to match feature names in the form of "rh_caudalanteriorcingulate", "lh_caudalanteriorcingulate"
df <- read.csv("")#[,-1]

# Get regions in atlas present in MS file
#!! adjust column numbers to get feature names of MS file
regions <- intersect(colnames(atlas), colnames(df[features])) 
atlas <- atlas[, c(1,match(regions, colnames(atlas)))] 

# get maximum per column (feature) and assign a class
i_max <- sapply(atlas[, -1], which.max)

if(atlas_name!="economo"){
  # assign each class a number (rh and lh symmetric classes)
  atlas$fsaverage_atlas_aparc_num <- rep(1:(nrow(atlas)/2), times = 2)
  atlasoverlap <- data.frame(class = atlas[i_max, 1])
  atlasoverlap$class_num <- atlas[i_max, ncol(atlas)]
}else{
  atlasoverlap <- data.frame(class = atlas[i_max, 1])
  atlasoverlap$class_num <- gsub(".*-CT(\\d+).*", "\\1", atlasoverlap$class)
  # !! to be noted: class 1_2 will be assigned number 6 because for posterior analyses we need each class to be a number
  atlasoverlap$class_num <- ifelse(atlasoverlap$class_num == 1 & !endsWith(atlasoverlap$class, "1"), "6", atlasoverlap$class_num)
}
atlasoverlap$region<-regions

write.csv(atlasoverlap,paste0("",atlas_name,"_overlap.csv"))

### create files for github upload ###
economo <- read.csv("")[,-1]
economo <- economo[,-1]
economo$class_num <- gsub(6, "1_2", economo$class_num)
economo$region <- factor(economo$region, levels = features)
economo <- economo[order(economo$region), ]
colnames(economo)[2] <- "class"
#economo <- economo[, c(2,1)]

yeo <- read.csv("")[,-1]
yeo <- yeo[-1]
colnames(yeo)[1] <- "network"
yeo$region <- factor(yeo$region, levels = features)
yeo <- yeo[order(yeo$region), ]
yeo <- yeo[, c(2,1)]

write.csv(economo, file = "" ,row.names=FALSE)
write.csv(yeo, file = "" ,row.names=FALSE)


#--------------------------------------------------------
# Differences in Z-scores or features within atlas classes 
#---------------------------------------------------------
library(effectsize)
library(rlang)
atlas_name <- "economo"
atlas_overlap <- read.csv("")[,-1] #change atlas !!

for (sample in clinical_samples) {
  df <- df_Z_scores  # change this by df_features if you want to see differences in features instead of diff in Z-scores
  df <- cbind(df,df_cov_te)
  df <- df[df$sample %in% sample,]
  tps <- unique(df$timepoint)
  i=0
  for (tp in tps) {
    # for all clinical samples, only look at baseline subjects
    # if you want to run it only on first timepoint subjects of all clinical samples uncomment the following line
    #if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    df_tp <- df[df$timepoint == tp, ]
    
    networks <- atlas_overlap$class_num
    z_scores <- df_tp[,features]
    classes_total <- max(atlas_overlap$class_num)
    p_value_list <- list()
    t_statistic_list <- list()
    cohens_list <- list()
    for (class in 1:classes_total) {
      myregions <- which(networks == class)
      if (is_empty(myregions)) {
        print(paste0("Sample ", sample_name," timepoint ",tp," Class ",class," not present"))
        next}  # Skip to the next iteration
      
      nsubs <- nrow(df_tp)
      classreg <- numeric(nsubs)
      # Average of Z-scores/features across all regions that belong to the same class, for each subject
      for (subj in 1:nsubs) {
        classreg[subj] <- mean(as.numeric(z_scores[subj, myregions]))}
      
      df_tp$dcode <- ifelse(df_tp$dcode == 0, 1, ifelse(df_tp$dcode == 1, 0, df_tp$dcode))
      dcode <- factor(df_tp$dcode, labels = c('Schizophrenia', 'Healthy'))
      a <- data.frame(cbind(classreg,dcode))
      groups <- compareGroups(dcode ~ .,method = NA, data = a,alpha = 0.05)
      normality <- attributes(groups[[1]])$method[2] # normal or non-normal
      if (normality == "normal") {
        # Parametric test (t-test)
        param <- "p" #parametric test
        t_test_result <- t.test(a$classreg ~ a$dcode)
        cohens <- cohens_d(a$classreg ~ a$dcode)$Cohens_d
        p_value <- t_test_result$p.value
        t_statistic <- t_test_result$statistic
      } else {
        # Non-parametric test (Man-Whitney U test)
        param <- "np" #non-parametric test
        wilcox_test_result <- wilcox.test(a$classreg ~ a$dcode)
        cohens <- cohens_d(a$classreg ~ a$dcode)$Cohens_d
        p_value <- wilcox_test_result$p.value
        t_statistic <- wilcox_test_result$statistic
      }
        
      p_value_list[[class]] <- p_value
      t_statistic_list[[class]] <- t_statistic
      cohens_list[[class]] <- cohens
      #print(paste0("P-value features (", sample_name," timepoint ",tp,") class ",class,": ",groups[["classreg"]][["p.overall"]]))
      #p_value <- groups[["classreg"]][["p.overall"]]
      #p_value_list[[class]]<-p_value
      
      if(class==1){class_name = "Agranular"}
      if(class==2){class_name = "Frontal"}
      if(class==3){class_name = "Parietal"}
      if(class==4){class_name = "Polar"}
      if(class==5){class_name = "Granulous"}
      if(class==6){class_name = "Transitional"}
      
      dcode_column <- a$dcode
      reg_column <- a$classreg
      b <- as.data.frame(cbind(dcode_column,reg_column))
      b$dcode_column <- factor(b$dcode_column, labels = c('Schizophrenia', 'Healthy'))
      colnames(b) <- c("dcode_column","network")
      if (tp==1){tp_text="Baseline"}else{tp_text="Follow-up"}
      plot <- ggbetweenstats(data=b, x=dcode_column, y=network,
                             ylab = "Z-scores",
                             xlab="",
                             type=param, 
                             effsize.type="biased",
                             bf.message = FALSE,
                             centrality.label.args = list(size = 4, nudge_x = 0.4, segment.linetype = 4,min.segment.length = 0),
                             title=paste0(class_name," Class"," (",tp_text,")"),
                             ggplot.component = list(ggplot2::scale_y_continuous(
                               breaks = seq(-1.5, 2, 0.5),
                               limits = (c(-1.5,2)))))+
        scale_color_manual(values = c("#1F78B4","#33A02C"))
                            #ggplot.component = list(theme(plot.subtitle = element_text(size = 24, face = "bold"))))
      sub_1 <- plot[["plot_env"]][["subtitle"]][[2]][3]
      sub_1 <- gsub("\\(|\\)", "", sub_1)
      sub_2 <- plot[["labels"]][["subtitle"]][[4]][3]
      sub_2 <- gsub("\\(|\\)", "", sub_2)
      sub_3 <- plot[["labels"]][["subtitle"]][[5]][3]
      sub_3 <- gsub("\\(|\\)", "", sub_3)
      sub_3 <- as.numeric(gsub("[^0-9.-]", "", sub_3))
      sub_4 <- plot[["labels"]][["subtitle"]][[6]][2]
      sub_4 <- gsub("\\(|\\)", "", sub_4)
      sub_5 <- plot[["plot_env"]][["subtitle"]][[2]][2]
      sub_5 <- as.numeric(gsub("[^0-9.]", "", sub_5))
      if (normality=="normal"){plot <- plot + 
        labs(subtitle = ggplot2::expr(paste(italic("t")["Welch"],"(",!!sub_5,")"," = ", !!sub_1 ,
                                            ", ", widehat(italic("d"))["Cohen"], " = ", !!sub_2 ,
                                            ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]")))
      } else{plot <- plot + 
        labs(subtitle = ggplot2::expr(paste(italic("W")["Mann-Whitney"]," = ", !!trunc(as.numeric(sub_1)) ,
                                            ", ", widehat(italic("r"))["biserial"]^"rank", " = ", !!sub_2 ,
                                            ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]"
        )))}
      plot <- plot + theme(axis.text = element_text(size = 15),  
                           axis.title = element_text(size = 16),  
                           plot.subtitle = element_text(size = 15, face="italic"),
                           plot.title = element_text(size = 18)) 
      png(file = paste0("",class,"_",sample_name,"_tp",tp,".png"), width=6, height=4.5, units="in", res=300)
      print(plot)
      dev.off()
    }
    fdr_p_value <- p.adjust(do.call(c, p_value_list), method = "fdr")
    p_value_list <- data.frame(class = 1:classes_total, p_value = do.call(c, p_value_list),statistic = do.call(c, t_statistic_list), fdr_p=fdr_p_value, cohens_d=do.call(c, cohens_list))
    write.csv(p_value_list, file = paste0("",atlas_name,"_atlas_",sample_name,"_tp",tp,".csv"),row.names=FALSE)
    p_value_list_signif <- p_value_list[p_value_list$p_value<0.05 & complete.cases(p_value_list),]
    if (nrow(p_value_list_signif) > 0){
      write.csv(p_value_list_signif, file = paste0("",atlas_name,"_atlas_",sample_name,"_tp",tp,".csv"),row.names=FALSE)}
  }
}
a2  <- read.csv("")
b4 <- read.csv("")
c4<- read.csv("")
d4 <- read.csv("")

atlas_name <- "yeo"
atlas_overlap <- read.csv("")[,-1] #change atlas !!

for (sample in clinical_samples) {
  df <- df_Z_scores  # change this by df_features if you want to see differences in features instead of diff in Z-scores
  df <- cbind(df,df_cov_te)
  df <- df[df$sample %in% sample,]
  tps <- unique(df$timepoint)
  i=0
  for (tp in tps) {
    # for all clinical samples, only look at baseline subjects
    # if you want to run it only on first timepoint subjects of all clinical samples uncomment the following line
    #if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    df_tp <- df[df$timepoint == tp, ]
    
    networks <- atlas_overlap$class_num
    z_scores <- df_tp[,features]
    classes_total <- max(atlas_overlap$class_num)
    p_value_list <- list()
    t_statistic_list <- list()
    cohens_list <- list()
    for (class in 1:classes_total) {
      myregions <- which(networks == class)
      if (is_empty(myregions)) {
        print(paste0("Sample ", sample_name," timepoint ",tp," network ",class," not present"))
        next}  # Skip to the next iteration
      
      nsubs <- nrow(df_tp)
      classreg <- numeric(nsubs)
      # Average of Z-scores/features across all regions that belong to the same class, for each subject
      for (subj in 1:nsubs) {
        classreg[subj] <- mean(as.numeric(z_scores[subj, myregions]))}
      
      df_tp$dcode <- ifelse(df_tp$dcode == 0, 2, ifelse(df_tp$dcode == 1, 1, df_tp$dcode))
      dcode <- factor(df_tp$dcode, labels = c('Schizophrenia', 'Healthy'))
      a <- data.frame(cbind(classreg,dcode))
      a$dcode <- factor(a$dcode, labels = c('Schizophrenia', 'Healthy'))
      groups <- compareGroups(dcode ~ .,method = NA, data = a,alpha = 0.05)
      normality <- attributes(groups[[1]])$method[2] # normal or non-normal
      if (normality == "normal") {
        # Parametric test (t-test)
        param <- "p" #parametric test
        t_test_result <- t.test(a$classreg ~ a$dcode)
        cohens <- cohens_d(a$classreg ~ a$dcode)$Cohens_d
        p_value <- t_test_result$p.value
        t_statistic <- t_test_result$statistic
      } else {
        # Non-parametric test (Man-Whitney U test)
        param <- "np" #non-parametric test
        wilcox_test_result <- wilcox.test(a$classreg ~ a$dcode)
        cohens <- cohens_d(a$classreg ~ a$dcode)$Cohens_d
        p_value <- wilcox_test_result$p.value
        t_statistic <- wilcox_test_result$statistic
      }
      
      p_value_list[[class]] <- p_value
      t_statistic_list[[class]] <- t_statistic
      cohens_list[[class]] <- cohens
      
      if(class==1){class_name = "Visual"}
      if(class==2){class_name = "Somato Motor"}
      if(class==3){class_name = "Dorsal Attention"}
      if(class==4){class_name = "Ventral Attention"}
      if(class==5){class_name = "Limbic"}
      if(class==6){class_name = "Fronto Parietal"}
      if(class==7){class_name = "Default Mode"}
      
      dcode_column <- a$dcode
      dcode_column <- factor(dcode_column, labels = c('Schizophrenia', 'Healthy'))
      reg_column <- a$classreg
      b <- as.data.frame(cbind(dcode_column,reg_column))
      b$dcode_column <- factor(b$dcode_column, labels = c('Schizophrenia', 'Healthy'))
      colnames(b) <- c("dcode_column","network")
      if (tp==1){tp_text="Baseline"}else{tp_text="Follow-up"}
      plot <- ggbetweenstats(data=b, x=dcode_column, y=network,
                             ylab = "Z-scores",
                             xlab="",
                             type=param, 
                             effsize.type="biased",
                             bf.message = FALSE,
                             centrality.label.args = list(size = 4, nudge_x = 0.4, segment.linetype = 4,min.segment.length = 0),
                             title=paste0(class_name," Network"," (",tp_text,")"),
                             ggplot.component = list(ggplot2::scale_y_continuous(
                               breaks = seq(-1.5, 1.5, 0.5),
                               limits = (c(-1.5, 1.5)))))+
        scale_color_manual(values = c("#1F78B4","#33A02C"))
                             #ggplot.component = list(theme(plot.subtitle = element_text(size = 24, face = "bold"))))
      sub_1 <- plot[["plot_env"]][["subtitle"]][[2]][3]
      sub_1 <- gsub("\\(|\\)", "", sub_1)
      sub_2 <- plot[["labels"]][["subtitle"]][[4]][3]
      sub_2 <- gsub("\\(|\\)", "", sub_2)
      sub_3 <- plot[["labels"]][["subtitle"]][[5]][3]
      sub_3 <- gsub("\\(|\\)", "", sub_3)
      sub_3 <- as.numeric(gsub("[^0-9.-]", "", sub_3))
      sub_4 <- plot[["labels"]][["subtitle"]][[6]][2]
      sub_4 <- gsub("\\(|\\)", "", sub_4)
      sub_5 <- plot[["plot_env"]][["subtitle"]][[2]][2]
      sub_5 <- as.numeric(gsub("[^0-9.]", "", sub_5))
      if (normality!="normal"){plot <- plot + 
      labs(subtitle = ggplot2::expr(paste(italic("W")["Mann-Whitney"]," = ", !!trunc(as.numeric(sub_1)) ,
                                            ", ", widehat(italic("r"))["biserial"]^"rank", " = ", !!sub_2 ,
                                            ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]"
      )))} else{plot <- plot + 
        labs(subtitle = ggplot2::expr(paste(italic("t")["Welch"],"(",!!sub_5,")"," = ", !!sub_1 ,
                                            ", ", widehat(italic("d"))["Cohen"], " = ", !!sub_2 ,
                                            ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]")))
      }
      plot <- plot + theme(axis.text = element_text(size = 15),  
                           axis.title = element_text(size = 16),  
                           plot.subtitle = element_text(size = 15, face="italic"),
                           plot.title = element_text(size = 18)) 
      png(file = paste0("",class,"_",sample_name,"_tp",tp,".png"), width=6, height=4.5, units="in", res=300)
      print(plot)
      dev.off()
    }
    fdr_p_value <- p.adjust(do.call(c, p_value_list), method = "fdr")
    p_value_list <- data.frame(class = 1:classes_total, p_value = do.call(c, p_value_list),statistic = do.call(c, t_statistic_list), fdr_p=fdr_p_value, cohens_d=do.call(c, cohens_list))
    write.csv(p_value_list, file = paste0("",atlas_name,"_atlas_",sample_name,"_tp",tp,".csv"),row.names=FALSE)
    p_value_list_signif <- p_value_list[p_value_list$p_value<0.05 & complete.cases(p_value_list),]
    if (nrow(p_value_list_signif) != 0){
      write.csv(p_value_list_signif, file = paste0("",atlas_name,"_atlas_",sample_name,"_tp",tp,".csv"),row.names=FALSE)}
  }
}

rm(a,atlas_name,atlas_overlap,df,df_tp,p_value_list,p_value_list_signif,z_scores,atlas_name,class,classes_total,classreg,
   dcode,i,myregions,networks,nsubs,sample,sample_name,subj,tp,tps)

w <- read.csv("")
f  <- read.csv("")
g <- read.csv("")
h <- read.csv("")


#--------------------------------------------------------
# Differences in MS or features within atlas classes GAM
#---------------------------------------------------------
library(mgcv)
atlas_name <- "economo"
atlas_overlap <- read.csv("")[,-1] #change atlas !!

for (sample in clinical_samples) {
  df <- df_features 
  df <- cbind(df,df_cov_te)
  df <- df[df$sample %in% sample,]
  tps <- unique(df$timepoint)
  i=0
  for (tp in tps) {
    # for all clinical samples, only look at baseline subjects
    # if you want to run it only on first timepoint subjects of all clinical samples uncomment the following line
    #if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    df_tp <- df[df$timepoint == tp, ]
    MS_raw <- df_tp[,features]
    networks <- atlas_overlap$class_num
    classes_total <- max(atlas_overlap$class_num)
    p_value_list <- list()
    t_statistic_list <- list()
    cohens_list <- list()
    
    # make scanner binary variables
    library(tidyverse)
    cov_te_scanners <- df_tp %>%
      mutate(value = 1) %>%
      pivot_wider(names_from=scanner, values_from=value, values_fill=0, names_prefix="scanner_")
    
    # Prepare covariates
    scanner_cols <- grep("^scanner", names(cov_te_scanners), value = TRUE)
    scanners <- cov_te_scanners[scanner_cols]
    colnames(scanners) <- paste0("scanner", 1:ncol(scanners))
    df_tp <- cbind(df_tp, scanners)
    df_tp <- subset(df_tp, select = -scanner) #remove scanner column
    
    columns_resid <- list()
    fitted_values_fit <- list()
    scanner_cols <- colnames(df_tp)[grepl("^scanner", colnames(df_tp))]
    formula <- paste("+ sex + euler_med +", paste(scanner_cols, collapse = " + "))
    for (class in 1:classes_total) {
      myregions <- which(networks == class)
      if (is_empty(myregions)) {
        print(paste0("Sample ", sample_name," timepoint ",tp," Class ",class," not present"))
        next}  # Skip to the next iteration
      
      nsubs <- nrow(df_tp)
      classreg <- numeric(nsubs)
      # Average features across all regions that belong to the same class, for each subject
      for (subj in 1:nsubs) {
        classreg[subj] <- mean(as.numeric(MS_raw[subj, myregions]))}
      
      cat("Class ",class)
      reg <- cbind(df_tp[,63:74], classreg)
      controls <- reg[reg$dcode==0,"classreg"]
      pats <- reg[reg$dcode==1,"classreg"]
      print("mean hc ")
      print(mean(controls$classreg))
      print("sd hc ")
      print(sd(controls$classreg))
      print("mean pats ")
      print(mean(pats$classreg))
      print("sd pats ")
      print(sd(pats$classreg))
      fit <- gam(as.formula(paste("classreg ~ s(age, bs = 'cr', k = 5)", formula)), data = reg)
      
      
      residuals<-fit$residuals
      fitted_values_f<-fit$fitted.values
      columns_resid<-residuals
      fitted_values_fit<-fitted_values_f
      
      #plot(df$age,fitted_values_fit[[32]])
      
      residuals<-as.data.frame(columns_resid)
      #names(residuals)<-paste0(features,"_resid")
      df_tp$dcode <- ifelse(df_tp$dcode == 0, 1, ifelse(df_tp$dcode == 1, 0, df_tp$dcode))
      dcode <- factor(df_tp$dcode, labels = c('Schizophrenia', 'Healthy'))
      residuals<- cbind(residuals,dcode)
      
      # using compare groups for t-test
      groups <- compareGroups(dcode ~ .,method = NA, data = residuals,alpha = 0.05)
      for(i in seq_along(groups)){
        p_value<-groups[[i]]$p.overall
        p_value_list[[i]]<-p_value
      }
      groups <- compareGroups(dcode ~ .,method = NA, data = residuals,alpha = 0.05)
      normality <- attributes(groups[[1]])$method[2] # normal or non-normal
      if (normality == "normal") {
        # Parametric test (t-test)
        param <- "p" #parametric test
        t_test_result <- t.test(residuals$columns_resid ~ residuals$dcode)
        cohens <- cohens_d(residuals$columns_resid ~ residuals$dcode)$Cohens_d
        p_value <- t_test_result$p.value
        t_statistic <- t_test_result$statistic
      } else {
        # Non-parametric test (Man-Whitney U test)
        param <- "np" #non-parametric test
        wilcox_test_result <- wilcox.test(residuals$columns_resid ~ residuals$dcode)
        cohens <- cohens_d(residuals$columns_resid ~ residuals$dcode)$Cohens_d
        p_value <- wilcox_test_result$p.value
        t_statistic <- wilcox_test_result$statistic
      }
      
      p_value_list[[class]] <- p_value
      t_statistic_list[[class]] <- t_statistic
      cohens_list[[class]] <- cohens
      
      if(class==1){class_name = "Agranular"}
      if(class==2){class_name = "Frontal"}
      if(class==3){class_name = "Parietal"}
      if(class==4){class_name = "Polar"}
      if(class==5){class_name = "Granulous"}
      if(class==6){class_name = "Transitional"}
      
      dcode_column <- residuals$dcode
      reg_column <- residuals$columns_resid
      b <- as.data.frame(cbind(dcode_column,reg_column))
      b$dcode_column <- factor(b$dcode_column, labels = c('Schizophrenia', 'Healthy'))
      colnames(b) <- c("dcode_column","network")
      if (tp==1){tp_text="Baseline"}else{tp_text="Follow-up"}
      plot <- ggbetweenstats(data=b, x=dcode_column, y=network,
                             ylab = "MS",
                             xlab="",
                             type=param, 
                             effsize.type="biased",
                             bf.message = FALSE,
                             centrality.label.args = list(size = 4, nudge_x = 0.4, segment.linetype = 4,min.segment.length = 0),
                             title=paste0(class_name," Class"," (",tp_text,")"),
                             ggplot.component = list(ggplot2::scale_y_continuous(
                               breaks = seq(-0.05, 0.05, 0.05),
                               limits = (c(-0.05, 0.05)))))+
        scale_color_manual(values = c("#1F78B4","#33A02C"))
                             #ggplot.component = list(theme(plot.subtitle = element_text(size = 24, face = "bold"))))
      sub_1 <- plot[["plot_env"]][["subtitle"]][[2]][3]
      sub_1 <- gsub("\\(|\\)", "", sub_1)
      sub_2 <- plot[["labels"]][["subtitle"]][[4]][3]
      sub_2 <- gsub("\\(|\\)", "", sub_2)
      sub_3 <- plot[["labels"]][["subtitle"]][[5]][3]
      sub_3 <- gsub("\\(|\\)", "", sub_3)
      sub_3 <- as.numeric(gsub("[^0-9.-]", "", sub_3))
      sub_4 <- plot[["labels"]][["subtitle"]][[6]][2]
      sub_4 <- gsub("\\(|\\)", "", sub_4)
      sub_5 <- plot[["plot_env"]][["subtitle"]][[2]][2]
      sub_5 <- as.numeric(gsub("[^0-9.]", "", sub_5))
      if (normality!="normal"){plot <- plot + 
        labs(subtitle = ggplot2::expr(paste(italic("W")["Mann-Whitney"]," = ", !!trunc(as.numeric(sub_1)) ,
                                            ", ", widehat(italic("r"))["biserial"]^"rank", " = ", !!sub_2 ,
                                            ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]"
        )))} else{plot <- plot + 
          labs(subtitle = ggplot2::expr(paste(italic("t")["Welch"],"(",!!sub_5,")"," = ", !!sub_1 ,
                                              ", ", widehat(italic("d"))["Cohen"], " = ", !!sub_2 ,
                                              ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]")))
        }
      plot <- plot + theme(axis.text = element_text(size = 15),  
                           axis.title = element_text(size = 16),  
                           plot.subtitle = element_text(size = 15, face="italic"),
                           plot.title = element_text(size = 18))  
      png(file = paste0("",class,"_",sample_name,"_tp",tp,"_MS.png"), width=6, height=4.5, units="in", res=300)
      print(plot)
      dev.off()
    }
    fdr_p_value <- p.adjust(do.call(c, p_value_list), method = "fdr")
    p_value_list <- data.frame(network = 1:classes_total, p_value = do.call(c, p_value_list),statistic = do.call(c, t_statistic_list),fdr_p = fdr_p_value, cohens_d=do.call(c, cohens_list))
    #write.csv(p_value_list, file = paste0("",sample_name,"_tp",tp,".csv"),row.names=FALSE)
    p_value_list_signif <- p_value_list[p_value_list$p_value<0.05,]
    if (nrow(p_value_list_signif) != 0){
      #write.csv(p_value_list_signif, file = paste0("",sample_name,"_tp",tp,".csv"),row.names=FALSE)}
  }
}}

r<- read.csv("")
s <- read.csv("")
t<- read.csv("")
u<- read.csv("")

atlas_name <- "yeo"
atlas_overlap <- read.csv("")[,-1] #change atlas !!

for (sample in clinical_samples) {
  df <- df_features 
  df <- cbind(df,df_cov_te)
  df <- df[df$sample %in% sample,]
  tps <- unique(df$timepoint)
  i=0
  for (tp in tps) {
    # for all clinical samples, only look at baseline subjects
    # if you want to run it only on first timepoint subjects of all clinical samples uncomment the following line
    #if (tp==2 & length(sample)>1){break}
    if (length(sample)>1){sample_name="all_clinical_sites"}else{sample_name=sample}
    df_tp <- df[df$timepoint == tp, ]
    MS_raw <- df_tp[,features]
    networks <- atlas_overlap$class_num
    classes_total <- max(atlas_overlap$class_num)
    p_value_list <- list()
    t_statistic_list <- list()
    cohens_list <- list()
    
    # make scanner binary variables
    library(tidyverse)
    cov_te_scanners <- df_tp %>%
      mutate(value = 1) %>%
      pivot_wider(names_from=scanner, values_from=value, values_fill=0, names_prefix="scanner_")
    
    # Prepare covariates
    scanner_cols <- grep("^scanner", names(cov_te_scanners), value = TRUE)
    scanners <- cov_te_scanners[scanner_cols]
    colnames(scanners) <- paste0("scanner", 1:ncol(scanners))
    df_tp <- cbind(df_tp, scanners)
    df_tp <- subset(df_tp, select = -scanner) #remove scanner column
    
    columns_resid <- list()
    fitted_values_fit <- list()
    scanner_cols <- colnames(df_tp)[grepl("^scanner", colnames(df_tp))]
    formula <- paste("+ sex + euler_med +", paste(scanner_cols, collapse = " + "))
    for (class in 1:classes_total) {
      myregions <- which(networks == class)
      if (is_empty(myregions)) {
        print(paste0("Sample ", sample_name," timepoint ",tp," Class ",class," not present"))
        next}  # Skip to the next iteration
      
      nsubs <- nrow(df_tp)
      classreg <- numeric(nsubs)
      # Average features across all regions that belong to the same class, for each subject
      for (subj in 1:nsubs) {
        classreg[subj] <- mean(as.numeric(MS_raw[subj, myregions]))}
      
      reg <- cbind(df_tp[,63:74], classreg)
      cat("Class ",class)
      reg <- cbind(df_tp[,63:74], classreg)
      controls <- reg[reg$dcode==0,"classreg"]
      pats <- reg[reg$dcode==1,"classreg"]
      print("mean hc ")
      print(round(mean(controls),3))
      print("sd hc ")
      print(round(sd(controls),3))
      print("mean pats ")
      print(round(mean(pats),3))
      print("sd pats ")
      print(round(sd(pats),3))
      fit <- gam(as.formula(paste("classreg ~ s(age, bs = 'cr', k = 5)", formula)), data = reg)
      
      
      residuals<-fit$residuals
      fitted_values_f<-fit$fitted.values
      columns_resid<-residuals
      fitted_values_fit<-fitted_values_f
      
      #plot(df$age,fitted_values_fit[[32]])
      
      residuals<-as.data.frame(columns_resid)
      #names(residuals)<-paste0(features,"_resid")
      df_tp$dcode <- ifelse(df_tp$dcode == 0, 2, ifelse(df_tp$dcode == 1, 1, df_tp$dcode))
      dcode <- factor(df_tp$dcode, labels = c('Schizophrenia', 'Healthy'))
      residuals<- data.frame(cbind(residuals,dcode))
      residuals$dcode <- factor(residuals$dcode, labels = c('Schizophrenia', 'Healthy'))
      
      # using compare groups for t-test
      groups <- compareGroups(dcode ~ .,method = NA, data = residuals,alpha = 0.05)
      for(i in seq_along(groups)){
        p_value<-groups[[i]]$p.overall
        p_value_list[[i]]<-p_value
      }
      groups <- compareGroups(dcode ~ .,method = NA, data = residuals,alpha = 0.05)
      normality <- attributes(groups[[1]])$method[2] # normal or non-normal
      if (normality == "normal") {
        # Parametric test (t-test)
        param <- "p" #parametric test
        t_test_result <- t.test(residuals$columns_resid ~ residuals$dcode)
        cohens <- cohens_d(residuals$columns_resid ~ residuals$dcode)$Cohens_d
        p_value <- t_test_result$p.value
        t_statistic <- t_test_result$statistic
      } else {
        # Non-parametric test (Man-Whitney U test)
        param <- "np" #non-parametric test
        wilcox_test_result <- wilcox.test(residuals$columns_resid ~ residuals$dcode)
        cohens <- cohens_d(residuals$columns_resid ~ residuals$dcode)$Cohens_d
        p_value <- wilcox_test_result$p.value
        t_statistic <- wilcox_test_result$statistic
      }
      
      p_value_list[[class]] <- p_value
      t_statistic_list[[class]] <- t_statistic
      cohens_list[[class]] <- cohens
      
      if(class==1){class_name = "Visual"}
      if(class==2){class_name = "Somato Motor"}
      if(class==3){class_name = "Dorsal Attention"}
      if(class==4){class_name = "Ventral Attention"}
      if(class==5){class_name = "Limbic"}
      if(class==6){class_name = "Fronto Parietal"}
      if(class==7){class_name = "Default Mode"}
      
      dcode_column <- residuals$dcode
      dcode_column <- factor(dcode_column, labels = c('Schizophrenia', 'Healthy'))
      reg_column <- residuals$columns_resid
      b <- as.data.frame(cbind(dcode_column,reg_column))
      b$dcode_column <- factor(b$dcode_column, labels = c('Schizophrenia', 'Healthy'))
      colnames(b) <- c("dcode_column","network")
      if (tp==1){tp_text="Baseline"}else{tp_text="Follow-up"}
      plot <- ggbetweenstats(data=b, x=dcode_column, y=network,
                             ylab = "MS",
                             xlab="",
                             type=param, 
                             effsize.type="biased",
                             bf.message = FALSE,
                             centrality.label.args = list(size = 4, nudge_x = 0.4, segment.linetype = 4,min.segment.length = 0),
                             title=paste0(class_name," Network"," (",tp_text,")"),
                             ggplot.component = list(ggplot2::scale_y_continuous(
                               breaks = seq(-0.04, 0.04, 0.04),
                               limits = (c(-0.04, 0.04)))))+
        scale_color_manual(values = c("#1F78B4","#33A02C"))
                             #ggplot.component = list(theme(plot.subtitle = element_text(size = 24, face = "bold"))))
      sub_1 <- plot[["plot_env"]][["subtitle"]][[2]][3]
      sub_1 <- gsub("\\(|\\)", "", sub_1)
      sub_2 <- plot[["labels"]][["subtitle"]][[4]][3]
      sub_2 <- gsub("\\(|\\)", "", sub_2)
      sub_3 <- plot[["labels"]][["subtitle"]][[5]][3]
      sub_3 <- gsub("\\(|\\)", "", sub_3)
      sub_3 <- as.numeric(gsub("[^0-9.-]", "", sub_3))
      sub_4 <- plot[["labels"]][["subtitle"]][[6]][2]
      sub_4 <- gsub("\\(|\\)", "", sub_4)
      sub_5 <- plot[["plot_env"]][["subtitle"]][[2]][2]
      sub_5 <- as.numeric(gsub("[^0-9.]", "", sub_5))
      if (normality!="normal"){plot <- plot + 
        labs(subtitle = ggplot2::expr(paste(italic("W")["Mann-Whitney"]," = ", !!trunc(as.numeric(sub_1)) ,
                                            ", ", widehat(italic("r"))["biserial"]^"rank", " = ", !!sub_2 ,
                                            ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]"
        )))} else{plot <- plot + 
          labs(subtitle = ggplot2::expr(paste(italic("t")["Welch"],"(",!!sub_5,")"," = ", !!sub_1 ,
                                              ", ", widehat(italic("d"))["Cohen"], " = ", !!sub_2 ,
                                              ", CI"["95%"]," [", !!sub_3, ", ", !!sub_4, "]")))
        }
      plot <- plot + theme(axis.text = element_text(size = 15),  
                           axis.title = element_text(size = 16),  
                           plot.subtitle = element_text(size = 15, face="italic"),
                           plot.title = element_text(size = 18)) 
      png(file = paste0("",class,"_",sample_name,"_tp",tp,"_MS.png"), width=6, height=4.5, units="in", res=300)
      print(plot)
      dev.off()
    }
    fdr_p_value <- p.adjust(do.call(c, p_value_list), method = "fdr")
    p_value_list <- data.frame(network = 1:classes_total, p_value = do.call(c, p_value_list),statistic = do.call(c, t_statistic_list),fdr_p = fdr_p_value, cohens_d=do.call(c, cohens_list))
    write.csv(p_value_list, file = paste0("",sample_name,"_tp",tp,".csv"),row.names=FALSE)
    p_value_list_signif <- p_value_list[p_value_list$p_value<0.05,]
    if (nrow(p_value_list_signif) != 0){
      write.csv(p_value_list_signif, file = paste0("",sample_name,"_tp",tp,".csv"),row.names=FALSE)}
  }
}
i <- read.csv("")
j  <- read.csv("")
k <- read.csv("")
l <- read.csv("")



#---------------------------------------------------
# Z-scores for each region belonging to a network --
#---------------------------------------------------
yeo <- read.csv("")

#----------------------------------------------
# --- Distribution p-values 19 regions --------
#----------------------------------------------
library(dplyr)

# Assuming your dataframes are named df_Z_scores and df_cov_te

# Function to perform one iteration
perform_iteration <- function(df_Z_scores, df_cov_te) {
  # Select 19 random regions
  selected_regions <- sample(1:ncol(df_Z_scores), 19, replace = FALSE)
  
  # Filter Z scores dataframe for selected regions
  selected_zscores <- df_Z_scores[, selected_regions]
  
  # Subset Z scores for control and patient groups
  control_zscores <- selected_zscores[df_cov_te$dcode == 0, ]
  patient_zscores <- selected_zscores[df_cov_te$dcode == 1, ]
  
  # Calculate mean Z scores for control and patient subjects
  control_mean <- rowMeans(control_zscores)
  patient_mean <- rowMeans(patient_zscores)
  
  # Perform t-test
  t_test_result <- t.test(control_mean, patient_mean)
  
  # Extract p-value from t-test result
  p_value <- t_test_result$p.value
  t_value <- t_test_result[["statistic"]][["t"]]
  return(t_value)
}

# Perform 1000 iterations
set.seed(123)  # For reproducibility
t_values <- replicate(1000, perform_iteration(df_Z_scores, df_cov_te))

# Plot distribution of p-values
hist(t_values, breaks = 40,xlim = c(-10,10), col = "lightgreen", main = "Distribution of t-statistics", xlab = "t-statistic", ylab = "Frequency")


# Calculate real distribution of p-value between Hc and PATs from default mode ntw
yeo <- read.csv("")
regions_default_ntw <- yeo[yeo$network == 7, "region"]

# Filter Z scores dataframe for selected regions
selected_zscores <- df_Z_scores[, regions_default_ntw]

# Subset Z scores for control and patient groups
control_zscores <- selected_zscores[df_cov_te$dcode == 0, ]
patient_zscores <- selected_zscores[df_cov_te$dcode == 1, ]

# Calculate mean Z scores for control and patient subjects
control_mean <- rowMeans(control_zscores)
patient_mean <- rowMeans(patient_zscores)

# Perform t-test
t_test_result <- t.test(control_mean, patient_mean)

# Extract p-value from t-test result
p_value <- t_test_result$p.value
t_value <- t_test_result[["statistic"]][["t"]]

abline(v = t_value, col = "red", lwd = 2)


hist(p_values, breaks = 400, xlim = c(0, 0.05), col = "skyblue", main = "Distribution of p-values [0,0.05]", xlab = "p-value", ylab = "Frequency")
abline(v = p_value, col = "red", lwd = 2)


#--------------------------
# Table atlases Z-scores
#---------------------------
library(flextable)

economo_tp1 <- read.csv("")
economo_tp2 <- read.csv("")
economo_tp1_Z <- read.csv("")
economo_tp2_Z <- read.csv("")
economo_Z <- rbind(economo_tp1_Z,economo_tp2_Z)
economo <- rbind(economo_tp1,economo_tp2)
econ <- cbind(economo_Z, economo[,2:4])
econ["timepoint"] <- c(1,1,1,1,1,1,2,2,2,2,2,2)
#econ[,2:7] <- round(econ[,2:7],3)
econ[,2:7] <- as.numeric(sapply(econ[,2:7], function(x) sprintf("%.3f", round(x, 3))))

colnames(econ) <- c("Class","p-value","statistic","fdr-p","p-value ","statistic ","fdr-p ","Timepoint")
econ[,1] <- c("Agranular","Frontal","Parietal","Polar","Granulous","Transitional","Agranular","Frontal","Parietal","Polar","Granulous","Transitional")
econ <- econ[,c(8,1:7)]

econ[3, 3] <- paste0(format(econ[3, 3], nsmall = 3), " *")
econ[5,3] <- paste0(econ[5,3]," *")
econ[8,3] <- paste0(econ[8,3],0)
econ <- econ[,c(-5,-8)]

format_number <- function(x) {
  if (x < 10) {
    return(format(round(x, 3), nsmall = 3))
  } else {
    return(format(round(x), nsmall = 0))
  }
}
econ$statistic <- sapply(econ$statistic, format_number)
econ$`statistic `<- sapply(econ$`statistic `, format_number)

set_flextable_defaults(font.family = "Arial", font.size = 10,padding = 6,background.color = "white")
ft <- flextable(econ)
ft <- merge_v(ft, j = "Timepoint")
ft <- ft %>% 
  add_header_row(top = TRUE, values = c("Timepoint","Class", "Z-scores","meanMS"), colwidths = c(1,1,2,2))
ft
ft <- merge_v(ft, part = "header")
ft
ft <- theme_vanilla(ft)
ft
ft <- align(ft,align = c("center"),part = "all")
ft

save_as_image(ft, "")
save_as_docx(table, path="")

#--------- Yeo ------
economo_tp1 <- read.csv("")
economo_tp2 <- read.csv("")
economo_tp1_Z <- read.csv("")
economo_tp2_Z <- read.csv("")
economo_Z <- rbind(economo_tp1_Z,economo_tp2_Z)
economo <- rbind(economo_tp1,economo_tp2)
econ <- cbind(economo_Z, economo[,2:4])
econ["timepoint"] <- c(1,1,1,1,1,1,1,2,2,2,2,2,2,2)

colnames(econ) <- c("Network","p-value","statistic","fdr-p","p-value ","statistic ","fdr-p ","Timepoint")
econ[,1] <- c("Visual","Somato Motor","Dorsal Attention","Ventral Attention","Limbic","Fronto Parietal","Default Mode","Visual","Somato Motor","Dorsal Attention","Ventral Attention","Limbic","Fronto Parietal","Default Mode")

econ <- econ[,c(8,1:7)]

econ[,c(3:7)] <- as.numeric(sapply(econ[,c(3:7)], function(x) sprintf("%.3f", round(x, 3))))
#econ[,c(3:7)] <- round(econ[,c(3:7)],3)


econ <- econ[,c(-5,-8)]

format_number <- function(x) {
  if (x < 10) {
    return(format(round(x, 3), nsmall = 3))
  } else {
    return(format(round(x), nsmall = 0))
  }
}
econ$statistic <- sapply(econ$statistic, format_number)
econ$`statistic `<- sapply(econ$`statistic `, format_number)

econ[2,3] <- paste0(format(econ[2,3], nsmall = 3), " *")
econ[7,3] <- paste0(format(econ[7,3], nsmall = 3), " *")
econ[10,3] <- paste0(econ[10,3],0)
econ[14,3] <- paste0(econ[14,3], 0)

set_flextable_defaults(font.family = "Arial", font.size = 10,padding = 6,background.color = "white")
ft <- flextable(econ)
ft <- merge_v(ft, j = "Timepoint")
ft <- ft %>% 
  add_header_row(top = TRUE, values = c("Timepoint","Network", "Z-scores","meanMS"), colwidths = c(1,1,2,2))
ft
ft <- merge_v(ft, part = "header")
ft
ft <- theme_vanilla(ft)
ft
ft <- align(ft,align = c("center"),part = "all")
ft


save_as_image(ft, "")
save_as_docx(table, path="")

#-----------------------
# Table clustering Yeo
#----------------------
atlas <- read.csv("")[,-1]
clusters_both <- read.csv("")

atlas <- atlas[order(atlas$region), ]
clusters_both <- clusters_both[order(clusters_both$regions), ]

bind = cbind(atlas, clusters_both)
group1 = bind %>% filter(group==1)
table(group1$class_num)
group2 = bind %>% filter(group==2)
table(group2$class_num)


library(grid)
library(flextable)
cluster1_counts <- bind %>%
  filter(group == 1) %>%
  group_by(class_num) %>%
  summarize(Count = n())

# Filtrar las regiones del Cluster 2 y contarlas por clase
cluster2_counts <- bind %>%
  filter(group == 2) %>%
  group_by(class_num) %>%
  summarize(Count = n())

final_df <- merge(cluster1_counts, cluster2_counts, by = "class_num", all = TRUE)
final_df[is.na(final_df)] <- 0  # Reemplazar NA con 0 en caso de que no haya recuentos para alguna clase en algún cluster

# Renombrar
colnames(final_df) <- c("Network", "Cluster 1", "Cluster 2")
final_df[,1] <- c("Visual", "Somato Motor", "Dorsal Attention", "Ventral Attention", "Limbic", "Fronto Parietal", "Default Mode")

set_flextable_defaults(font.family = "Arial", font.size = 10,padding = 6,background.color = "white")
table <- flextable(final_df)
table <- theme_vanilla(table)
table

#flextable(head(final_df)) %>% bold(part = "header")
#final_df <- set_table_properties(final_df, align = "center", layout = "autofit")
#table <- border_inner_h(table, border = NULL, part = "all")

table <- border_inner_v(table, border = NULL, part = "all")
table
table <-border_outer(table, border = NULL, part = "all")
table <- align(table,align = c("center"),part = "all")
table
table <- set_table_properties(table, layout="autofit")
table <- autofit(table)
table
save_as_image(table, "",res=300)
save_as_docx(table, path="")

