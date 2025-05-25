library(ggplot2)
doubletdecon_bench<-read.csv("./datasets/result/benchmark_doubletdecon.csv")
doubletfinder_bench<-read.csv("./datasets/result/doubletfinder_benchmark.csv")
scrublet_bench<-read.csv("./datasets/result/scrublet_benchmark.csv")
doubletdecon_bench<-doubletdecon_bench[,-1]
doubletdecon_bench$type<-"DoubleDecon"
doubletfinder_bench$type<-"DoubletFinder"
scrublet_bench$type<-"Scrublet"
data<-rbind(doubletdecon_bench,doubletfinder_bench)
data<-rbind(data,scrublet_bench)
# Specificity
data %>%
  ggplot( aes(x=type, y=specificity, fill=type)) +
  geom_boxplot(width=0.5,outlier.shape = NA) +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=13),
    panel.grid=element_blank()
  ) +
  ggtitle("Specificty Benchmark") +
  xlab("")

# Precision
data %>%
  ggplot( aes(x=type, y=precision, fill=type)) +
  geom_boxplot(width=0.5,outlier.shape = NA) +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=13),
    panel.grid=element_blank()
  ) +
  ggtitle("Precision Benchmark") +
  xlab("")


# F1 Score
data %>%
  ggplot( aes(x=type, y=f1_score, fill=type)) +
  geom_boxplot(width=0.5,outlier.shape = NA) +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=13),
    panel.grid=element_blank()
  ) +
  ggtitle("F1 score Benchmark") +
  xlab("")

# Recall
data %>%
  ggplot( aes(x=type, y=recall, fill=type)) +
  geom_boxplot(width=0.5,outlier.shape = NA) +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=13),
    panel.grid=element_blank()
  ) +
  ggtitle("Recall Benchmark") +
  xlab("")

# Comparison between number and expected rate
doubletfinder_th<-read.table('./datasets/result/doubletfinder_threshold.txt')
scrublet_th<-read.table('./datasets/result/srublet_threshold.txt')
doubletfinder_iderate<-read.table('./datasets/result/srublet_benchmark_threshold.txt')
scrublet_iderate<-read.table('./datasets/result/srublet_benchmark_threshold.txt')
colnames(doubletfinder_iderate)<-doubletfinder_iderate[1,]
colnames(scrublet_iderate)<-scrublet_iderate[1,]
doubletfinder_iderate<-doubletfinder_iderate[-1,]         
scrublet_iderate<-scrublet_iderate[-1,]
doubletfinder_iderate_0.1<-doubletfinder_iderate[,c(1:4)]
colnames(doubletfinder_iderate_0.1)<-c("precision","recall","specificity","f1_score")
doubletfinder_iderate_0.1$type<-"doubletfinder_0.1"
doubletfinder_iderate_0.2<-doubletfinder_iderate[,c(5:8)]
colnames(doubletfinder_iderate_0.2)<-c("precision","recall","specificity","f1_score")
doubletfinder_iderate_0.2$type<-"doubletfinder_0.2"
doubletfinder_iderate_0.4<-doubletfinder_iderate[,c(9:12)]
colnames(doubletfinder_iderate_0.4)<-c("precision","recall","specificity","f1_score")
doubletfinder_iderate_0.4$type<-"doubletfinder_0.4"


scrublet_iderate_0.1<-scrublet_iderate[,c(1:4)]
colnames(scrublet_iderate_0.1)<-c("precision","recall","specificity","f1_score")
scrublet_iderate_0.1$type<-"scrublet_0.1"

scrublet_iderate_0.2<-scrublet_iderate[,c(5:8)]
colnames(scrublet_iderate_0.2)<-c("precision","recall","specificity","f1_score")
scrublet_iderate_0.2$type<-"scrublet_0.2"

scrublet_iderate_0.4<-scrublet_iderate[,c(9:12)]
colnames(scrublet_iderate_0.4)<-c("precision","recall","specificity","f1_score")
scrublet_iderate_0.4$type<-"scrublet_0.4"

# Precision
plot_data<-rbind(scrublet_iderate_0.1,scrublet_iderate_0.2)
plot_data<-rbind(plot_data,scrublet_iderate_0.4)
plot_data<-rbind(plot_data,scrublet_bench[,-1])
plot_data$precision<-as.numeric(plot_data$precision)
plot_data %>%
  ggplot( aes(x=type, y=precision, fill=type)) +
  geom_boxplot(width=0.5,outlier.shape = NA) +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    panel.grid=element_blank()
  ) +
  ggtitle("Scrublet precision under different conditions") +
  xlab("")


# F1 Score
plot_data$f1_score<-as.numeric(plot_data$f1_score)
plot_data %>%
  ggplot( aes(x=type, y=f1_score, fill=type)) +
  geom_boxplot(width=0.5,outlier.shape = NA) +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    panel.grid=element_blank()
  ) +
  ggtitle("Scrublet F1 score under different conditions") +
  xlab("")

# Recall
plot_data<-rbind(doubletfinder_iderate_0.1,doubletfinder_iderate_0.2)
plot_data<-rbind(plot_data,doubletfinder_iderate_0.4)
plot_data<-rbind(plot_data,doubletfinder_bench[,-1])
plot_data$recall<-as.numeric(plot_data$recall)
plot_data %>%
  ggplot( aes(x=type, y=recall, fill=type)) +
  geom_boxplot(width=0.5,outlier.shape = NA) +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    panel.grid=element_blank()
  ) +
  ggtitle("DoubletFinder recall under different conditions") +
  xlab("")
