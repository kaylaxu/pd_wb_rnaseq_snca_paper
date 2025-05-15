library(ggbeeswarm)
library(ggpubr)
library(patchwork)

mi_eval <- read.csv( "neutrophil_percent_prediction/mi_eval_split.csv")
bio_eval <- read.csv('neutrophil_percent_prediction/bio_eval_split.csv')
combined_eval <- read.csv( 'neutrophil_percent_prediction/combined_eval_split.csv')
xgb_eval <- read.csv('neutrophil_percent_prediction/xgb_eval_split.csv')
mi_eval$X <- NULL
bio_eval$X <- NULL
combined_eval$X <- NULL
xgb_eval$X <- NULL

mi_eval$Model <- "MI Linear"
bio_eval$Model <- "Cell-Based Linear"
combined_eval$Model <- "Combined Linear"
xgb_eval$Model <- "XGBoost"

eval_values <- rbind(mi_eval, bio_eval, combined_eval, xgb_eval)
eval_values$Model <- factor(eval_values$Model,levels=c('Cell-Based Linear', 'MI Linear', 'Combined Linear', 'XGBoost'))
                            
r2_plot <- ggplot(eval_values, aes(x=Model, y=R2, fill=Model)) +
  geom_violin(alpha=0.5) +
  geom_beeswarm(cex = 1.5) +
  geom_boxplot(width = 0.1, alpha=0.5) +
  ggtitle("R-squared in Test Sets ") +
  xlab('') +
  ylab('R-squared') + 
  theme_classic() +
  theme(text = element_text(size=20, color = 'black'),
        axis.text = element_text(size=20, color = 'black'),
        legend.text = element_text(size=20, color = 'black'),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_signif(comparisons = list(c("Cell-Based Linear", "MI Linear"),
                                 c("Cell-Based Linear", "Combined Linear"),
                                 c("MI Linear", "Combined Linear"),
                                 c("Cell-Based Linear", "XGBoost"),
                                 c("MI Linear", "XGBoost"),
                                 c("Combined Linear", "XGBoost")),
              test = "wilcox.test", step_increase = 0.075, textsize = 8,
              map_signif_level = TRUE, tip_length = 0) +
  guides(fill="none")


rmse_plot <- ggplot(eval_values, aes(x=Model, y=RMSE, fill=Model)) +
  geom_violin(alpha=0.5) +
  geom_beeswarm(cex = 1.5) +
  geom_boxplot(width = 0.1, alpha=0.5) +
  ggtitle("RMSE in Test Sets ") +
  xlab('') +
  ylab('RMSE') + 
  theme_classic() +
  theme(text = element_text(size=20, color = 'black'),
        axis.text = element_text(size=20, color = 'black'),
        legend.text = element_text(size=20, color = 'black'),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_signif(comparisons = list(c("Cell-Based Linear", "MI Linear"),
                                 c("Cell-Based Linear", "Combined Linear"),
                                 c("MI Linear", "Combined Linear"),
                                 c("Cell-Based Linear", "XGBoost"),
                                 c("MI Linear", "XGBoost"),
                                 c("Combined Linear", "XGBoost")),
              test = "wilcox.test", step_increase = 0.075, textsize = 8,
              map_signif_level = TRUE, tip_length = 0) +
  guides(fill="none")



mae_plot <- ggplot(eval_values, aes(x=Model, y=MAE, fill=Model)) +
  geom_violin(alpha=0.5) +
  geom_beeswarm(cex = 1.5) +
  geom_boxplot(width = 0.1, alpha=0.5) +
  ggtitle("MAE in Test Sets ") +
  xlab('') +
  ylab('MAE') + 
  theme_classic() +
  theme(text = element_text(size=20, color = 'black'),
        axis.text = element_text(size=20, color = 'black'),
        legend.text = element_text(size=20, color = 'black'),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_signif(comparisons = list(c("Cell-Based Linear", "MI Linear"),
                                 c("Cell-Based Linear", "Combined Linear"),
                                 c("MI Linear", "Combined Linear"),
                                 c("Cell-Based Linear", "XGBoost"),
                                 c("MI Linear", "XGBoost"),
                                 c("Combined Linear", "XGBoost")),
              test = "wilcox.test", step_increase = 0.075, textsize = 8,
              map_signif_level = TRUE, tip_length = 0) +
  scale_fill_discrete(name = "Model")




ggsave('neutrophil_percent_prediction/figure2.jpeg',plot = r2_plot + rmse_plot + mae_plot, dpi=300, width = 20, height=10)











