plot_QC_paper <- function(seu_object, feature, ytext, xtext, paneltitle, colorviolin,...){
  VlnPlot(seu_object, pt.size = 0.01, features = feature, ncol = 1, group.by = "orig.ident", cols =rep(colorviolin, 10)) +
    stat_summary(fun.y = median, geom='point', size = 3, colour = "red", shape = 95)+
    labs(y = ytext, title = paneltitle, x=xtext) +
    #geom_hline(yintercept = 4000, size = 0.3) +
    #scale_x_discrete(labels = c("0", "2", "4", "6", "6 \n+ Ibru", "60", "180", "180 \n+ Ibru"))+
    theme_half_open()+
    theme(axis.text.x = element_text(colour = 'black',angle = 0, hjust=0.4),
          axis.text.y = element_text(colour = 'black'),
          # axis.title.x = element_blank(),
          plot.title = element_text(size=7),
          legend.position = "none",
          text = element_text(size=7),
          axis.text=element_text(size=7),
          axis.line = element_line(colour = 'black', size = 0.3),
          axis.ticks = element_line(colour = 'black', size = 0.3))
}
