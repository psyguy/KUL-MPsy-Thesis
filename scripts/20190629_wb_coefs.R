name <- "Family 1 "
t3 <- Sys.time
print(name)
coefs.this.round <- coefs.wb.b %>%
  filter(`a Proportion` == "(1, 5, 0)" &
         `Epsilon Proportion` == "(0, 6, 0)")
owner.name <- coefs.this.round$Owner[1] %>% as.character()
title <- paste0("Evolution of coefficients for ",
                name,
                "eps ",
                coefs.this.round$`Epsilon Proportion`,
                ", a ",
                coefs.this.round$`a Proportion`)

p.cl <- ggplot(data = coefs.this.round,
             aes(x = Rewiring,
                 y = Clustering,
                 colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) +
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.cl) + 
  theme(legend.title = element_blank(),#element_text(size=12, color = "thistle4"),
        legend.justification=c(0,1), 
        legend.key.width = unit(1.5, "line"),
        legend.position= c(0, 1),#c(0.05, 0.95),
        legend.background = element_blank(),
        legend.text = element_text(size=8, color = "black", face="bold"),
        legend.key = element_blank())
  # + labs(subtitle="Legend: Top-Left Inside the Plot")

p.cl










p.sw <- ggplot(data = coefs.this.round,
             aes(x = Rewiring,
                 y = `Small World`,
                 colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) +
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.sw) + 
  theme(legend.title = element_blank(),#element_text(size=12, color = "thistle4"),
        legend.justification=c(0,1), 
        legend.key.width = unit(1.5, "line"),
        legend.position= c(0, 1),#c(0.05, 0.95),
        legend.background = element_blank(),
        legend.text = element_text(size=8, color = "black", face="bold"),
        legend.key = element_blank())


p.mo <- ggplot(data = coefs.this.round,
             aes(x = Rewiring,
                 y = Modularity,
                 colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) + 
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.mo) + 
  theme(legend.title = element_blank(),#element_text(size=12, color = "thistle4"),
        legend.justification=c(0,1), 
        legend.key.width = unit(1.5, "line"),
        legend.position= c(0, 1),#c(0.05, 0.95),
        legend.background = element_blank(),
        legend.text = element_text(size=8, color = "black", face="bold"),
        legend.key = element_blank())


p.ed <- ggplot(data = coefs.this.round,
             aes(x = Rewiring,
                 y = `Edge Density`,
                 colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) + 
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.de) + 
  theme(legend.title = element_blank(),#element_text(size=12, color = "thistle4"),
        legend.justification=c(0,1), 
        legend.key.width = unit(1.5, "line"),
        legend.position= c(0, 1),#c(0.05, 0.95),
        legend.background = element_blank(),
        legend.text = element_text(size=8, color = "black", face="bold"),
        legend.key = element_blank())


p.pl <- ggplot(data = coefs.this.round,
             aes(x = Rewiring,
                 y = `Avg Path Length`,
                 colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) + 
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.pl) + 
  theme(legend.title = element_blank(),#element_text(size=12, color = "thistle4"),
        legend.justification=c(0,1), 
        legend.key.width = unit(1.5, "line"),
        legend.position= c(0, 1),#c(0.05, 0.95),
        legend.background = element_blank(),
        legend.text = element_text(size=8, color = "black", face="bold"),
        legend.key = element_blank())


p.ef <- ggplot(data = coefs.this.round,
             aes(x = Rewiring,
                 y = Efficiency,
                 colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) + 
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.ef) + 
  theme(legend.title = element_blank(),#element_text(size=12, color = "thistle4"),
        legend.justification=c(0,1), 
        legend.key.width = unit(1.5, "line"),
        legend.position= c(0, 1),#c(0.05, 0.95),
        legend.background = element_blank(),
        legend.text = element_text(size=8, color = "black", face="bold"),
        legend.key = element_blank())


figure <- ggarrange(p.cl, p.ed, p.sw, p.pl, p.mo, p.ef,
                    ncol=2, nrow=3,
                    common.legend = F
                    # legend="bottom"
)

annotate_figure(figure, top = textGrob(title,
                                       gp=gpar(fontsize=35,font=8)
)
)# %>% print()
paste0("wb_", name, ".png") %>% ggsave(width = 14,
                                       height = 21,
                                       dpi = "retina")

Sys.time() - t3



# denom <- c(300, 50, 250)
# denom <- (denom*(denom-1)/2) %>% c(50*250)
# coefs.wb.b$Degree <- coefs.wb.b$Degree / rep(denom, nrow(coefs.wb.b)/4)

