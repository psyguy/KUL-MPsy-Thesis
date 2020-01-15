name <- "Logan Wauters"#title# <- coefs.wb.b$Owner %>% as.character() %>% unique()
t3 <- Sys.time()
print(name)
coefs.this.round <- coefs.wb.b %>% filter(Owner == name)
owner.name <- coefs.this.round$Owner[1] %>% as.character()
title <- paste0("Evolution of coefficients for ",
                name,
                "eps ",
                coefs.this.round$`Epsilon Proportion`,
                ", a ",
                coefs.this.round$`a Proportion`)

p1 <- ggplot(data = coefs.this.round,
             aes(x = Rewiring,
                 y = Clustering,
                 colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) +
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.cl)

p2 <- ggplot(data = coefs.this.round,
             aes(x = Rewiring,
                 y = `Small World`,
                 colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) +
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.sw)

p3 <- ggplot(data = coefs.this.round,
             aes(x = Rewiring,
                 y = Modularity,
                 colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) + 
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.mo)

p4 <- ggplot(data = coefs.this.round,
             aes(x = Rewiring,
                 y = `Edge Density`,
                 colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) + 
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.de)

p5 <- ggplot(data = coefs.this.round,
             aes(x = Rewiring,
                 y = `Avg Path Length`,
                 colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) + 
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.pl)

p6 <- ggplot(data = coefs.this.round,
             aes(x = Rewiring,
                 y = Efficiency,
                 colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) + 
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.ef)


figure <- ggarrange(p1, p2, p3, p4, p5, p6,
                    ncol=1, #nrow=2,
                    common.legend = TRUE,
                    legend="bottom"
)

annotate_figure(figure, top = title)# %>% print()
paste0("wb_", name, ".png") %>% ggsave(width = 7,
                                       height = 42,
                                       dpi = "retina")

Sys.time() - t3
