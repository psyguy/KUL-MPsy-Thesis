library(grid)
library(gridBase)
library(ggplot2)

# start new page
plot.new() 

# setup layout
gl <- grid.layout(nrow=1, ncol=2)
# grid.show.layout(gl)

# setup viewports
vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1) 
vp.2 <- viewport(layout.pos.col=2, layout.pos.row=1) 
# init layout
pushViewport(viewport(layout=gl))
# access the first position
pushViewport(vp.1)

# start new base graphics in first viewport
par(new=TRUE, mar = rep(0.05, 4))# fig=gridFIG())

plot.new()
plot(x = 1:10, y = 10:1)

# done with the first viewport
popViewport()

# move to the next viewport
pushViewport(vp.2)

ggplotted <- qplot(x=1:10,y=10:1, 'point')
# print our ggplot graphics here
print(ggplotted, newpage = FALSE)

# done with this viewport
popViewport(1)



# trying pimage -----------------------------------------------------------




## put several pimages on a page (uses viewports and newpage = FALSE)
library(grid)
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow = 1, ncol = 2),
                      x=.5, y=.5,
                      width = 1.1, height = 1.1
                       ))

library(grid)
library(seriation)

set.seed(1)
x <- matrix(rnorm(90000), nrow = 300)

grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow = 1, ncol = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))

pimage(x, 
       axes = "both",
       key = FALSE,
       prop = TRUE,
       newpage = FALSE)

upViewport(1)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))


o <- seriate(x)
pimage(x, o,
       axes = "both",
       key = FALSE,
       prop = TRUE,
       newpage = FALSE)

upViewport(1)
popViewport(0)

hw <- 15
graph2pdf(height = hw, width=2*hw)

# igraph ------------------------------------------------------------------
library(grid)
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow = 2, ncol = 1)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))

extract_plotcon(m, save = T)

upViewport(1)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))

plot.new()

extract_plotnet(m, save = FALSE, first.add = TRUE)





# s -----------------------------------------------------------------------


