library(ggplot2)
library(viridis)
library(season)
library(gridExtra)
library(mapproj)
library(dichromat)
library(scales)
library(ggpubr)

# lie factor

data_df <- data.frame( size = c(rep("big", 100), rep("small", 150)))

p1 <- ggplot(data_df, aes(x = size)) + geom_bar() + theme_minimal()
p2 <- ggplot(data_df, aes(x = size)) + geom_bar() +
  scale_y_continuous(limits=c(100,150),oob = rescale_none) + theme_minimal()

p3 <- ggarrange(p1, p2) + theme_void()

ggplot2::ggsave("lie_factor.png", p3)


data_df <- data.frame(count = c(25, 75),
                      size = c("big", "small"))

p1 <- ggplot(data_df, aes(x = factor(1), y = count, fill = as.factor(size))) +
  geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") +theme_void()

p2 <- ggplot(data_df, aes(x = factor(1), y = count, fill = as.factor(size))) +
  geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y", start = -5) + theme_void()

p3 <- ggarrange(p1, p2) + theme_void()

ggplot2::ggsave("bad_pie.png", p3)

pa<-ggplot(schz, aes(year, month, fill = SczBroad)) + 
  geom_tile(colour="gray20", size=1.5, stat="identity") + 
  scale_fill_viridis(option="A") +
  scale_y_continuous(breaks=1:12, labels=month.abb[1:12])+
  xlab("") + 
  ylab("") +
  ggtitle("Total Australian Schizophrenics Born By Month and Year") +
  theme(
    plot.title = element_text(color="white",hjust=0,vjust=1, size=rel(2)),
    plot.background = element_rect(fill="gray20"),
    panel.background = element_rect(fill="gray20"),
    panel.border = element_rect(fill=NA,color="gray20", size=0.5, linetype="solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(), 
    axis.text = element_text(color="white", size=rel(1.5)),
    axis.text.y  = element_text(hjust=1),
    legend.text = element_text(color="white", size=rel(1.3)),
    legend.background = element_rect(fill="gray20"),
    legend.position = "bottom",
    legend.title=element_blank()
  )


pb<-ggplot(schz, aes(year, month, fill = SczBroad)) + 
  geom_tile(colour="gray20", size=1.5, stat="identity") + 
  scale_fill_viridis(option="B") +
  scale_y_continuous(breaks=1:12, labels=month.abb[1:12])+
  xlab("") + 
  ylab("") +
  ggtitle("Total Australian Schizophrenics Born By Month and Year") +
  theme(
    plot.title = element_text(color="white",hjust=0,vjust=1, size=rel(2)),
    plot.background = element_rect(fill="gray20"),
    panel.background = element_rect(fill="gray20"),
    panel.border = element_rect(fill=NA,color="gray20", size=0.5, linetype="solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(), 
    axis.text = element_text(color="white", size=rel(1.5)),
    axis.text.y  = element_text(hjust=1),
    legend.text = element_text(color="white", size=rel(1.3)),
    legend.background = element_rect(fill="gray20"),
    legend.position = "bottom",
    legend.title=element_blank()
  )


pc<-ggplot(schz, aes(year, month, fill = SczBroad)) + 
  geom_tile(colour="gray20", size=1.5, stat="identity") + 
  scale_fill_viridis(option="C") +
  scale_y_continuous(breaks=1:12, labels=month.abb[1:12])+
  xlab("") + 
  ylab("") +
  ggtitle("Total Australian Schizophrenics Born By Month and Year") +
  theme(
    plot.title = element_text(color="white",hjust=0,vjust=1, size=rel(2)),
    plot.background = element_rect(fill="gray20"),
    panel.background = element_rect(fill="gray20"),
    panel.border = element_rect(fill=NA,color="gray20", size=0.5, linetype="solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(), 
    axis.text = element_text(color="white", size=rel(1.5)),
    axis.text.y  = element_text(hjust=1),
    legend.text = element_text(color="white", size=rel(1.3)),
    legend.background = element_rect(fill="gray20"),
    legend.position = "bottom",
    legend.title=element_blank()
  )


pd<-ggplot(schz, aes(year, month, fill = SczBroad)) + 
  geom_tile(colour="gray20", size=1.5, stat="identity") + 
  scale_fill_viridis(option="D") +
  scale_y_continuous(breaks=1:12, labels=month.abb[1:12])+
  xlab("") + 
  ylab("") +
  ggtitle("Total Australian Schizophrenics Born By Month and Year") +
  theme(
    plot.title = element_text(color="white",hjust=0,vjust=1, size=rel(2)),
    plot.background = element_rect(fill="gray20"),
    panel.background = element_rect(fill="gray20"),
    panel.border = element_rect(fill=NA,color="gray20", size=0.5, linetype="solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(), 
    axis.text = element_text(color="white", size=rel(1.5)),
    axis.text.y  = element_text(hjust=1),
    legend.text = element_text(color="white", size=rel(1.3)),
    legend.background = element_rect(fill="gray20"),
    legend.position = "bottom",
    legend.title=element_blank()
  )




data(unemp, package = "viridis")

county_df <- map_data("county", projection = "albers", parameters = c(39, 45))
names(county_df) <- c("long", "lat", "group", "order", "state_name", "county")
county_df$state <- state.abb[match(county_df$state_name, tolower(state.name))]
county_df$state_name <- NULL

state_df <- map_data("state", projection = "albers", parameters = c(39, 45))

choropleth <- merge(county_df, unemp, by = c("state", "county"))
choropleth <- choropleth[order(choropleth$order), ]

ggplot(choropleth, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = rate), colour = alpha("white", 1 / 2), linewidth = 0.2) +
  geom_polygon(data = state_df, colour = "white", fill = NA) +
  coord_fixed() +
  theme_minimal() +
  ggtitle("US unemployment rate by county") +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) +
  scale_fill_viridis(option="magma")



pal <- function(col, ...)
  image(seq_along(col), 1, matrix(seq_along(col), ncol = 1),
        col = col, axes = FALSE, ...)
opar <- par(mar = c(1, 2, 1, 1))
layout(matrix(1:6, ncol = 1))
pal(colorschemes$BrowntoBlue.10, main = "Brown to Blue (10)")
pal(colorRampPalette(colorschemes$BrowntoBlue.10, space = "Lab")(100),
    main = "Brown to Blue Ramp")
pal(dichromat(colorschemes$BrowntoBlue.10),
    main = "Brown to Blue (10) -- deuteranopia")
pal(colorschemes$Categorical.12, main = "Categorical (12)")
pal(dichromat(colorschemes$Categorical.12),
    main = "Categorical (12) -- deuteranopia")
pal(dichromat(colorschemes$Categorical.12, "protan"),
    main = "Categorical (12) -- protanopia")
par(opar)


x <- y <- seq(-8*pi, 8*pi, len = 40)
r <- sqrt(outer(x^2, y^2, "+"))
filled.contour(cos(r^2)*exp(-r/(2*pi)), 
               axes=FALSE,
               color.palette=viridis,
               asp=1)
