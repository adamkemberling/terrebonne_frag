#using plotly

library(plotly)

#basic histogram
p <- plot_ly(x = ~rnorm(50), type = "histogram")
p

#normalized histogram
p <- plot_ly(x = ~rnorm(50),
             type = "histogram",
             histnorm = "probability")
p


#create a shareable link to your chart by setting up plotly account
#set up API credentials: https://plot.ly/r/getting-started

#using plotly with ggplot
set.seed(100)
d <- diamonds[sample(nrow(diamonds), 1000), ]
p <- ggplot(data = d, aes(x = carat, y = price)) +
  geom_point(aes(text = paste("Clarity:", clarity)), size = 4) +
  geom_smooth(aes(colour = cut, fill = cut)) + 
  facet_wrap(~ cut)

(gg <- ggplotly(p))
