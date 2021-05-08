library(ggplot2)

google.red <- '#EA4335'
google.yellow <- '#FBBC05'
google.green <- '#34A853'
google.blue <- '#4285F4'

google.colors <- c(google.blue, google.red, google.yellow, google.green)

default.colors <- c('#3266cc','#dc3812','#fe9900','#109619','#990099',
                    '#0099c5','#dd4578','#66aa00','#b82e2e','#316394',
                    '#994499','#21aa98','#aaab12','#6633cc','#e67300',
                    '#329262','#5474a5','#3c3ead','#8b0607','#641066',
                    '#f8756b','#e76af2','#02b0f7','#02bf7d','#8c564a',
                    '#e377c2','#9467bc','#7f7f7f','#bcbd23','#17bed0',
                    '#aec6e8','#ffbc78','#97df89','#ff9897','#c4b0d5',
                    '#c49c94','#f7b7d2','#dadb8d','#9edae5','#ffed6f')

test <- data.frame(x=1:40, y=1:40, stringsAsFactors = F)
ggplot() + geom_point(data = test, aes(x, y), color=default.colors, size=5)

colors <- c(google.colors, default.colors[5:40])
ggplot() + geom_point(data = test, aes(x, y), color=colors, size=5)

