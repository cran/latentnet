#
pause <- function(){readline(prompt="Pause. Press <Enter> to continue...");invisible()}
#
# Use 'data(package = "statnet")' to list the data sets in it
#
data(package="latentnet")
#
# load the Sampson's Monks network
#
data(sampson)
pause()
#
# create a plot of the social network
#
plot(samplike)
pause()
#
# Set colors
#
#
# Create the plot symbols for the groups
#
oneL <- samplike %v% "group"
oneL[oneL=="Turks"] <- "T"
oneL[oneL=="outcasts"] <- "O"
oneL[oneL=="loyal"] <- "L"
oneL[c(1,7,15)] <- "W"
oneL
#
# Create colors
oneLcolors <- c("red","blue","black","green")[match(oneL,c("T","O","L","W"))]
#
# now make the vertex color reflect the groups
#
plot(samplike, vertex.col=oneLcolors, main="Liking Ties")
