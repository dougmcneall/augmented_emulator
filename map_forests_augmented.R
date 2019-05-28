# map_forests.R
# Map the forest regions used in the augmented emulator paper.
# 

pdf(height=4, width=8,file='graphics/map_forests_augmented.pdf')
#dev.new(height = 4, width = 8)
map('world', fill = TRUE,  ylim = c(-40,40), col = 'darkgrey',
    mar = c(0.1, 0.1, 0.1, 0.1), border = 'darkgrey')

polygon(x = c(-90, -45, -45, -90), y = c(-15, -15, 15, 15)) # Amazon
text(-94,11, 'Amazon', cex = 0.6, pos = 4, font = 2)

polygon(x = c(7.5, 30, 30, 7.5), y = c(-15,-15,10,10))
text(3.5,2, 'Central\nAfrica', cex = 0.6, pos = 4, font = 2)

polygon(x = c(90,150,150, 90), y = c(-12,-12,10,10))
text(86,6, 'SE Asia', cex = 0.6, pos = 4, font = 2)

dev.off()