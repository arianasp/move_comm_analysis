#data frame of call types and their respective colors and symbols

types <- c(	'CC',		'LEAD CC',	'LEAD',		'MOVE',		'SOCIAL',	'SN',		'AGGRESS',	'ALARM',	'OTHER')
colors <- c('#CCCCCC',	'#006600' ,	'#00DD00',	'#00DD00',	'purple',	'orange',	'red',		'yellow',	'white')
symbs <-  c(1,			0,			0,			2,			12,			5,			4,			8,			12)	

call.types <- data.frame(type=types,col=colors,sym=symbs)

save(list=c('call.types'),file='~/Dropbox/meerkats/data/Kalahari2017/RData/call_colors.R')				