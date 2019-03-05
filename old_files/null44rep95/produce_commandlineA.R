# to produce commandline arguments
for(i in 1:50) {
	tmp <- sample(1:44,44)
	cat(c(formatC(i, width = 2, format = "d", flag = "0"),tmp,'\n'))
}
