time1 <- strptime('15 16:56:38 2020', "%d %H:%M:%S %Y")
time2 <- strptime('15 17:04:57 2020', "%d %H:%M:%S %Y")

mins <- as.integer(as.numeric(difftime(time2, time1, units = 'secs'))/60)
secs <- as.numeric(difftime(time2, time1, units = 'secs')) %% 60

message(paste0(mins, "'", secs, '"'))
