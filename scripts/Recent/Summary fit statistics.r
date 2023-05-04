#author:
#Francis Durnin-Vermette (francis.durnin-vermette@agr.gc.ca)

validation_calculate_stats <- function(simulated, measured) {
	# Takes: Two vectors
	# Returns: Tibble of all validation stats including root mean square error (RMSE),
	# 				 normalized RMSE (NRMSE), and Nash-Sutcliffe Efficiency (NSE)
	
	# Root mean square error (RMSE)
	P <- simulated %>%
		na.omit()
	O <- measured %>%
		na.omit()
	PminusO <- P - O
	PminusOsquared <- PminusO^2
	PminusOsquaredsums <- sum(PminusOsquared)
	n <- length(PminusO)
	RMSE <- sqrt((PminusOsquaredsums)/n)
	
	# Normalized RMSE (NRMSE)
	Obar <- mean(O)
	NRMSE <- RMSE / Obar
	
	# Nash-Sutcliffe Efficiency (NSE)
	OiminusObar <- O - Obar
	OiminusObarsquared <- OiminusObar^2
	OiminusObarsquaredsums <- sum(OiminusObarsquared, na.rm = TRUE)
	NSE <- ((OiminusObarsquaredsums) - (PminusOsquaredsums))/OiminusObarsquaredsums
	
	# Relative mean percent error (PE)
	OminusPoverO <- (O - P)/O
	OminusPoverOsums <- sum(OminusPoverO, na.rn=TRUE)
	PE <- 100 / (n * OminusPoverOsums)
	
	# d-index
	Pbar <- mean(P)
	PiminusPbar <- P - Pbar
	PiminusPbarsquared <- PiminusPbar^2
	PiminusPbarsquaredsums <- sum(PiminusPbarsquared, na.rm = TRUE)
	d_indexdenominator <- sum((abs(PiminusPbar) + abs(OiminusObar))^2, na.rm = TRUE)
	d_index = 1 - (PiminusPbarsquaredsums / d_indexdenominator)
	
	# R2
	lm_result <- summary(lm(O~P))
	r2 <- lm_result$r.squared
	r2_adj <- lm_result$adj.r.squared
	
	# Return all
	all_stats <- c(rmse = RMSE,
								 nrmse = NRMSE,
								 nse = NSE,
								 pe = PE,
								 d_index = d_index,
								 r2 = r2,
								 r2_adj = r2_adj)
	
	return(all_stats)
}
