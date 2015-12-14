# Reference:
# @ARTICLE{bootsvd,
#   author = {{Fisher}, A. and {Caffo}, B. and {Schwartz}, B. and {Zipunnikov}, V.},
#   title = "{Fast, Exact Bootstrap Principal Component Analysis for p$\gt$1 million}",
#   journal = {ArXiv e-prints},
#   year = 2014,
#   month = may
# }

# This function accepts truncated svds as well. 
getBootSVDs <- function(u,v,d, B=100){
	n <- nrow(v)
	svds <- pblapply(1:B, function(b){
		# Boot indices:
		Pb <- matrix(0, n,n)
		Pb[cbind(sample(n, replace=TRUE), 1:n)] <- 1

		# Calculating bootstrapped svd:
		DVPb <- diag(d) %*% t(v) %*% Pb
		Psvd <- svd(DVPb)
		ub <- u %*% ab
		vb <- Psvd$v

		# Taking care of the arbitrary signs:
		ab <- Psvd$u
		signs <- sign(diag(ab))
		ub <- t(t(ub) * signs)
		vb <- vb * signs

		list(u=ub, v=vb, d=db)
	})
	svds
}
