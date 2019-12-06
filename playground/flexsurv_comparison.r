require(flexsurv)
require(elfDistr)
require(microbenchmark)

flexsurv.test = function(n = 10000){
	x = runif(0, 10, n=n);
	a = rep(1, 13);
	b = rep(1, 64);
	c = rep(1, 37);
	return( sum(dgengamma(x, a, b, c, log=T)) );
}

flexsurv.test.orig = function(n = 10000){
	x = runif(0, 10, n=n);
	a = rep(1, 13);
	b = rep(1, 64);
	c = rep(1, 37);
	return( sum(dgengamma.orig(x, a, b, c, log=T)) );
}

elf.test = function(n = 10000){
	x = runif(0, 10, n=n);
	a = rep(1, 13);
	b = rep(1, 64);
	c = rep(1, 37);
	return( sum(dggamma(x, a, b, c, log=T)) );
}

result = microbenchmark(
	flexsurv.test(),
	flexsurv.test.orig(),
	elf.test(),
	unit="ms",
	times=1000,
	control=list(warmup=100)
)

print(summary(result))


# Now we test if they give the same values
check = function(a, b){
	if(all(abs(a - b) < 1e-5)){
		print("OK!");
	} else {
		print("CRAP!");
	};
}
x = seq(1e-3, 10, length=100);
check(dggamma(x, 1, 1, 1), dgengamma.orig(x, 1, 1, 1))
check(dggamma(x, 0.5, 1, 5), dgengamma.orig(x, 0.5, 1, 5))
check(dggamma(x, 1, 10, 1), dgengamma.orig(x, 1, 10, 1))
check(dggamma(x, 3, 3, 1), dgengamma.orig(x, 3, 3, 1))
check(dggamma(x, 5, 0.3, 1), dgengamma.orig(x, 5, 0.3, 1))

df = NULL
seqN = seq(64, 10000, length=100)
for(n in seqN){
	cat("n = ", n, "\n");
	result = microbenchmark(
		flexsurv.test(n),
		elf.test(n),
		unit="ms",
		times=400,
		control=list(warmup=10)
	);
	
	df = rbind(df, summary(result)$mean);
}

print(df)
