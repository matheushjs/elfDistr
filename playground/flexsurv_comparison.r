require(flexsurv)
require(elfDistr)
require(microbenchmark)

flexsurv.test = function(){
	x = runif(0, 10, n=10000);
	return( sum(dgengamma(x, 1, 1, 1, log=T)) );
}

flexsurv.test.orig = function(){
	x = runif(0, 10, n=10000);
	return( sum(dgengamma.orig(x, 1, 1, 1, log=T)) );
}

elf.test = function(){
	x = runif(0, 10, n=10000);
	return( sum(dggamma(x, 1, 1, 1, log=T)) );
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
