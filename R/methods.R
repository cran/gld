print.starship <- function(x,digits = max(3, getOption("digits") - 3), ...)
{
# I should include the call
# Add names to the vector in starship
cat(paste(x$param,"parameterisation\n"))
print.default(format(x$lambda,digits=digits), print.gap = 2,quote=FALSE)
}

summary.starship <- function(object,...)
{
cat("\nAdaptive Grid estimates:\n")
fake.lambda <- object$grid.results$lambda
names(fake.lambda) <- paste("lambda",1:length(fake.lambda),sep="")
fake.starship.object <- list(param=object$param,lambda=fake.lambda)
print.starship(fake.starship.object)
cat(paste("internal g-o-f measure at grid minimum:",
format(object$grid.results$response),"\n"))
cat("\nOptim (final) estimates (starting from grid estimates):\n")
print.starship(object)
cat(paste("internal g-o-f measure at optim minimum:",
format(object$optim.results$value),"\n"))
cat("optim.details:\nCounts: ")
print(object$optim.results$counts)
cat("Convergence: ")
print(object$optim.results$convergence)
cat("Message: ")
print(object$optim.results$message)
}

plot.starship <- function(x,ask=FALSE,one.page=TRUE,breaks="Sturges",histogram.title=NULL,...)
{
if (ask) {par(ask) <- TRUE}
if (one.page) {opar <- par(mfrow=c(2,1))}
qqgl(y=x$data,lambda.pars1=x$lambda,param=x$param,xlab="Fitted Theoretical Quantiles")
hist(x$data,prob=TRUE,xlab="Data",breaks=breaks,main=histogram.title,...)
plotgld(lambda1=x$lambda,param=x$param,new=FALSE,...)
if (one.page) {par(opar)} # Return to previous par
}
