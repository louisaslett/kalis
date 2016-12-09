MakeForwardTable <- function(from_recipient, to_recipient, num_donors) {
  if(from_recipient>to_recipient) {
    stop("from_recipient must be smaller than to_recipient.")
  }
  if(from_recipient < 1) {
    from_recipient <- 1
  }
  if(to_recipient > N) {
    to_recipient <- N
  }
  delN <- to_recipient-from_recipient+1

  list(alpha          = matrix(0, num_donors, delN),
       alpha.f        = rep(0, delN),
       alpha.f2       = rep(0, delN),
       l              = c(0),
       from_recipient = from_recipient,
       to_recipient   = to_recipient)
}

MakeBackwardTable <- function(from_recipient, to_recipient, num_donors) {
  if(from_recipient>to_recipient) {
    stop("from_recipient must be smaller than to_recipient.")
  }
  if(from_recipient < 1) {
    from_recipient <- 1
  }
  if(to_recipient > N) {
    to_recipient <- N
  }
  delN <- to_recipient-from_recipient+1

  list(beta           = matrix(0, num_donors, delN),
       beta.g         = rep(0, delN),
       beta.g2        = rep(0, delN),
       l              = c(2147483647),
       from_recipient = from_recipient,
       to_recipient   = to_recipient)
}
