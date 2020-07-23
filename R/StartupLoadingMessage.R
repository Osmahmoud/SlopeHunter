SlopeHunterStartupMessage <- function()
{
  # Startup message can be obtained as
  # > figlet -f slant SlopeHunter
  msg <- c(paste(cat(
    "
    ______                        __   __
   / ____/__    ____  _________  / /  / /_  ____   __________________
  / /__  / /   / _  \\/ _  / __/ / /__/ / / / / \\  / /_  __/ ___/ _  /
  \\___ \\/ /   / / / / ___/ /_  /  __  / / / / /\ \\/ / / / / /_ /   _/
 ____/ / /___/ /_/ / /  / /_  / /  / / /_/ / /  \\\ / / / / /__/ /\ \\
/_____/\\____/\\____/_/  /___/ /_/  /_/\\____/_/  \\_/ /_/ /____/_/ \\_\\  version"), packageVersion("SlopeHunter"),

    "\nType 'citation(\"SlopeHunter\")' for citing this R package in publications."
))
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # unlock .SlopeHunter variable allowing its modification
  # unlockBinding(".SlopeHunter", asNamespace("SlopeHunter"))
  # startup message
  msg <- SlopeHunterStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'SlopeHunter' version", packageVersion("SlopeHunter"))
  packageStartupMessage(msg)
  invisible()
}
