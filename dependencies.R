deps <- c(
  "tidyverse",
  "magrittr",
  "devtools"
)

diff = deps[
  !(deps %in% installed.packages()[,"Package"])
]

if (length(diff)) {
  install.packages(diff)
}

install.packages(
  "cecs_0.1.0.tar.gz",
  repos = NULL,
  type = "source"
)

install.packages(
  "cec2017_0.2.0.tar.gz",
  repos = NULL,
  type = "source"
)

install.packages(
  "cec2013_0.1.5.tar.gz",
  repos = NULL,
  type = "source"
)

install.packages(
  "cec2013_0.0.0.9000.tar.gz",
  repos = NULL,
  type = "source"
)
