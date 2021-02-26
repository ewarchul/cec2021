deps <- c(
  "tidyverse",
  "magrittr",
  "devtools",
  "zeallot",
  "iterators",
  "doParallel",
  "foreach",
  "here",
  "pracma",
  "future",
  "furrr",
  "BBmisc"
)

diff = deps[
  !(deps %in% installed.packages()[, "Package"])
]

if (length(diff)) {
  install.packages(diff)
}

install.packages(
  "/deps/cecs_0.1.0.tar.gz",
  repos = NULL,
  type = "source"
)

install.packages(
  "/deps/cec2017_0.2.0.tar.gz",
  repos = NULL,
  type = "source"
)

install.packages(
  "/deps/cec2013_0.1-5.tar.gz",
  repos = NULL,
  type = "source"
)

install.packages(
  "/deps/cecb_0.0.0.9000.tar.gz",
  dependencies = TRUE,
  repos = NULL,
  type = "source"
)
