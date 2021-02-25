library(tidyverse)
library(magrittr)

get_aocs = function(cec, dim, prob, rep) {
interp_ = stringr::str_interp
 if (cec %in% c(13, 14, 17)) {
    paths = c(
      interp_("~/cec2020/data/final/cec${cec}/csa-cec${cec}-best"),
      interp_("~/cec2020/data/final/cec${cec}/csa-cec${cec}-mean"),
      interp_("~/cec2020/data/final/cec${cec}/msr-cec${cec}-best"),
      interp_("~/cec2020/data/final/cec${cec}/msr-cec${cec}-mean"),
      interp_("~/cec2020/data/final/cec${cec}/ppmf-cec${cec}-ave")
    )
  } else {
    paths = c(
      interp_("~/cec2020/data/final/cec${cec}/csa-bsr-best"),
      interp_("~/cec2020/data/final/cec${cec}/csa-bsr-mean"),
      interp_("~/cec2020/data/final/cec${cec}/msr-bsr-best"),
      interp_("~/cec2020/data/final/cec${cec}/msr-bsr-mean"),
      interp_("~/cec2020/data/final/cec${cec}/ppmf-bsr-ave")
    )
  }
  df = cecb::get_dfr(
    paths,
    list(dim = dim, prob = prob, rep = rep)
  )
  df %>%
    set_paper_names(cec) %>%
    cecb::compute_aoc()
}


set_paper_names = function(dfx, cec) {
  interp_ = stringr::str_interp
  if (cec != 21) {
    dfx %>%
      dplyr::mutate(
        Method = 
          ifelse(Method == interp_("csa-cec${cec}-best"), "IPOP-CSA",
                 ifelse(Method == interp_("csa-cec${cec}-mean"), "RB-IPOP-CSA",
                        ifelse(Method == interp_("msr-cec${cec}-best"), "IPOP-MSR",
                               ifelse(Method == interp_("msr-cec${cec}-mean"), "RB-IPOP-MSR",
                                      ifelse(Method == interp_("ppmf-cec${cec}-ave"), "RB-IPOP-PPMF", "X")))))
      )
  } else {
  dfx %>%
    dplyr::mutate(
      Method = 
        ifelse(Method == interp_("csa-bsr-best"), "IPOP-CSA",
               ifelse(Method == interp_("csa-bsr-mean"), "RB-IPOP-CSA",
                      ifelse(Method == interp_("msr-bsr-best"), "IPOP-MSR",
                             ifelse(Method == interp_("msr-bsr-mean"), "RB-IPOP-MSR",
                                    ifelse(Method == interp_("ppmf-bsr-ave"), "RB-IPOP-PPMF", "X")))))
    )
  }

}

get_cec_plot = function(cec, dim, prob, rep) {
  interp_ = stringr::str_interp
  if (cec %in% c(13, 14, 17)) {
    paths = c(
      interp_("~/cec2020/data/cec${cec}/csa-cec${cec}-best"),
      interp_("~/cec2020/data/cec${cec}/csa-cec${cec}-mean"),
      interp_("~/cec2020/data/cec${cec}/msr-cec${cec}-best"),
      interp_("~/cec2020/data/cec${cec}/msr-cec${cec}-mean"),
      interp_("~/cec2020/data/cec${cec}/ppmf-cec${cec}-ave")
    )
  } else {
    paths = c(
      interp_("~/cec2020/data/cec${cec}/csa-bsr-best"),
      interp_("~/cec2020/data/cec${cec}/csa-bsr-mean"),
      interp_("~/cec2020/data/cec${cec}/msr-bsr-best"),
      interp_("~/cec2020/data/cec${cec}/msr-bsr-mean"),
      interp_("~/cec2020/data/cec${cec}/ppmf-bsr-ave")
    )
  }
  df = 
    cecb::get_dfr(
      paths,
      list(dim = dim, problems = 1:prob, rep = rep)
    ) %>%
  set_paper_names(cec)
  df %>%
    cecb::ecdf_plot() + 
    ggplot2::theme(
      legend.position = c(0.8, 0.15),
      legend.key = 
        ggplot2::element_rect(colour = NA, fill = NA)
    )
}

#' Generate table with benchmark results
#' 
#' @param idpath path to benchmark results files :: String
#' @param problems list of function indices :: [Int] 
#' @param dim dimensionality
#' @param type type of budget step counter: {M, m} :: String
#' @param ... extra params for `kbl()` 
#' @export

get_resultTable = function(idpath, problems, dim, type = "M", ...) {
    format_ = purrr::partial(base::formatC, format = "e", digits = 4)
    table = problems %>%
     purrr::map(function(p){
         datapath = stringr::str_glue("{idpath}/{type}/{type}-{p}-D-{dim}.txt")
         df = readr::read_delim(datapath, ",", col_names = FALSE, col_types = cols(.default = col_double()))
         df %>%
            dplyr::last() %>%
            tibble::enframe() %>%
            dplyr::transmute(
                Function = p,
                Best = format_(base::min(value)),
                Worst = format_(base::max(value)),
                Median = format_(stats::median(value)),
                Mean = format_(base::mean(value)),
                Std = format_(stats::sd(value))
            ) %>%
            dplyr::slice(dplyr::n())
     })  %>%
     purrr::reduce(dplyr::bind_rows)
   table
}


get_cecResults = function(problems, dim, type = "m", cec = 21) {
interp_ = stringr::str_interp
  paths = c(
    interp_("~/cec2020/data/final/cec${cec}/ppmf-r-ave"),
    interp_("~/cec2020/data/final/cec${cec}/ppmf-s-ave"),
    interp_("~/cec2020/data/final/cec${cec}/ppmf-b-ave"),
    interp_("~/cec2020/data/final/cec${cec}/ppmf-sr-ave"),
    interp_("~/cec2020/data/final/cec${cec}/ppmf-bs-ave"),
    interp_("~/cec2020/data/final/cec${cec}/ppmf-br-ave"),
    interp_("~/cec2020/data/final/cec${cec}/ppmf-bsr-ave")
  )
  df = paths %>% purrr::map(function(path) {
    table = get_resultTable(
      idpath = path,
      problems = problems,
      dim = dim,
      type = type
    ) %>%
    dplyr::mutate(Path = path)
    sufix = stringr::str_extract(path, "ppmf-.*")
    table %>%
    readr::write_csv(
      file = interp_("~/cec2020/doc/csv/${sufix}-${dim}.csv")
    )
    })
  table
}
