
library(tidyverse)

common_lh <- 
  read_delim("commonx1a.out.out", delim = " ", col_names = FALSE) %>%
  select(-X11)

common_fsh <- 
  read_delim("commonx1a.out.out", delim = " ", col_names = FALSE) %>%
  select(-X11)

common_lh <- common_lh %>%
  gather(key = varname, value = value, -X1) %>%
  mutate_at(vars(varname), factor) %>%
  mutate_at(vars(value), as.numeric)

lh_trace <- common_lh %>%
  ggplot(aes(x = X1, y = value)) +
    geom_path() +
    facet_wrap(~ varname, scales = "free")

lh_posterior <- common_lh %>%
  ggplot(aes(x = value)) +
    geom_histogram() +
    facet_wrap(~ varname, scales = "free")



common_fsh <- common_fsh %>%
  gather(key = varname, value = value, -X1) %>%
  mutate_at(vars(varname), factor, ordered = TRUE) %>%
  mutate_at(vars(value), as.numeric)

fsh_trace <- common_fsh %>%
  ggplot(aes(x = X1, y = value)) +
    geom_path() +
    facet_wrap(~ varname, scales = "free")

fsh_posterior <- common_fsh %>%
  ggplot(aes(x = value)) +
    geom_histogram() +
    facet_wrap(~ varname, scales = "free")

