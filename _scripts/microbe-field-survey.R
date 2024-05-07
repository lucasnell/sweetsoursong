
#'
#' Plotting field survey data from Chappell et al. (2022; doi: 10.7554/eLife.79647)
#'
#' Download `2015_survey_samples_edit_cfu_master.csv` from here:
#'
#' https://archive.softwareheritage.org/browse/content/sha1_git:8dd550db96e00018a553fe2bfe0ca65ba5e42714/?origin_url=https://gitlab.com/teamnectarmicrobe/n06_nectarmicrobes_ecoevo&path=45281c8e25c1d9c6f9876b0680f42e884b606586/field_survey/2015_survey_samples_edit_cfu_master.csv&revision=023e6c17cafa08d701e56dbc9415dfdb75122a8c&snapshot=4d650c31be6067fb8741f028a60aaddde200e6fb
#'
#' and put it into the `_data` folder.
#' They're not my data to share, so I've left them out of this repo.
#'

suppressPackageStartupMessages({
    library(vegan)
    library(sweetsoursong)
    library(tidyverse)
    library(patchwork)
})







# (Starting with `read_lines` to remove empty columns 32:40 before reading.)

survey_df <- read_lines("_data/2015_survey_samples_edit_cfu_master.csv") |>
    str_split(",") |>
    map_chr(\(x) paste0(x[1:31], collapse = ",")) |>
    paste(collapse = "\n") |>
    read_csv(col_types = cols()) |>
    rename_with(tolower) |>
    rename(log1p_fcfu = `log(fcfu+1)`,
           log1p_bcfu = `log(bcfu+1)`,
           bcfu_ul = bcfuul) |>
    filter(! is.na(plant_id), !grepl("^75", plant_id)) |>
    mutate(site = str_sub(plant_id, 1, 2),
           plant = str_sub(plant_id, 3, 4),
           flower_id = paste(flower_id),
           log1p_fcfu = ifelse(str_detect(log1p_fcfu, "[:alpha:]"), "0",
                               log1p_fcfu),
           log1p_bcfu = ifelse(str_detect(log1p_bcfu, "[:alpha:]"), "0",
                               log1p_bcfu),
           fcfu_ul = ifelse(str_detect(fcfu_ul, "[:alpha:]"), "1", fcfu_ul),
           bcfu_ul = ifelse(str_detect(bcfu_ul, "[:alpha:]"), "0", bcfu_ul),
           log1p_fcfu = as.numeric(log1p_fcfu),
           log1p_bcfu = as.numeric(log1p_bcfu),
           fcfu_ul = as.numeric(fcfu_ul),
           bcfu_ul = as.numeric(bcfu_ul),
           fcfu_ul1 = fcfu_ul + 1,
           bcfu_ul1 = bcfu_ul + 1)




class_df <- survey_df |>
    split(~ site) |>
    map_dfr(\(dd) {
        comm_mat_i <- dd[,c("fcfu_ul1", "bcfu_ul1")] |>
            set_names(c("fungi", "bacteria")) |>
            as.matrix() |>
            t()
        suppressWarnings({clam_i <- clamtest(comm_mat_i, c("fungi", "bacteria"),
                                             coverage.limit = 10,
                                             specialization = 2/3,
                                             npoints = 5, alpha = 0.05)})
        dd <- dd |>
            select(site, plant, flower_id) |>
            mutate(class = fct_recode(clam_i$Classes,
                                      "Fungi-dominated" = "Specialist_fungi",
                                      "Co-dominated" = "Generalist",
                                      "Bacteria-dominated" = "Specialist_bacteria",
                                      "Uncolonized" = "Too_rare"))
        return(dd)
    })


# Summarize proportion dominated by yeast / bacteria by site and plant:
class_summ_df <- class_df |>
    group_by(site, plant) |>
    summarize(Y = mean(class == "Fungi-dominated"),
              B = mean(class == "Bacteria-dominated"),
              .groups = "drop") |>
    mutate(n_site = factor(site, levels = c("BS", "CH", "JP", "OH", "SA", "SV",
                                            "SG", "LH", "SB", "SR", "MW", "BB"),
                           labels = 1:12))


class_summ_df |>
    filter(site == "JP") |>
    ggplot(aes(Y * 100)) +
    geom_histogram(binwidth = 12.5, fill = "#FFCC33") +
    scale_x_continuous("Percent yeast-dominated", breaks = 0:3 * 20) +
    coord_cartesian(xlim = c(-6.25, 75)) +
    facet_wrap(~ n_site) +
    NULL


class_summ_df |>
    ggplot(aes(Y, B)) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray60") +
    geom_vline(xintercept = 0, linetype = 2, color = "gray60") +
    stat_bin_2d(aes(fill = after_stat(density)), binwidth = rep(0.125, 2)) +
    scale_x_continuous("Yeast proportion", limits = c(NA, 1)) +
    scale_y_continuous("Bacteria proportion", limits = c(NA, 1)) +
    scale_fill_viridis_c(option = "rocket", end = 0.9) +
    coord_equal()





