################################################################################
##
##  Annual Ozone, Temperature, Mortality, Population data in Japan 
##  Chris Fook Sheng Ng, chrisng@m.u-tokyo.ac.jp
##  Matthew Lam
##
################################################################################

##### Ozone data #####
# Source: National Institude for Environmental Studies, Japan
library(readr)
library(dplyr)
library(stringr)
library(janitor)
library(rlang)

# JP -> EN map (same as before)
jp_to_en <- c(
  "測定年度" = "year_recorded",
  "都道府県名" = "prefecture",
  "市区町村コード" = "city_code",
  "市区町村名" = "city_name",
  "測定局名" = "station_name",
  "用途地域名" = "use_zone",
  "昼間の１時間値が0.06ppmを超えた日数(日)" = "days_over_0_06",
  "昼間の１時間値が0.06ppmを超えた時間数(時間)" = "hours_over_0_06",
  "昼間の１時間値の最高値(ppm)" = "max_hourly_ppm",
  "昼間の日最高１時間値の年平均値(ppm)" = "avg_daily_peak_ppm"
)

read_tokyo_file <- function(path) {
  readr::read_csv(
    file = path,
    locale = readr::locale(encoding = "CP932"),
    show_col_types = FALSE,
    guess_max = 10000,
    progress = FALSE
  )
}

year_from_filename <- function(path) {
  fn <- basename(path)
  as.integer(stringr::str_sub(fn, 3, 6))
}

summarise_wards <- function(path) {
  df_raw <- read_tokyo_file(path)
  
  keep_cols <- intersect(names(df_raw), names(jp_to_en))
  df <- df_raw |>
    select(all_of(keep_cols)) |>
    janitor::clean_names()
  
  # Map to English names
  present_map <- jp_to_en[names(jp_to_en) %in% names(df_raw)]
  names(present_map) <- janitor::make_clean_names(names(present_map))
  
  df <- df |>
    rename_with(~ dplyr::recode(.x, !!!present_map, .default = .x))
  
  # Ensure city_code is character
  df <- mutate(df, city_code = as.character(city_code))
  
  # Filter: 3rd char == "1" (wards)
  df <- filter(df, str_sub(city_code, 3, 3) == "1")
  
  # Filter by use_zone only if it exists AND has some "住"
  if ("use_zone" %in% names(df)) {
    if (any(!is.na(df$use_zone) & df$use_zone == "住")) {
      df <- filter(df, use_zone == "住")
    }
  }
  
  # If no rows remain, return NA metrics
  if (nrow(df) == 0L) {
    return(tibble::tibble(
      year = year_from_filename(path),
      mean_days_over_0_06       = NA_real_,
      mean_hours_over_0_06      = NA_real_,
      mean_max_hourly_ppm       = NA_real_,
      mean_avg_daily_peak_ppm   = NA_real_
    ))
  }
  
  summarise(
    df,
    year = year_from_filename(path),
    mean_days_over_0_06       = mean(days_over_0_06, na.rm = TRUE),
    mean_hours_over_0_06      = mean(hours_over_0_06, na.rm = TRUE),
    mean_max_hourly_ppm       = mean(max_hourly_ppm, na.rm = TRUE),
    mean_avg_daily_peak_ppm   = mean(avg_daily_peak_ppm, na.rm = TRUE),
    .groups = "drop"
  )
}

# ---- single file ----
# summarise_wards("Air Pollutant/Oxidant txt/TD19760613.txt")

# ---- apply to many files (1976..2022, etc.) ----
files <- list.files("Air Pollutant/Oxidant txt", 
                    pattern = "^TD\\d{8}\\.txt$", full.names = TRUE)[-c(1:5)]
results_ox <- purrr::map_dfr(files, summarise_wards) |> arrange(year)
results_ox

# lag yearly ozone by one year
results_ox_lagged <- results_ox %>%
  arrange(year) %>%
  mutate(
    mean_days_over_0_06_lag1  = lag(mean_days_over_0_06, 1),
    mean_hours_over_0_06_lag1 = lag(mean_hours_over_0_06, 1),
    mean_max_hourly_ppm_lag1  = lag(mean_max_hourly_ppm, 1),
    mean_avg_daily_peak_ppm_lag1 = lag(mean_avg_daily_peak_ppm, 1)
  )

glimpse(results_ox_lagged)

##### Temperature #####
library(readxl)
library(dplyr)
library(stringr)

yearly_temps <- read_excel("Weather/Weather_tokyo.xlsx", 
  skip = 6,
  .name_repair = "minimal") |>
  select(1, 8:12) |>
  setNames(c("year", "daily_mean", "daily_max", "daily_min",
             "max_temp", "min_temp")) |>
  mutate(
    year = as.integer(year),
    across(
      -year,
      ~ .x |>
        as.character() |>                        # treat as text
        str_remove_all("\\]") |>                 # remove trailing ]
        as.numeric() |>                          # convert to numeric
        round(1)                                 # 1 decimal place
    )
  ) |>
  filter(!is.na(year) & year <= 2024)

print(yearly_temps)

# lag yearly temperatures by one year
yearly_temps_lagged <- yearly_temps %>%
  arrange(year) %>%
  mutate(
    daily_mean_lag1 = lag(daily_mean, n = 1),
    daily_max_lag1 = lag(daily_max, n = 1),
    daily_min_lag1 = lag(daily_min, n = 1),
    max_temp_lag1 = lag(max_temp, n = 1),
    min_temp_lag1 = lag(min_temp, n = 1)
  )

print(yearly_temps_lagged)

##### Mortality & Population #####
library(readr)
library(stringr)
library(dplyr)
library(tidyr)

read_ipss_table <- function(path) {
  # Peek to find header / first data row
  ln <- readr::read_lines(path, n_max = 1000)
  hdr_line  <- which(str_detect(ln, "^\\s*Year\\s+Age\\s+Female\\s+Male\\s+Total\\s*$"))[1]
  data_line <- which(str_detect(ln, "^\\s*\\d{4}\\s+\\S+\\s+[-0-9\\.]+\\s+[-0-9\\.]+\\s+[-0-9\\.]+\\s*$"))[1]
  
  start <- if (!is.na(hdr_line)) hdr_line else data_line
  if (is.na(start)) stop("Could not detect the start of the data in: ", path)
  
  if (!is.na(hdr_line)) {
    skip_n <- hdr_line - 1L
    out <- readr::read_table(
      file = path,
      skip = skip_n,
      col_names = TRUE,
      col_types = readr::cols(
        Year   = readr::col_integer(),
        Age    = readr::col_character(),
        Female = readr::col_double(),
        Male   = readr::col_double(),
        Total  = readr::col_double()
      ),
      comment = "#"
    )
  } else {
    skip_n <- data_line - 1L
    out <- readr::read_table(
      file = path,
      skip = skip_n,
      col_names = c("Year", "Age", "Female", "Male", "Total"),
      col_types = "icddd",
      comment = "#"
    )
  }
  
  # Manually trim/squish whitespace (since no trim_ws arg)
  out |>
    mutate(
      Age = stringr::str_squish(Age),
      Age = stringr::str_replace_all(Age, "\\s+", "")
    )
}

collapse_age <- function(age) {
  age <- stringr::str_replace_all(age, "\\s+", "")
  dplyr::case_when(
    age == "0"                       ~ "0",
    age == "1-4"                     ~ "1-4",
    age %in% c("5-9","10-14")        ~ "5-14",
    age %in% c("15-19","20-24")      ~ "15-24",
    age %in% c("25-29","30-34","35-39","40-44") ~ "25-44",
    age %in% c("45-49","50-54","55-59","60-64") ~ "45-64",
    age %in% c("65-69","70-74")      ~ "65-74",
    age %in% c("75-79","80-84")      ~ "75-84",
    age %in% c("85-89","90-94")      ~ "85-94",
    age %in% c("95-99","100-104","105-109","110+") ~ "95+",
    TRUE ~ NA_character_
  )
}

aggregate_by_year_sex <- function(path, year_min = 1947, year_max = 2023) {
  dat <- read_ipss_table(path) |>
    filter(Year >= year_min, Year <= year_max) |>
    mutate(age_group = collapse_age(Age)) |>
    filter(!is.na(age_group)) |>
    select(Year, age_group, Female, Male)
  
  list(
    female = dat |>
      group_by(Year, age_group) |>
      summarise(value = sum(Female, na.rm = TRUE), .groups = "drop") |>
      pivot_wider(names_from = age_group, values_from = value) |>
      arrange(Year),
    male = dat |>
      group_by(Year, age_group) |>
      summarise(value = sum(Male, na.rm = TRUE), .groups = "drop") |>
      pivot_wider(names_from = age_group, values_from = value) |>
      arrange(Year)
  )
}

# 1) Deaths (counts) -> yearly time series (female & male separate)
deaths_ts <- aggregate_by_year_sex("Mortality/IPSS/Deaths_5x1.txt")

# 2) Population (counts) -> yearly time series (female & male separate)
pop_ts <- aggregate_by_year_sex("Mortality/IPSS/Population5.txt")

# peek
dplyr::glimpse(deaths_ts)
dplyr::glimpse(pop_ts)

# Combine data
# Define consistent ordering
age_levels <- c("0","1-4","5-14","15-24","25-44",
                "45-64","65-74","75-84","85-94","95+")
sex_levels <- c("Female", "Male")

# helper: convert wide df (e.g. deaths_ts$female) into long with sex/metric
make_long <- function(df, sex, metric) {
  df |>
    rename(year = Year) |>
    pivot_longer(
      cols = -year,
      names_to = "age_group",
      values_to = metric
    ) |>
    mutate(sex = sex)
}

combine_four <- function(deaths_ts, pop_ts) {
  # deaths
  d_f <- make_long(deaths_ts$female, "Female", "Deaths")
  d_m <- make_long(deaths_ts$male,   "Male",   "Deaths")
  # population
  p_f <- make_long(pop_ts$female,    "Female", "Population")
  p_m <- make_long(pop_ts$male,      "Male",   "Population")
  
  # join deaths + pop for each sex
  df_f <- left_join(d_f, p_f,
                    by = c("year", "age_group", "sex"))
  df_m <- left_join(d_m, p_m,
                    by = c("year", "age_group", "sex"))
  
  bind_rows(df_f, df_m) |>
    mutate(
      age_group = factor(age_group, levels = age_levels, ordered = TRUE),
      sex = factor(sex, levels = sex_levels)
    ) |>
    arrange(year, sex, age_group)
}

# ---- usage ----
long_ts <- combine_four(deaths_ts, pop_ts)
dplyr::glimpse(long_ts)

# Merge mortality, temperature, and ozone data
library(dplyr)

# long_ts       = mortality dataset (year, age_group, sex, Deaths, Population)
# yearly_temps  = temperature dataset (year, daily_mean, daily_max, daily_min, max_temp, min_temp)
# results_ox = ozone dataset (year, ozone variables...)

# Merge step by step
mort_env <- long_ts |>
  left_join(yearly_temps_lagged, by = "year") |>
  left_join(results_ox_lagged, by = "year") |>
  filter(year >= 1978, year <= 2022)  # remove 1976-1977 (start of Ox data)

# check
glimpse(mort_env)

# save(mort_env, file="mort_tky_year.RData")





##### Draft #####


##### Mortality by specific causes / diseases #####
library(readxl)
library(dplyr)
library(tidyr)

# Read Excel
df <- read_excel("Mortality/Hokeniryo/a0512.xlsx", sheet = 1, col_names = FALSE)

# Extract rows 43–45 (肺炎, 慢性閉塞性肺疾患, 喘息), cols 2–23
diseases_df <- df[43:45, 2:23] |>
  mutate(Disease = c("肺炎", "慢性閉塞性肺疾患", "喘息"), .before = 1)

# Assign Western year labels directly (2002–2023)
colnames(diseases_df) <- c("Disease", as.character(2002:2023))

# Replace Japanese disease names with English
disease_map <- c(
  "肺炎" = "Pneumonia",
  "慢性閉塞性肺疾患" = "COPD",
  "喘息" = "Asthma"
)

diseases_df$Disease <- recode(diseases_df$Disease, !!!disease_map)

# Reshape to long format
ts_long <- diseases_df |>
  pivot_longer(
    cols = -Disease,
    names_to = "year",
    values_to = "value"
  ) |>
  mutate(
    year = as.integer(year),
    value = as.numeric(value)
  ) |>
  arrange(Disease, year)

print(ts_long)
# the disease-specific mortality has no age group information
