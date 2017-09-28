#' Identify Acute Kidney Injury
#'
#' Uses KDIGO Clinical practice guidline for acute kidney injury to
#' detect AKI events.
#'
#' Using Kidney Disease: Improving Global Outcomes (KDIGO)
#' guidelines, Acute Kidney Injury (AKI) are detected using this
#' function. Guidelines implemented can be found at:
#'
#' Selby N. M., Hill R., & Fluck R. J. (2015). Standardizing the Early
#' Identification of Acute Kidney Injury: The NHS England National
#' Patient Safety Alert. \emph{Nephron, 131, 113-117.}
#'
#' @param creatinine Data frame containing specific columns:
#'
#' \strong{ID:} Unique patient identification.
#'
#' \strong{COLLECTED_DTS:} Dates when corresponding samples were
#' collected. Must be R datetime object.
#'
#' \strong{RESULT:} Creatinine value cooresponding to collected sample.
#'
#' \strong{BIRTH_DTS:} Birth date of patient. Must be R date or
#' datetime object.
#'
#' \strong{RI_HIGH:} Reference index for high side of normal expected
#' creatinine levels.
#'
#' @return Returns a data frame identifying acute kidney injury with specific
#' colums:
#'
#' \strong{AKI_FLAG:} TRUE if some AKI event is detected, FALSE if
#' AKI event was not detected.
#'
#' \strong{STAGE:} Classifies the stage (1,2,3) of AKI event. Returns
#' NA if classification is unknown or AKI_FLAG is FALSE.
#'
#' \strong{BASELINE:} Baseline used in detection of AKI_FLAG. Returns
#' NA if no baseline was established.
#'
#' @keywords kdigo, kidney
#' @export
detect_aki <- function(creatinine) {
  creatinine <- tibble::as.tibble(creatinine)
  creatinine <- dplyr::arrange(creatinine, ID, COLLECTED_DTS)
  creatinine <- dplyr::mutate(creatinine,
                              C_DTS = as.numeric(readr::parse_datetime(COLLECTED_DTS)) / (24 * 3600),
                              B_DTS = as.numeric(readr::parse_datetime(BIRTH_DTS)) / (24 * 3600),
                              AKI_FLAG = as.logical(NA),
                              STAGE = as.numeric(NA),
                              BASELINE = as.numeric(NA))
  creatinine <- dplyr::group_by(creatinine, ID)
  # function global variable
  two_day_low <- NA

  # option to determine aki where there is no baseline
  no_baseline <- function(index) {
    if (creatinine$RESULT[[index]] < creatinine$RI_HIGH[[index]]) {
      creatinine$AKI_FLAG[[index]] <<- FALSE
    } else {
      creatinine$AKI_FLAG[[index]] <<- TRUE
    }
  }

  # determines whether there is a baseline and what that should be
  find_baseline <- function(previous, current) {
    differences <- creatinine$C_DTS[current] -1 * creatinine$C_DTS[previous]
    two_days <- creatinine$RESULT[previous[differences < 2]]
    week <- creatinine$RESULT[previous[differences <= 7]]
    year <- creatinine$RESULT[previous[differences <= 365]]
    two_day_low <<- NA
    if (length(two_days) > 0) {
      two_day_low <<- min(two_days)
    }
    if (length(week) > 0 && length(two_days) > 0) {
      return(min(min(week), median(year)))
    } else if (length(week) > 0) {
      return(min(week))
    } else if (length(year) > 0) {
      return(median(year))
    } else {
      return(NA)
    }
  }

  # use calculated baseline to identify aki
  use_baseline <- function(index, baseline) {
    cr <- creatinine$RESULT[[index]]
    ratio <- cr / baseline
    age <- creatinine$C_DTS[[index]] - creatinine$B_DTS[[index]]
    if (ratio >= 1.5) {
      creatinine$AKI_FLAG[[index]] <<- TRUE
      if ((age < 6570 && cr > 3 * creatinine$RI_HIGH[[index]]) || (age >= 6570 && cr > 4) || ratio > 3) {
        creatinine$STAGE[[index]] <<- 3
      } else if (ratio >= 2) {
        creatinine$STAGE[[index]] <<- 2
      } else {
        creatinine$STAGE[[index]] <<- 1
      }
    } else if (!is.na(two_day_low) && cr - two_day_low > .3) {
      creatinine$AKI_FLAG[[index]] <<- TRUE
      creatinine$STAGE[[index]] <<- 1
    } else {
      creatinine$AKI_FLAG[[index]] <<- FALSE
    }
  }

  groups <- attributes(creatinine)$indices
  for (i in seq_along(groups)) {
    for (j in seq_along(groups[[i]])) {
      index <- groups[[i]][[j]] + 1
      if (j == 1) {
        no_baseline(index)
      } else {
        baseline <- find_baseline(groups[[i]][1:(j - 1)] + 1, index)
        creatinine$BASELINE[[index]] <- baseline
        if (is.na(baseline)) {
          no_baseline(index)
        } else {
          use_baseline(index, baseline)
        }
      }
    }
  }
  return(dplyr::select(creatinine, -C_DTS, -B_DTS))
}


