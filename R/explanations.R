#' Check whether explanations are enabled
#'
#' Explanations are on by default. The setting is stored in the
#' `tidyseq.explain` option.
#'
#' @return TRUE if explanations are enabled globally, FALSE otherwise
#' @export
tidyseq_explanations_enabled <- function() {
  isTRUE(getOption("tidyseq.explain", default = TRUE))
}

#' Turn explanations on globally
#'
#' @return The previous setting, invisibly
#' @export
tidyseq_explanations_on <- function() {
  previous_setting <- tidyseq_explanations_enabled()
  options(tidyseq.explain = TRUE)
  message("TidySeq explanations are on.")
  invisible(previous_setting)
}

#' Turn explanations off globally
#'
#' @return The previous setting, invisibly
#' @export
tidyseq_explanations_off <- function() {
  previous_setting <- tidyseq_explanations_enabled()
  options(tidyseq.explain = FALSE)
  message("TidySeq explanations are off.")
  invisible(previous_setting)
}

#' Display an explanation of what a function does
#'
#' @param explanation_text Plain-language description of the step
#' @param example_code Example R code showing the underlying process
#' @param function_name Name of the calling function, used as a heading
#' @param explain TRUE or FALSE to override the global setting for this call;
#'   NULL to use the global setting
#'
#' @return NULL, invisibly
#' @keywords internal
display_explanation <- function(explanation_text, example_code, function_name, explain = NULL) {
  show_explanation <- if (is.null(explain)) tidyseq_explanations_enabled() else isTRUE(explain)
  if (!show_explanation) {
    return(invisible(NULL))
  }

  message(
    "\n== ", function_name, "() ==\n",
    explanation_text, "\n\n",
    "Example code:\n",
    example_code, "\n\n",
    "(Turn these explanations off with tidyseq_explanations_off() or explain = FALSE.)\n"
  )
  invisible(NULL)
}
