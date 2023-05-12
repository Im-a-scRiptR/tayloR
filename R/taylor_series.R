# Functions

#' @title Taylor Remainder
#' @description
#' The Taylor Remainder Function is integrated with the make_taylor_series
#' function to provide an error estimate when calculating the Taylor series
#' for a given function.
#'
#' @param df_np1 The n+1th derivative of the last derivative present in the
#' Taylor polynomial
#' @param x the value of x to plug into the absolute value of the
#' n+1th derivative
#' @param n The integer number of derivatives taken in total
#' @param a The center being tested
#' @returns Returns the absolute value of n+1th derivative multiplied by
#' x - a raised to a power of n+1 all divided by the factorial of n+1
#' @export
R_nx <-
  function(df_np1, x, n, a) {
    (abs(eval(df_np1, list(x = x))) / factorial(n)) * (abs(x - a) ^ (n))
  }

#' @title Taylor Series Converter
#' @description
#' The Taylor Series Converter Function performs a single variable Taylor
#' Expansion of the given function. It returns the call in full, the error
#' approximation, and the Taylor Remainder graph.
#'
#' @param f The function you want to convert into a Taylor Series
#' @param num_derivs How derivatives you want to take of f
#' @param center The center to test. In Stewart Calculus, this is notated "a"
#' @param intr_start The lower bound to try for the remainder interval
#' @param intr_end The upper bound to try for the remainder interval
#' @param var_name Differentiating f w.r.t var_name. Can only be one variable name
#' at the moment. The default variable name is "x".
#' @param open_closed Whether or not the interval defined by sequence
#' from intr_start to intr_end
#' @importFrom dplyr `%>%`
#' @export
make_taylor_series <-
  function(f,
           num_derivs,
           center,
           intr_start,
           intr_end,
           var_name = "x",
           open_closed = "closed") {

    # Redefine terms
    a <- center
    c <- intr_start
    d <- intr_end

    # Format the function for processing
    # If f is already a function, it will pull out the body else. Stop if the
    # supplied expression is not of class character or function
    if(is.function(f)) {
      if(length(body(f)) > 1) {
        f <- body(f)[[2L]]
      } else {
        f <- body(f)
      }
    } else if(is.character(f)) {
      f <- str2lang(f)
    } else {
      stop("Supplied expression must be of class character or function")
    }

    # Make an empty derivatives list
    d_list <- list()
    # Pre-fill the list
    d_list[[1]] <- f
    d_list[[2]] <- stats::D(d_list[[1]],var_name)

    # Finish up the symbolic successive differentiation(s)
    if(num_derivs >= 3) {
      for(i in 3:(num_derivs + 1)) {
        d_list[[i]] <- stats::D(d_list[[i - 1]],var_name)
      }
    }

    # Evaluate f at the center to check for NA, Inf, or NaN and make a mild correction
    test <- eval(f,list(x = a))
    center_check <- is.na(test) | is.infinite(test) | is.nan(test)
    if (center_check) {
      a <- a + 1e-6
    }

    # Calculate the terms for the Taylor series
    call_tbl <-
      tibble::tibble(n = 0:(length(d_list) - 1)) %>%
      dplyr::mutate(a = a, call = d_list) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(term_1 = eval(call, list(x = a)) / factorial(n)) %>%
      dplyr::select(-call) %>%
      dplyr::mutate(
        term_2 =
          glue::glue("(", "(x-a)^", "{n})") %>%
          stringr::str_replace("\\ba\\b", as.character(a)),
        T_n = paste0(as.character(term_1), " * ", term_2)
      ) %>%
      dplyr::select(-starts_with("term")) %>%
      dplyr::ungroup()

    # Combine terms into one Taylor Series
    T_n <- paste(call_tbl$T_n,collapse = " + ")

    # Make the n+1th derivative for the R_nx function
    df_np1 <- stats::D(d_list[[length(d_list)]], name = var_name)

    # Correct for interval openess if needed
    test_interval <- seq(c, d, 0.01)
    if (open_closed == "open") {
      test_interval <-
        test_interval[-which(as.character(test_interval) %in% as.character(c(c, d)))]
    }
    if (center_check) {
      test_interval <- test_interval[-which(test_interval == center)]
    }

    # Run the (n+1)th derivative's max value
    vals <-
      tibble::tibble(x = test_interval) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(f = abs(eval(df_np1, list(x = x)))) %>%
      dplyr::ungroup() %>%
      dplyr::slice_max(f) %>%
      dplyr::slice_max(x)

    # Plot the remainder function
    R_plot <-
      tibble::tibble(x = test_interval) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(R = R_nx(df_np1, x, length(d_list), a)) %>%
      dplyr::ungroup() %>%
      ggplot2::ggplot(ggplot2::aes(x, R)) + ggplot2::geom_line() +
      ggplot2::ggtitle(glue::glue("Remainder for f^{num_derivs}({var_name})"))

    # Calculate the error for the Taylor Expansion
    error <-
      (vals$f / factorial(length(d_list))) * abs(vals$x - a) ^ (length(d_list))

    # Return Taylor Series, its error, and the remainder plot
    list(
      T_n             = T_n,
      taylor_error    = error,
      remainder_graph = R_plot
    )
  }

# Testing

# Works with

# f <- \(x) {sin(x)}
# f <- \(x) sin(x)
# f <- "sin(x)"

# taylor <-
#   make_taylor_series(
#     f,
#     num_derivs   = 4,
#     center       = pi/6,
#     intr_start   = 0,
#     intr_end     = pi/3,
#     var_name     = "x",
#     open_closed  = "closed"
#   )

# taylor$taylor_error
